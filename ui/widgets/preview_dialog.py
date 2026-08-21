"""
Preview Dialogs - slider-based data browsing before Select Data.

ProfilePreviewDialog:      Browse time points, showing R-space profile at each time.
TimeTracePreviewDialog:    Browse channels, showing param1(t) and param2(t) at each R.
TranspProfilePreviewDialog:  Browse TRANSP CDF time points with variable selector.
TranspTimeTracePreviewDialog: Browse TRANSP CDF runs with variable selector.
"""

import numpy as np
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QSlider, QPushButton,
    QLineEdit, QComboBox, QStyle, QGroupBox, QWidget,
)
from PySide6.QtCore import Qt, QTimer

from ui.ui_constants import apply_dark_figure_style, get_icon
from ui.theme import ThemeManager

_ZERO_BOTTOM = {'Ti', 'Te', 'ne', 'q', 'j'}
_VT_YMIN = -30.0  # Fixed lower bound for vT [km/s] to avoid edge-noise outliers


def _robust_zmax(arr, edge_skip_frac=0.10, percentile=99.0):
    """Compute a robust upper bound for z-axis ignoring edge noise.

    arr: 2D ndarray with the radial axis as axis 0.
    edge_skip_frac: drop this fraction of points from the outer end of axis 0.
    percentile: take this percentile of remaining finite values.
    """
    a = np.asarray(arr, dtype=float)
    if a.ndim == 2 and a.shape[0] > 4:
        n_keep = max(1, int(round(a.shape[0] * (1.0 - edge_skip_frac))))
        a = a[:n_keep, :]
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return None
    return float(np.percentile(finite, percentile))


def _shrink_3d_box(ax, zoom=0.78):
    """Set 3D cube zoom relative to its axes bbox."""
    try:
        ax.set_box_aspect(None, zoom=zoom)        # matplotlib >= 3.6
    except TypeError:
        try:
            ax.dist = 10.0 / max(zoom, 0.1)        # legacy fallback
        except Exception:
            pass


def _layout_3d_axes(ax, col, n_cols, top_frac=0.92, bottom_frac=0.04,
                    shift_frac=0.30, shrink_frac=0.10):
    """Position a 3D axes so the cube anchors to the left of its panel
    region, leaving room on the right for the z-axis label.

    shift_frac: shift bbox left by this fraction of panel width.
    shrink_frac: reduce bbox width by this fraction so zlabel has clear
                 margin from the figure right edge.
    """
    panel_w = 1.0 / max(n_cols, 1)
    new_left = col * panel_w - panel_w * shift_frac
    new_w = panel_w * (1.0 - shrink_frac)
    height = max(0.05, top_frac - bottom_frac)
    ax.set_position([new_left, bottom_frac, new_w, height])


def _robust_zrange(arr, edge_skip_frac=0.10, low_pct=1.0, high_pct=99.0):
    """Compute robust (min, max) range for z-axis ignoring edge noise.

    Uses percentiles to suppress outliers (e.g. vT noise that produces
    extreme negative values at the plasma edge).
    """
    a = np.asarray(arr, dtype=float)
    if a.ndim == 2 and a.shape[0] > 4:
        n_keep = max(1, int(round(a.shape[0] * (1.0 - edge_skip_frac))))
        a = a[:n_keep, :]
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return None, None
    return (float(np.percentile(finite, low_pct)),
            float(np.percentile(finite, high_pct)))


# ============================================================
# Shared mixin for slider, playback, fix-axes, add-to-selected
# ============================================================

class _PreviewBase(QDialog):
    """Common infrastructure for preview dialogs."""

    def __init__(self, parent, n_items, title, selected_listbox=None):
        super().__init__(parent)
        self.selected_listbox = selected_listbox
        self.added_count = 0
        self._artists = []
        self._axes = []
        self._fix_axes = False
        self._n_items = max(1, n_items)
        self._view_mode = '2D'

        self._play_timer = QTimer()
        self._play_timer.timeout.connect(self._play_step)
        self._playing = False

        self.setWindowTitle(title)
        self.resize(1000, 550)

    # ---- Shared UI builder pieces ----

    def _build_left(self, main_layout, hint_text):
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)

        self.figure = Figure((8, 4.5))
        apply_dark_figure_style(self.figure)
        self.canvas = FigureCanvasQTAgg(self.figure)
        toolbar = NavigationToolbar2QT(self.canvas, self)
        left_layout.addWidget(toolbar)
        left_layout.addWidget(self.canvas, stretch=1)

        slider_row = QHBoxLayout()
        prev_btn = QPushButton()
        prev_btn.setIcon(get_icon(QStyle.SP_ArrowBack))
        prev_btn.setFixedWidth(28)
        prev_btn.clicked.connect(lambda: self._step(-1))
        slider_row.addWidget(prev_btn)

        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(self._n_items - 1)
        self.slider.valueChanged.connect(self._on_slider)
        slider_row.addWidget(self.slider, stretch=1)

        next_btn = QPushButton()
        next_btn.setIcon(get_icon(QStyle.SP_ArrowForward))
        next_btn.setFixedWidth(28)
        next_btn.clicked.connect(lambda: self._step(1))
        slider_row.addWidget(next_btn)
        left_layout.addLayout(slider_row)

        hint = QLabel(hint_text)
        hint.setStyleSheet("color: gray; font-size: 8pt;")
        hint.setAlignment(Qt.AlignCenter)
        left_layout.addWidget(hint)

        main_layout.addWidget(left, stretch=1)

    def _build_right(self, main_layout, nav_widgets_func):
        from ui.widgets.toggle_switch import ToggleSwitch

        ctrl = QWidget()
        ctrl.setFixedWidth(240)
        cl = QVBoxLayout(ctrl)
        cl.setContentsMargins(4, 0, 4, 4)

        # View mode (top): 2D / 3D
        from PySide6.QtWidgets import QRadioButton, QButtonGroup
        view_group = QGroupBox("View")
        view_layout = QHBoxLayout(view_group)
        self._view_btn_group = QButtonGroup(self)
        self.view_2d_radio = QRadioButton("2D")
        self.view_3d_radio = QRadioButton("3D")
        self.view_2d_radio.setChecked(True)
        self._view_btn_group.addButton(self.view_2d_radio)
        self._view_btn_group.addButton(self.view_3d_radio)
        view_layout.addWidget(self.view_2d_radio)
        view_layout.addWidget(self.view_3d_radio)
        view_layout.addStretch()
        self.view_2d_radio.toggled.connect(
            lambda c: c and self._on_view_changed("2D"))
        self.view_3d_radio.toggled.connect(
            lambda c: c and self._on_view_changed("3D"))
        cl.addWidget(view_group)

        # Navigation (subclass-specific)
        nav_group = QGroupBox("Navigation")
        nav_layout = QVBoxLayout(nav_group)
        nav_widgets_func(nav_layout)
        cl.addWidget(nav_group)

        # Playback
        play_group = QGroupBox("Playback")
        play_layout = QHBoxLayout(play_group)
        self.speed_combo = QComboBox()
        self.speed_combo.addItems(['1 fps', '2 fps', '5 fps', '10 fps', '20 fps', '50 fps'])
        self.speed_combo.setCurrentText('5 fps')
        play_layout.addWidget(self.speed_combo)
        self.play_btn = QPushButton()
        self.play_btn.setIcon(get_icon(QStyle.SP_MediaPlay))
        self.play_btn.clicked.connect(self._toggle_play)
        play_layout.addWidget(self.play_btn)
        stop_btn = QPushButton()
        stop_btn.setIcon(get_icon(QStyle.SP_MediaStop))
        stop_btn.clicked.connect(self._stop_play)
        play_layout.addWidget(stop_btn)
        cl.addWidget(play_group)

        # Dα signal panel (diagnostic Browse only — TRANSP/Derived opt out)
        if getattr(self, '_HAS_DALPHA', True):
            from PySide6.QtWidgets import QDoubleSpinBox
            da_group = QGroupBox("Dα")
            da_layout = QVBoxLayout(da_group)

            # Filter (TRANSP-style)
            filter_row = QHBoxLayout()
            filter_row.addWidget(QLabel("Filter"))
            self.dalpha_signal_filter = QLineEdit()
            self.dalpha_signal_filter.setPlaceholderText("Type to filter...")
            self.dalpha_signal_filter.textChanged.connect(
                self._populate_dalpha_combo)
            filter_row.addWidget(self.dalpha_signal_filter, 1)
            da_layout.addLayout(filter_row)

            # Signal selector (TRANSP-style: native popup + always-on scrollbar)
            self.dalpha_signal_combo = QComboBox()
            self.dalpha_signal_combo.setMaxVisibleItems(20)
            self.dalpha_signal_combo.setSizeAdjustPolicy(
                QComboBox.AdjustToMinimumContentsLengthWithIcon)
            self.dalpha_signal_combo.setMinimumContentsLength(10)
            self.dalpha_signal_combo.setStyleSheet("QComboBox { combobox-popup: 0; }")
            self.dalpha_signal_combo.view().setVerticalScrollBarPolicy(
                Qt.ScrollBarAlwaysOn)
            self.dalpha_signal_combo.currentTextChanged.connect(
                self._on_dalpha_signal_changed)
            self._populate_dalpha_combo()
            da_layout.addWidget(self.dalpha_signal_combo)

            # Raw data toggle: signal → -<NODE>:FOO
            raw_row = QHBoxLayout()
            self.dalpha_raw_toggle = ToggleSwitch()
            self.dalpha_raw_toggle.toggled.connect(self._on_dalpha_raw_changed)
            raw_row.addWidget(self.dalpha_raw_toggle)
            raw_row.addWidget(QLabel("raw data (:FOO)"))
            raw_row.addStretch()
            da_layout.addLayout(raw_row)

            # Zoom
            dz_row = QHBoxLayout()
            self.dalpha_zoom_toggle = ToggleSwitch()
            self.dalpha_zoom_toggle.toggled.connect(self._on_dalpha_zoom)
            dz_row.addWidget(self.dalpha_zoom_toggle)
            dz_row.addWidget(QLabel("Zoom  ±"))
            self.dalpha_zoom_spin = QDoubleSpinBox()
            self.dalpha_zoom_spin.setRange(0.05, 5.0)
            self.dalpha_zoom_spin.setDecimals(2)
            self.dalpha_zoom_spin.setSingleStep(0.05)
            self.dalpha_zoom_spin.setSuffix(" s")
            self.dalpha_zoom_spin.setFixedWidth(80)
            self.dalpha_zoom_spin.valueChanged.connect(self._on_dalpha_zoom_changed)
            dz_row.addWidget(self.dalpha_zoom_spin)
            dz_row.addStretch()
            da_layout.addLayout(dz_row)
            cl.addWidget(da_group)

        # Options (Fix Axes — generic, kept at the bottom)
        opt_group = QGroupBox("Options")
        opt_layout = QHBoxLayout(opt_group)
        self.fix_axes_toggle = ToggleSwitch()
        self.fix_axes_toggle.toggled.connect(self._on_fix_axes)
        opt_layout.addWidget(self.fix_axes_toggle)
        opt_layout.addWidget(QLabel("Fix Axes"))
        opt_layout.addStretch()
        cl.addWidget(opt_group)

        cl.addStretch()

        self.add_btn = QPushButton("Add to Selected")
        self.add_btn.clicked.connect(self._add_current)
        cl.addWidget(self.add_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self._on_close)
        cl.addWidget(close_btn)

        main_layout.addWidget(ctrl)

    # ---- Navigation / Playback (shared) ----

    def _on_slider(self, value):
        if self._view_mode == '3D':
            self._render_3d_highlight(value)
        else:
            self._update_plot(value)

    def _step(self, delta):
        self.slider.setValue(
            max(0, min(self.slider.maximum(), self.slider.value() + delta)))

    def _on_fix_axes(self, checked):
        self._fix_axes = checked
        if not checked:
            self._update_plot(self.slider.value())

    def _on_dalpha_zoom(self, checked):
        self._dalpha_zoom = bool(checked)
        ax = getattr(self, '_ax_dalpha', None)
        if ax is not None and not checked:
            full = getattr(self, '_dalpha_full_xlim', None)
            if full is not None:
                ax.set_xlim(*full)
        if self._view_mode == '2D':
            self._update_plot(self.slider.value())

    def _on_dalpha_zoom_changed(self, value):
        self._dalpha_zoom_half = float(value)
        if self._dalpha_zoom and self._view_mode == '2D':
            self._update_plot(self.slider.value())

    _DALPHA_CATEGORIES = (
        [f"\\TOR_HA{n:02d}" for n in range(1, 21)],
        [f"\\POL_HA{n:02d}" for n in range(1, 62)],
        [f"\\DIV_KHA{n:02d}" for n in range(1, 10)],
        [f"\\DIV_GHA{n:02d}" for n in range(1, 10)],
    )

    def _populate_dalpha_combo(self):
        """Refill combo based on filter text. Keep category separators only
        when no filter is active (otherwise filtered list is flat)."""
        filt = self.dalpha_signal_filter.text().strip().upper()
        current = self.dalpha_signal_combo.currentText() or self._dalpha_signal

        self.dalpha_signal_combo.blockSignals(True)
        self.dalpha_signal_combo.clear()
        first = True
        for cat in self._DALPHA_CATEGORIES:
            items = [s for s in cat if not filt or filt in s.upper()]
            if not items:
                continue
            if not filt and not first:
                self.dalpha_signal_combo.insertSeparator(
                    self.dalpha_signal_combo.count())
            for s in items:
                self.dalpha_signal_combo.addItem(s)
            first = False

        idx = self.dalpha_signal_combo.findText(current)
        if idx >= 0:
            self.dalpha_signal_combo.setCurrentIndex(idx)
        self.dalpha_signal_combo.blockSignals(False)

    def _on_dalpha_signal_changed(self, text):
        if not text or text.startswith('---'):
            return
        if text == self._dalpha_signal:
            return
        self._dalpha_signal = text
        self._reload_dalpha()

    def _on_dalpha_raw_changed(self, checked):
        self._dalpha_raw = bool(checked)
        self._reload_dalpha()

    def _reload_dalpha(self):
        self._dalpha_t, self._dalpha_v = self._load_dalpha(
            self._dalpha_signal, self._dalpha_raw)
        if self._view_mode == '2D':
            self._init_plot()
            self._update_plot(self.slider.value())

    def _get_fps(self):
        return int(self.speed_combo.currentText().split()[0])

    def _toggle_play(self):
        if self._playing:
            self._stop_play()
        else:
            self._playing = True
            self.play_btn.setIcon(get_icon(QStyle.SP_MediaPause))
            self._play_timer.start(max(1, 1000 // self._get_fps()))

    def _stop_play(self):
        self._playing = False
        self.play_btn.setIcon(get_icon(QStyle.SP_MediaPlay))
        self._play_timer.stop()

    def _play_step(self):
        if self.slider.value() >= self.slider.maximum():
            self._stop_play()
            return
        self.slider.setValue(self.slider.value() + 1)
        self._play_timer.setInterval(max(1, 1000 // self._get_fps()))

    def _on_close(self):
        self._stop_play()
        self.accept()

    # ---- View mode (2D / 3D) ----

    def _on_view_changed(self, mode):
        self._view_mode = mode
        # Fix Axes only meaningful in 2D
        if hasattr(self, 'fix_axes_toggle'):
            self.fix_axes_toggle.setEnabled(mode == '2D')

        if mode == '3D':
            self._render_3d()
            # Show highlight at current slider position
            if hasattr(self, 'slider'):
                self._render_3d_highlight(self.slider.value())
        else:
            self._init_plot()
            self._update_plot(self.slider.value() if hasattr(self, 'slider') else 0)

    def _render_3d(self):
        """Render full dataset as 3D surface. Override in subclasses."""
        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)
        ax.text(0.5, 0.5, "3D view not available for this dialog",
                ha='center', va='center', fontsize=14, color='gray',
                transform=ax.transAxes)
        ax.axis('off')
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()

    def _render_3d_highlight(self, idx):
        """Draw a highlight on the 3D surface at slider index. Override."""
        pass

    # ---- Abstract (override in subclasses) ----

    def _update_plot(self, idx):
        raise NotImplementedError

    def _add_current(self):
        raise NotImplementedError


# ============================================================
# ProfilePreviewDialog
# ============================================================

class ProfilePreviewDialog(_PreviewBase):
    """Non-modal dialog for browsing profile data (slider = time)."""

    def __init__(self, parent, data, shot_number, source,
                 param1_name, param2_name, param1_label, param2_label,
                 selected_listbox=None, extra_data=None):
        self.data = data
        self.extra_data = extra_data
        self.shot_number = shot_number
        self.source = source
        self.param1_name = param1_name
        self.param2_name = param2_name
        self.param1_label = param1_label
        self.param2_label = param2_label
        self._is_ece_only = (source == 'ECE')
        self.time_arr = np.asarray(getattr(data, 'time_prof', data.time), dtype=float)
        # D-alpha context (signal selectable from the Browse dialog)
        self._dalpha_signal = '\\TOR_HA11'
        self._dalpha_raw = False
        self._dalpha_t, self._dalpha_v = self._load_dalpha(
            self._dalpha_signal, self._dalpha_raw)

        super().__init__(parent, len(self.time_arr),
                         f"Browse  #{shot_number} ({source})",
                         selected_listbox)
        self._build_ui()
        self._init_plot()
        self._update_plot(0)

    # --- D-alpha helpers ---

    @staticmethod
    def _load_ip_window(mds, shot, threshold_a=1.0e5):
        """Return (t_start, t_end) of |Ip| > threshold, or (None, None)."""
        try:
            ip = np.asarray(mds.get('\\PCRC03').data(), dtype=np.float32)
            ip_t = np.asarray(mds.get('dim_of(\\PCRC03)').data(), dtype=np.float32)
            m = np.abs(ip) > threshold_a
            if np.any(m):
                return float(ip_t[m][0]), float(ip_t[m][-1])
        except Exception:
            pass
        return None, None

    def _load_dalpha(self, node='\\TOR_HA11', raw=False):
        """Load the selected Dα-like signal + \\PCRC03 (Ip window).

        node : MDS+ tag (e.g. '\\TOR_HA11', '\\POL_HA01', '\\DIV_KHA05').
        raw  : if True, fetch '-<node>:FOO' (sign-flipped raw photodiode
               counts); otherwise the calibrated <node>.
        """
        tag = f"{node}:FOO" if raw else node
        value_expr = f"-{tag}" if raw else tag
        try:
            from MDSplus import Connection
            from config.app_config import AppConfig
            mds = Connection(AppConfig().MDS_IP)
            mds.openTree('kstar', int(self.shot_number))
            try:
                v = np.asarray(mds.get(value_expr).data(), dtype=np.float32)
                t = np.asarray(mds.get(f'dim_of({tag})').data(), dtype=np.float32)
                ip_start, ip_end = self._load_ip_window(mds, int(self.shot_number))
            finally:
                mds.closeTree('kstar', int(self.shot_number))
            if t.size == 0 or v.size == 0:
                self._dalpha_window = (None, None)
                return None, None
            self._dalpha_window = (ip_start, ip_end)
            wstr = (f" (Ip window {ip_start:.2f}-{ip_end:.2f}s)"
                    if ip_start is not None else "")
            print(f"[Dalpha] #{self.shot_number} {value_expr} loaded{wstr}")
            return t, v
        except Exception as e:
            print(f"[Dalpha] not available ({value_expr}) for #{self.shot_number}: {e}")
            self._dalpha_window = (None, None)
            return None, None

    def _build_ui(self):
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(4)

        self._build_left(main_layout, "Use mouse wheel or arrow keys to navigate frames.")

        def _nav(layout):
            row1 = QHBoxLayout()
            row1.addWidget(QLabel("Frame"))
            self.frame_entry = QLineEdit("1")
            self.frame_entry.returnPressed.connect(self._goto_frame)
            row1.addWidget(self.frame_entry)
            row1.addWidget(QLabel(f"/ {len(self.time_arr)}"))
            go = QPushButton()
            go.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
            go.setFixedSize(24, 24)
            go.clicked.connect(self._goto_frame)
            row1.addWidget(go)
            layout.addLayout(row1)

            row2 = QHBoxLayout()
            row2.addWidget(QLabel("Time [s]"))
            self.time_entry = QLineEdit("0.0")
            self.time_entry.returnPressed.connect(self._goto_time)
            row2.addWidget(self.time_entry)
            go2 = QPushButton()
            go2.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
            go2.setFixedSize(24, 24)
            go2.clicked.connect(self._goto_time)
            row2.addWidget(go2)
            layout.addLayout(row2)

        self._build_right(main_layout, _nav)

    def _goto_frame(self):
        try:
            idx = max(0, min(len(self.time_arr) - 1, int(self.frame_entry.text()) - 1))
        except ValueError:
            return
        # Refresh even when the target equals the current frame (setValue would
        # not emit valueChanged, so the plot/entries would not update otherwise)
        if self.slider.value() == idx:
            self._update_plot(idx)
        else:
            self.slider.setValue(idx)

    def _goto_time(self):
        try:
            idx = int(np.argmin(np.abs(self.time_arr - float(self.time_entry.text()))))
        except (ValueError, TypeError):
            return
        if self.slider.value() == idx:
            self._update_plot(idx)
        else:
            self.slider.setValue(idx)

    # ---- Plot ----

    def _init_plot(self):
        self.figure.clear()
        self._axes = []
        self._ax_dalpha = None
        self._dalpha_vline = None
        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'

        has_dalpha = self._dalpha_t is not None
        if has_dalpha:
            gs = self.figure.add_gridspec(
                2, 2, height_ratios=[3, 1], hspace=0.35)
        else:
            gs = self.figure.add_gridspec(1, 2)

        for col, (pname, plabel) in enumerate([
            (self.param1_name, self.param1_label),
            (self.param2_name, self.param2_label),
        ]):
            ax = self.figure.add_subplot(gs[0, col])
            ax.set_xlabel('R [m]', fontsize=12)
            ax.set_ylabel(plabel, fontsize=12)
            ax.tick_params(labelsize=10)
            ax.grid(ls='--', lw=0.3, color='#444444')
            if 'v' in pname.lower():
                ax.axhline(y=0, color=zc, ls='--', gid='zero_ref')
            self._axes.append(ax)

        if has_dalpha:
            ax_da = self.figure.add_subplot(gs[1, :])
            ax_da.plot(self._dalpha_t, self._dalpha_v,
                       '-', color='#888', lw=0.8)
            ax_da.set_xlabel('Time [s]', fontsize=10)
            ylabel = self._dalpha_signal
            if self._dalpha_raw:
                ylabel = f"-{ylabel}:FOO"
            ax_da.set_ylabel(ylabel, fontsize=10)
            ax_da.tick_params(labelsize=9)
            ax_da.grid(ls='--', lw=0.3, color='#444444')
            # Default x-range = Ip-active window if available
            ip_start, ip_end = getattr(self, '_dalpha_window', (None, None))
            if ip_start is not None and ip_end is not None:
                ax_da.set_xlim(ip_start, ip_end)
            else:
                ax_da.set_xlim(float(self.time_arr[0]),
                               float(self.time_arr[-1]))
            self._ax_dalpha = ax_da
            self._dalpha_full_xlim = ax_da.get_xlim()
            # Adaptive zoom window: 5% of Ip-active duration, clamped to a
            # reasonable absolute range so neither very short nor very long
            # shots produce a useless window.
            if ip_start is not None and ip_end is not None:
                duration = max(0.0, float(ip_end - ip_start))
                self._dalpha_zoom_half = float(np.clip(0.05 * duration, 0.1, 1.0))
            else:
                self._dalpha_zoom_half = 0.3
            self._dalpha_zoom = True
            try:
                self.dalpha_zoom_toggle.setChecked(True, animate=False)
            except Exception:
                pass
            # Sync spinbox to the adaptive default (without re-triggering
            # the valueChanged signal during initialization).
            if hasattr(self, 'dalpha_zoom_spin'):
                self.dalpha_zoom_spin.blockSignals(True)
                self.dalpha_zoom_spin.setValue(self._dalpha_zoom_half)
                self.dalpha_zoom_spin.blockSignals(False)

        self.figure.subplots_adjust(
            left=0.10, right=0.97, top=0.92, bottom=0.10, wspace=0.30)
        apply_dark_figure_style(self.figure)

    def _get_param_data(self, pname, t_actual):
        data = self.data
        if pname not in data.measurements:
            return None, None, None, None
        pdata, eu, el = data.get_parameter_asymmetric(pname)
        if hasattr(data, 'time_prof') and pname in ('q', 'j'):
            t_idx = np.argmin(np.abs(data.time_prof - t_actual))
            meas = data.measurements[pname]
            R = meas['roa'][:, t_idx] if 'roa' in meas else data.radius[:pdata.shape[0]]
        else:
            t_idx = np.argmin(np.abs(data.time - t_actual))
            R = data.radius[:pdata.shape[0]]
        return R, pdata[:, t_idx], eu[:, t_idx], el[:, t_idx]

    def _update_plot(self, t_idx):
        t_actual = self.time_arr[t_idx]
        ms = int(round(t_actual * 1000))

        self.figure.suptitle(
            f"#{self.shot_number}  {ms:06d}ms  ({self.source})", fontsize=12)
        self.frame_entry.setText(str(t_idx + 1))
        self.time_entry.setText(f"{t_actual:.4f}")

        # Sync D-alpha vertical line to the current time
        if self._ax_dalpha is not None:
            if self._dalpha_vline is not None:
                try:
                    self._dalpha_vline.remove()
                except Exception:
                    pass
            self._dalpha_vline = self._ax_dalpha.axvline(
                t_actual, color='red', lw=1.5, zorder=20)
            # In zoom mode, slide the x-axis to follow current time
            if getattr(self, '_dalpha_zoom', False):
                half = getattr(self, '_dalpha_zoom_half', 0.5)
                self._ax_dalpha.set_xlim(t_actual - half, t_actual + half)

        for a in self._artists:
            a.remove()
        self._artists = []

        if self._is_ece_only:
            marker, mfc, mew = 's', 'none', 1.5
        else:
            marker, mfc, mew = 'o', None, 1

        for col, pname in enumerate([self.param1_name, self.param2_name]):
            ax = self._axes[col]
            R, profile, eu, el = self._get_param_data(pname, t_actual)
            if R is None:
                continue
            valid = np.isfinite(profile)
            if not np.any(valid):
                continue

            kwargs = dict(fmt=marker, markersize=5, capsize=5,
                          color='#1f77b4', markeredgewidth=mew, zorder=10)
            if mfc is not None:
                kwargs['markerfacecolor'] = mfc
            c = ax.errorbar(R[valid], profile[valid],
                            yerr=[el[valid], eu[valid]], **kwargs)
            self._artists.append(c)

            if self.extra_data is not None and pname == 'Te':
                self._overlay_extra(ax, pname, t_actual)

            if not self._fix_axes:
                y_vals = profile[valid]
                y_min = np.nanmin(y_vals - el[valid])
                y_max = np.nanmax(y_vals + eu[valid])
                margin = (y_max - y_min) * 0.1 if y_max > y_min else 0.1
                if pname in _ZERO_BOTTOM:
                    ax.set_ylim(0, y_max + margin)
                elif pname == 'vT':
                    ax.set_ylim(_VT_YMIN, y_max + margin)
                else:
                    ax.set_ylim(y_min - margin, y_max + margin)

        self.canvas.draw_idle()

    def _render_3d(self):
        """Surface plot of param1 and param2 vs (R, time)."""
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        self.figure.clear()
        self._axes = []
        self._3d_highlight = []   # store highlight artists per axis
        self._3d_data = []        # per-axis (R, time, Z, pname) for highlight redraw
        time = np.asarray(self.time_arr)
        R = np.asarray(self.data.radius)
        for col, (pname, plabel) in enumerate([
            (self.param1_name, self.param1_label),
            (self.param2_name, self.param2_label),
        ]):
            ax = self.figure.add_subplot(1, 2, col + 1, projection='3d')
            self._axes.append(ax)
            if pname not in self.data.measurements:
                ax.text2D(0.5, 0.5, f"{pname} not available",
                          ha='center', va='center', transform=ax.transAxes,
                          color='gray')
                self._3d_data.append(None)
                self._3d_highlight.append(None)
                continue
            try:
                pdata, _, _ = self.data.get_parameter_asymmetric(pname)
            except Exception:
                self._3d_data.append(None)
                self._3d_highlight.append(None)
                continue
            n_R = min(pdata.shape[0], len(R))
            n_t = min(pdata.shape[1], len(time))
            if n_R == 0 or n_t == 0:
                self._3d_data.append(None)
                self._3d_highlight.append(None)
                continue
            Rv = R[:n_R]
            Tv = time[:n_t]
            Z = pdata[:n_R, :n_t]
            R_grid, T_grid = np.meshgrid(Rv, Tv, indexing='ij')
            try:
                ax.plot_surface(R_grid, T_grid, Z, cmap='viridis',
                                edgecolor='none', alpha=0.6)
            except Exception:
                ax.plot_wireframe(R_grid, T_grid, Z, color='#1f77b4', lw=0.5)
            ax.set_xlabel('R [m]', fontsize=10, labelpad=2)
            ax.set_ylabel('Time [s]', fontsize=10, labelpad=2)
            ax.set_zlabel(plabel, fontsize=10, labelpad=2)
            ax.tick_params(labelsize=8, pad=2)
            _shrink_3d_box(ax)
            _layout_3d_axes(ax, col, 2)
            # Robust z-axis bounds (suppress edge-noise outliers)
            z_lo, z_hi = _robust_zrange(Z)
            if z_hi is not None and np.isfinite(z_hi):
                if pname in _ZERO_BOTTOM:
                    z_lo_use = 0.0
                elif pname == 'vT':
                    z_lo_use = _VT_YMIN
                else:
                    margin = (z_hi - z_lo) * 0.05 if z_lo is not None else 0
                    z_lo_use = (z_lo - margin) if z_lo is not None else 0
                ax.set_zlim(z_lo_use, z_hi * 1.05 if z_hi > 0 else z_hi * 0.95)
            self._3d_data.append((Rv, Tv, Z, pname))
            self._3d_highlight.append(None)
        self.figure.suptitle(
            f"#{self.shot_number}  ({self.source})  3D view", fontsize=12)
        # Wider margins so 3D axis labels/ticks don't clip
        self.figure.subplots_adjust(left=0.0, right=1.0, top=0.94, bottom=0.0,
                                    wspace=0.05)
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()

    def _render_3d_highlight(self, t_idx):
        """Draw a profile slice at the current time index on each 3D axis."""
        if not getattr(self, '_3d_data', None):
            return
        time = np.asarray(self.time_arr)
        if t_idx >= len(time):
            return
        t_actual = float(time[t_idx])
        ms = int(round(t_actual * 1000))
        # Sync navigation widgets
        if hasattr(self, 'frame_entry'):
            self.frame_entry.setText(str(t_idx + 1))
        if hasattr(self, 'time_entry'):
            self.time_entry.setText(f"{t_actual:.4f}")
        for col, dataset in enumerate(self._3d_data):
            ax = self._axes[col]
            old = self._3d_highlight[col]
            if old is not None:
                try:
                    old.remove()
                except Exception:
                    pass
            if dataset is None:
                continue
            Rv, Tv, Z, _ = dataset
            ti = min(t_idx, Z.shape[1] - 1)
            t_val = float(Tv[ti])
            line, = ax.plot(Rv, np.full_like(Rv, t_val), Z[:, ti],
                            color='red', lw=2.5, zorder=20)
            self._3d_highlight[col] = line
        self.figure.suptitle(
            f"#{self.shot_number}  {ms:06d}ms  ({self.source})  3D view",
            fontsize=12)
        self.canvas.draw_idle()

    def _overlay_extra(self, ax, pname, t_actual):
        ed = self.extra_data
        if pname not in ed.measurements:
            return
        t_idx = np.argmin(np.abs(ed.time - t_actual))
        te_data = ed.measurements[pname]['data'][:, t_idx]
        R_ece = ed.radius
        valid = np.isfinite(te_data) & (te_data > 0)
        if np.any(valid):
            line, = ax.plot(R_ece[valid], te_data[valid], 's',
                            color='#1f77b4', markersize=5,
                            markerfacecolor='none', markeredgewidth=1.5,
                            zorder=5)
            self._artists.append(line)

    # ---- Add to Selected ----

    def _add_current(self):
        t_val = self.time_arr[self.slider.value()]
        ms = int(round(t_val * 1000))
        entry = f'{self.shot_number:06d}_{t_val*1e3:06.0f} ({self.source})'

        if self.selected_listbox is not None:
            existing = {self.selected_listbox.item(i).text()
                        for i in range(self.selected_listbox.count())}
            if entry in existing:
                self.add_btn.setText(f"{ms:06d}ms already selected")
                return
            self.selected_listbox.addItem(entry)

        self.added_count += 1
        self.add_btn.setText(f"Added {ms:06d}ms  ({self.added_count} total)")


# ============================================================
# TimeTracePreviewDialog
# ============================================================

class TimeTracePreviewDialog(_PreviewBase):
    """Non-modal dialog for browsing time trace channels (slider = channel).

    channels: list of dicts with keys:
        'entry'  - listbox entry string
        'label'  - display label (e.g. 'mod01 R=1820mm')
        'time'   - 1D time array [s]
        'param1' - 1D data array for param1
        'param2' - 1D data array for param2 (or None)
    """

    _HAS_DALPHA = False  # traces are already vs time; no Dα subplot here

    def __init__(self, parent, shot_number, source,
                 param1_label, param2_label,
                 channels, selected_listbox=None):
        self.shot_number = shot_number
        self.source = source
        self.param1_label = param1_label
        self.param2_label = param2_label
        self.channels = channels

        super().__init__(parent, len(channels),
                         f"Browse  #{shot_number} ({source})",
                         selected_listbox)
        self._build_ui()
        self._init_plot()
        if self.channels:
            self._update_plot(0)

    def _build_ui(self):
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(4)

        self._build_left(main_layout, "Use mouse wheel or arrow keys to navigate channels.")

        def _nav(layout):
            row = QHBoxLayout()
            row.addWidget(QLabel("Channel"))
            self.idx_entry = QLineEdit("0")
            self.idx_entry.returnPressed.connect(self._goto_idx)
            row.addWidget(self.idx_entry)
            row.addWidget(QLabel(f"/ {self._n_items - 1}"))
            go = QPushButton()
            go.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
            go.setFixedSize(24, 24)
            go.clicked.connect(self._goto_idx)
            row.addWidget(go)
            layout.addLayout(row)

            self.ch_info_label = QLabel("")
            self.ch_info_label.setStyleSheet("color: #aaa; font-size: 9pt;")
            layout.addWidget(self.ch_info_label)

        self._build_right(main_layout, _nav)

    def _goto_idx(self):
        try:
            self.slider.setValue(
                max(0, min(len(self.channels) - 1, int(self.idx_entry.text()))))
        except ValueError:
            pass

    # ---- Plot ----

    def _init_plot(self):
        self.figure.clear()
        self._axes = []
        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'

        ax1 = self.figure.add_subplot(2, 1, 1)
        ax1.set_ylabel(self.param1_label, fontsize=12)
        ax1.tick_params(labelsize=10)
        ax1.grid(ls='--', lw=0.3, color='#444444')
        self._axes.append(ax1)

        ax2 = self.figure.add_subplot(2, 1, 2, sharex=ax1)
        ax2.set_xlabel('Time [s]', fontsize=12)
        ax2.set_ylabel(self.param2_label, fontsize=12)
        ax2.tick_params(labelsize=10)
        ax2.grid(ls='--', lw=0.3, color='#444444')
        if 'v' in self.param2_label.lower():
            ax2.axhline(y=0, color=zc, ls='--', lw=0.8, alpha=0.5)
        self._axes.append(ax2)

        self.figure.subplots_adjust(
            left=0.12, right=0.97, top=0.90, bottom=0.14, hspace=0.15)
        apply_dark_figure_style(self.figure)

    def _update_plot(self, ch_idx):
        if ch_idx >= len(self.channels):
            return
        ch = self.channels[ch_idx]

        self.figure.suptitle(
            f"#{self.shot_number}  {ch['label']}  ({self.source})", fontsize=12)
        self.idx_entry.setText(str(ch_idx))
        self.ch_info_label.setText(ch['label'])

        for a in self._artists:
            a.remove()
        self._artists = []

        time = ch['time']
        color = '#1f77b4'

        # Param1
        ax1 = self._axes[0]
        p1 = ch['param1']
        if p1 is not None:
            valid = np.isfinite(p1)
            if np.any(valid):
                line, = ax1.plot(time[valid], p1[valid], '-', color=color, lw=1)
                self._artists.append(line)
                if not self._fix_axes:
                    y = p1[valid]
                    margin = (np.nanmax(y) - np.nanmin(y)) * 0.1 or 0.1
                    param_key = ''.join(c for c in self.param1_label if c.isalpha())
                    if any(k in param_key for k in ('Ti', 'Te', 'ne')):
                        ax1.set_ylim(0, np.nanmax(y) + margin)
                    elif 'vT' in param_key or 'vt' in param_key:
                        ax1.set_ylim(_VT_YMIN, np.nanmax(y) + margin)
                    else:
                        ax1.set_ylim(np.nanmin(y) - margin, np.nanmax(y) + margin)

        # Param2
        ax2 = self._axes[1]
        p2 = ch['param2']
        if p2 is not None:
            valid = np.isfinite(p2)
            if np.any(valid):
                line, = ax2.plot(time[valid], p2[valid], '-', color=color, lw=1)
                self._artists.append(line)
                if not self._fix_axes:
                    y = p2[valid]
                    margin = (np.nanmax(y) - np.nanmin(y)) * 0.1 or 0.1
                    param_key2 = ''.join(c for c in self.param2_label if c.isalpha())
                    if 'vT' in param_key2 or 'vt' in param_key2:
                        ax2.set_ylim(_VT_YMIN, np.nanmax(y) + margin)
                    else:
                        ax2.set_ylim(np.nanmin(y) - margin, np.nanmax(y) + margin)

        if not self._fix_axes:
            ax1.set_xlim(time[0], time[-1])

        self.canvas.draw_idle()

    def _render_3d(self):
        """Surface plot: param vs (time, channel index/R)."""
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        self.figure.clear()
        self._axes = []
        self._3d_highlight = []
        self._3d_data = []
        if not self.channels:
            return

        # Prefer R-based y-axis. Use only channels that have R; fall back to
        # channel index only if NO channel has an R value.
        any_R = any(ch.get('R') is not None for ch in self.channels)
        if any_R:
            channels_use = [(i, ch) for i, ch in enumerate(self.channels)
                            if ch.get('R') is not None]
            channels_use.sort(key=lambda p: float(p[1]['R']))
            y_label = 'R [m]'
        else:
            channels_use = list(enumerate(self.channels))
            y_label = 'Channel index'

        ref_time = np.asarray(channels_use[0][1]['time'])
        n_t = len(ref_time)
        n_ch = len(channels_use)
        y_vals = np.asarray([float(ch.get('R')) if ch.get('R') is not None else j
                             for j, (_, ch) in enumerate(channels_use)])
        # Map sorted-position → original channel index for slider highlight
        self._3d_ch_order = [orig_i for orig_i, _ in channels_use]

        for col, (key, plabel) in enumerate([
            ('param1', self.param1_label), ('param2', self.param2_label)
        ]):
            ax = self.figure.add_subplot(1, 2, col + 1, projection='3d')
            self._axes.append(ax)
            Z = np.full((n_ch, n_t), np.nan)
            for i, (_orig_i, ch) in enumerate(channels_use):
                arr = ch.get(key)
                t = np.asarray(ch['time'])
                if arr is None:
                    continue
                arr = np.asarray(arr, dtype=float)
                if len(t) == n_t:
                    Z[i, :] = arr[:n_t]
                else:
                    try:
                        Z[i, :] = np.interp(ref_time, t, arr)
                    except Exception:
                        pass
            T_grid, Y_grid = np.meshgrid(ref_time, y_vals, indexing='xy')
            try:
                ax.plot_surface(T_grid, Y_grid, Z, cmap='viridis',
                                edgecolor='none', alpha=0.6)
            except Exception:
                ax.plot_wireframe(T_grid, Y_grid, Z, color='#1f77b4', lw=0.5)
            ax.set_xlabel('Time [s]', fontsize=10, labelpad=2)
            ax.set_ylabel(y_label, fontsize=10, labelpad=2)
            ax.set_zlabel(plabel, fontsize=10, labelpad=2)
            ax.tick_params(labelsize=8, pad=2)
            _shrink_3d_box(ax)
            _layout_3d_axes(ax, col, 2)
            # Robust z-axis bounds. Z shape: (n_ch_along_R, n_time);
            # transpose so that the radial axis becomes axis 0 for edge skipping.
            param_key = ''.join(c for c in plabel if c.isalpha())
            is_zero_bottom = any(k in param_key for k in _ZERO_BOTTOM)
            is_vt = ('vT' in param_key) or ('vt' in param_key)
            z_lo, z_hi = _robust_zrange(Z)
            if z_hi is not None and np.isfinite(z_hi):
                if is_zero_bottom:
                    z_lo_use = 0.0
                elif is_vt:
                    z_lo_use = _VT_YMIN
                else:
                    z_lo_use = z_lo if z_lo is not None else 0
                ax.set_zlim(z_lo_use, z_hi * 1.05 if z_hi > 0 else z_hi * 0.95)
            self._3d_data.append((ref_time, y_vals, Z))
            self._3d_highlight.append(None)

        self.figure.suptitle(
            f"#{self.shot_number}  ({self.source})  3D view", fontsize=12)
        self.figure.subplots_adjust(left=0.0, right=1.0, top=0.94, bottom=0.0,
                                    wspace=0.05)
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()

    def _render_3d_highlight(self, ch_idx):
        """Highlight time-trace at the current channel slice on each 3D axis."""
        if not getattr(self, '_3d_data', None):
            return
        if ch_idx >= len(self.channels):
            return
        label = self.channels[ch_idx].get('label', f'idx {ch_idx}')
        # Sync navigation widgets
        if hasattr(self, 'idx_entry'):
            self.idx_entry.setText(str(ch_idx))
        if hasattr(self, 'ch_info_label'):
            self.ch_info_label.setText(label)
        # Map slider index (original channel order) → sorted-by-R position
        order = getattr(self, '_3d_ch_order', None)
        if order is None:
            sorted_pos = ch_idx
        else:
            try:
                sorted_pos = order.index(ch_idx)
            except ValueError:
                sorted_pos = None
        for col, dataset in enumerate(self._3d_data):
            ax = self._axes[col]
            old = self._3d_highlight[col]
            if old is not None:
                try:
                    old.remove()
                except Exception:
                    pass
            if dataset is None or sorted_pos is None:
                continue
            ref_time, y_vals, Z = dataset
            ci = min(sorted_pos, Z.shape[0] - 1)
            y_val = float(y_vals[ci])
            line, = ax.plot(ref_time, np.full_like(ref_time, y_val), Z[ci, :],
                            color='red', lw=2.5, zorder=20)
            self._3d_highlight[col] = line
        self.figure.suptitle(
            f"#{self.shot_number}  {label}  ({self.source})  3D view",
            fontsize=12)
        self.canvas.draw_idle()

    # ---- Add to Selected ----

    def _add_current(self):
        if not self.channels:
            return
        ch = self.channels[self.slider.value()]
        entry = ch['entry']

        if self.selected_listbox is not None:
            existing = {self.selected_listbox.item(i).text()
                        for i in range(self.selected_listbox.count())}
            if entry in existing:
                self.add_btn.setText(f"{ch['label']} already selected")
                return
            self.selected_listbox.addItem(entry)

        self.added_count += 1
        self.add_btn.setText(f"Added {ch['label']}  ({self.added_count} total)")


# ============================================================
# TranspProfilePreviewDialog
# ============================================================

class TranspProfilePreviewDialog(_PreviewBase):
    """Preview TRANSP CDF profiles: filter + variable selector + time slider."""

    _HAS_DALPHA = False  # TRANSP CDF data has no Dα plot

    def __init__(self, parent, cdf, shot, run, selected_listbox=None):
        self.cdf = cdf
        self.shot = shot
        self.run = run
        self._prof_vars = sorted(cdf['profiles'].keys())
        self._all_var_items = []  # [(var_name, display_text), ...]
        for vn in self._prof_vars:
            vd = cdf['profiles'][vn]
            display = vn
            if vd['long_name'] and vd['long_name'] != vn:
                display += f" - {vd['long_name']}"
            self._all_var_items.append((vn, display))
        self._current_prof = None

        super().__init__(parent, 1,
                         f"Browse  #{shot} ({run})",
                         selected_listbox)
        self._build_ui()
        self._init_plot()
        if self._prof_vars:
            self.var_combo.setCurrentIndex(0)
            self._on_var_changed()

    def _build_ui(self):
        from PySide6.QtWidgets import QComboBox as _QCB

        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(4)

        self._build_left(main_layout, "Use mouse wheel or arrow keys to navigate time points.")

        def _nav(layout):
            # --- Variable Selection group ---
            var_group = QGroupBox("Variable")
            var_layout = QVBoxLayout(var_group)

            frow = QHBoxLayout()
            frow.addWidget(QLabel("Filter"))
            self.var_filter = QLineEdit()
            self.var_filter.setPlaceholderText("Type to filter...")
            self.var_filter.textChanged.connect(self._update_filter)
            frow.addWidget(self.var_filter, 1)
            var_layout.addLayout(frow)

            self.var_combo = _QCB()
            self.var_combo.setMaxVisibleItems(20)
            self.var_combo.setSizeAdjustPolicy(_QCB.AdjustToMinimumContentsLengthWithIcon)
            self.var_combo.setMinimumContentsLength(8)
            self.var_combo.setStyleSheet("QComboBox { combobox-popup: 0; }")
            self.var_combo.view().setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
            for vn, display in self._all_var_items:
                self.var_combo.addItem(display, vn)
            self.var_combo.currentIndexChanged.connect(self._on_var_changed)
            var_layout.addWidget(self.var_combo)
            layout.addWidget(var_group)

            # --- Navigation group ---
            nav_group = QGroupBox("Navigation")
            nav_layout = QVBoxLayout(nav_group)

            row1 = QHBoxLayout()
            row1.addWidget(QLabel("Frame"))
            self.frame_entry = QLineEdit("1")
            self.frame_entry.returnPressed.connect(self._goto_frame)
            row1.addWidget(self.frame_entry)
            self.frame_max_label = QLabel("/ 0")  # updated in _on_var_changed
            row1.addWidget(self.frame_max_label)
            go = QPushButton()
            go.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
            go.setFixedSize(24, 24)
            go.clicked.connect(self._goto_frame)
            row1.addWidget(go)
            nav_layout.addLayout(row1)

            row2 = QHBoxLayout()
            row2.addWidget(QLabel("Time [s]"))
            self.time_entry = QLineEdit("0.0")
            self.time_entry.returnPressed.connect(self._goto_time)
            row2.addWidget(self.time_entry)
            go2 = QPushButton()
            go2.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
            go2.setFixedSize(24, 24)
            go2.clicked.connect(self._goto_time)
            row2.addWidget(go2)
            nav_layout.addLayout(row2)
            layout.addWidget(nav_group)

        self._build_right(main_layout, _nav)

    def _update_filter(self):
        filt = self.var_filter.text().strip().upper()
        current_var = self.var_combo.currentData()
        self.var_combo.blockSignals(True)
        self.var_combo.clear()
        for vn, display in self._all_var_items:
            if filt and filt not in vn.upper() and filt not in display.upper():
                continue
            self.var_combo.addItem(display, vn)
        if current_var:
            for i in range(self.var_combo.count()):
                if self.var_combo.itemData(i) == current_var:
                    self.var_combo.setCurrentIndex(i)
                    break
        self.var_combo.blockSignals(False)
        if self.var_combo.count() > 0:
            self._on_var_changed()

    def _on_var_changed(self):
        vn = self.var_combo.currentData()
        if vn is None or vn not in self.cdf['profiles']:
            return
        self._current_prof = self.cdf['profiles'][vn]
        n = len(self._current_prof['time'])
        self._n_items = n
        self.slider.setMaximum(max(0, n - 1))
        self.frame_max_label.setText(f"/ {n}")
        self.slider.setValue(0)
        # Re-render surface in 3D mode; otherwise just update the 2D slice
        if getattr(self, '_view_mode', '2D') == '3D':
            self._render_3d()
            self._render_3d_highlight(0)
            self.canvas.draw_idle()
        else:
            self._update_plot(0)

    def _goto_frame(self):
        try:
            self.slider.setValue(
                max(0, min(self.slider.maximum(), int(self.frame_entry.text()) - 1)))
        except ValueError:
            pass

    def _goto_time(self):
        if self._current_prof is None:
            return
        try:
            t = float(self.time_entry.text())
            self.slider.setValue(
                int(np.argmin(np.abs(self._current_prof['time'] - t))))
        except ValueError:
            pass

    def _init_plot(self):
        self.figure.clear()
        self._axes = []
        ax = self.figure.add_subplot(1, 1, 1)
        ax.tick_params(labelsize=10)
        ax.grid(ls='--', lw=0.3, color='#444444')
        self._axes.append(ax)
        self.figure.subplots_adjust(
            left=0.12, right=0.97, top=0.90, bottom=0.14)
        apply_dark_figure_style(self.figure)

    def _update_plot(self, t_idx):
        prof = self._current_prof
        if prof is None:
            return
        if t_idx >= len(prof['time']):
            return

        t_actual = prof['time'][t_idx]
        ms = int(round(t_actual * 1000))
        vn = self.var_combo.currentData()

        # Use per-time shot labels if available
        shot_labels = self.cdf.get('shot_labels')
        if shot_labels and t_idx < len(shot_labels):
            title_shot = shot_labels[t_idx]
        else:
            title_shot = f"#{self.shot} ({self.run})"
        self.figure.suptitle(
            f"{title_shot}  {ms:06d}ms  —  {vn}", fontsize=12)
        self.frame_entry.setText(str(t_idx + 1))
        self.time_entry.setText(f"{t_actual:.4f}")

        for a in self._artists:
            a.remove()
        self._artists = []

        ax = self._axes[0]
        y = prof['data'][t_idx, :]
        x = prof['x']
        valid = np.isfinite(y)
        if np.any(valid):
            line, = ax.plot(x[valid], y[valid], '-', color='#1f77b4', lw=1.5, marker='.', markersize=4)
            self._artists.append(line)

        ax.set_xlabel(r'$\rho_{tor}$', fontsize=12)
        units = prof.get('units', '')
        ax.set_ylabel(f"{vn} [{units}]" if units else vn, fontsize=12)

        if not self._fix_axes and np.any(valid):
            yv = y[valid]
            margin = (np.nanmax(yv) - np.nanmin(yv)) * 0.1 or 0.1
            ax.set_ylim(np.nanmin(yv) - margin, np.nanmax(yv) + margin)
            ax.set_xlim(x[0], x[-1])

        self.canvas.draw_idle()

    def _render_3d(self):
        """Surface plot of variable vs (rho, time) for the selected variable."""
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        self.figure.clear()
        self._axes = []
        self._3d_highlight = None
        self._3d_data = None
        prof = self._current_prof
        if prof is None:
            return

        ax = self.figure.add_subplot(1, 1, 1, projection='3d')
        self._axes.append(ax)
        x = np.asarray(prof['x'])
        time = np.asarray(prof['time'])
        Z = np.asarray(prof['data'])  # (n_time, n_x)
        X_grid, T_grid = np.meshgrid(x, time, indexing='xy')
        try:
            ax.plot_surface(X_grid, T_grid, Z, cmap='viridis',
                            edgecolor='none', alpha=0.6)
        except Exception:
            ax.plot_wireframe(X_grid, T_grid, Z, color='#1f77b4', lw=0.5)
        ax.set_xlabel(r'$\rho_{tor}$', fontsize=10, labelpad=2)
        ax.set_ylabel('Time [s]', fontsize=10, labelpad=2)
        vn = self.var_combo.currentData() or ''
        units = prof.get('units', '')
        ax.set_zlabel(f"{vn} [{units}]" if units else vn,
                      fontsize=10, labelpad=2)
        ax.tick_params(labelsize=8, pad=2)
        _shrink_3d_box(ax)
        _layout_3d_axes(ax, 0, 1)
        # Robust z-axis bounds (TRANSP profiles: time on axis 0; transpose for edge skip)
        is_zero_bottom = any(k in vn for k in ('NE', 'TE', 'TI', 'NH', 'ND', 'NT', 'PRES', 'Q'))
        z_lo, z_hi = _robust_zrange(Z.T)  # radial axis becomes axis 0
        if z_hi is not None and np.isfinite(z_hi):
            if is_zero_bottom:
                z_lo_use = 0.0
            else:
                z_lo_use = z_lo if z_lo is not None else 0
            ax.set_zlim(z_lo_use, z_hi * 1.05 if z_hi > 0 else z_hi * 0.95)
        self._3d_data = (x, time, Z)
        self.figure.suptitle(f"#{self.shot} ({self.run})  3D view  —  {vn}",
                             fontsize=12)
        self.figure.subplots_adjust(left=0.0, right=1.0, top=0.94, bottom=0.0)
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()

    def _render_3d_highlight(self, t_idx):
        """Highlight profile slice at time t_idx on the 3D surface."""
        if not self._axes or self._3d_data is None:
            return
        ax = self._axes[0]
        if self._3d_highlight is not None:
            try:
                self._3d_highlight.remove()
            except Exception:
                pass
        x, time, Z = self._3d_data
        ti = min(t_idx, Z.shape[0] - 1)
        t_val = float(time[ti])
        line, = ax.plot(x, np.full_like(x, t_val), Z[ti, :],
                        color='red', lw=2.5, zorder=20)
        self._3d_highlight = line
        ms = int(round(t_val * 1000))
        vn = self.var_combo.currentData() or ''
        # Sync navigation widgets
        if hasattr(self, 'frame_entry'):
            self.frame_entry.setText(str(t_idx + 1))
        if hasattr(self, 'time_entry'):
            self.time_entry.setText(f"{t_val:.4f}")
        self.figure.suptitle(
            f"#{self.shot} ({self.run})  {ms:06d}ms  3D view  —  {vn}",
            fontsize=12)
        self.canvas.draw_idle()

    def _add_current(self):
        prof = self._current_prof
        if prof is None:
            return
        t_idx = self.slider.value()
        t_val = prof['time'][t_idx]
        ms = int(round(t_val * 1000))
        entry = f"{self.shot}_{ms:06d} ({self.run})"

        if self.selected_listbox is not None:
            existing = {self.selected_listbox.item(i).text()
                        for i in range(self.selected_listbox.count())}
            if entry in existing:
                self.add_btn.setText(f"{ms:06d}ms already selected")
                return
            self.selected_listbox.addItem(entry)

        self.added_count += 1
        self.add_btn.setText(f"Added {ms:06d}ms  ({self.added_count} total)")


# ============================================================
# TranspTimeTracePreviewDialog
# ============================================================

class TranspTimeTracePreviewDialog(QDialog):
    """Preview TRANSP CDF time traces: filter + variable selector.
    Plots all loaded runs for the selected variable at once. No slider needed.
    """

    def __init__(self, parent, cdf_cache, run_labels):
        super().__init__(parent)
        self.cdf_cache = cdf_cache
        self.run_labels = run_labels
        self._all_var_items = self._collect_vars()

        import re
        title_parts = []
        for label in run_labels:
            m = re.match(r'^(\d+)(.*)', label)
            title_parts.append(f"#{m.group(1)} ({m.group(2)})" if m else label)

        self.setWindowTitle(f"Browse  {', '.join(title_parts)}")
        self.resize(1000, 550)
        self._build_ui()
        self._init_plot()
        if self._all_var_items:
            self.var_combo.setCurrentIndex(0)
            self._on_var_changed()

    def _collect_vars(self):
        all_vars = {}
        for cdf in self.cdf_cache.values():
            for vn, vd in cdf['timetraces'].items():
                if vn not in all_vars:
                    parts = [vn]
                    if vd['long_name'] and vd['long_name'] != vn:
                        parts.append(f"- {vd['long_name']}")
                    if vd['units']:
                        parts.append(f"[{vd['units']}]")
                    all_vars[vn] = ' '.join(parts)
        return sorted(all_vars.items())

    def _build_ui(self):
        from PySide6.QtWidgets import QComboBox as _QCB
        from ui.widgets.toggle_switch import ToggleSwitch

        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(4)

        # Left: Canvas
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        self.figure = Figure((8, 4.5))
        apply_dark_figure_style(self.figure)
        self.canvas = FigureCanvasQTAgg(self.figure)
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
        toolbar = NavigationToolbar2QT(self.canvas, self)
        left_layout.addWidget(toolbar)
        left_layout.addWidget(self.canvas, stretch=1)
        main_layout.addWidget(left, stretch=1)

        # Right: Controls
        ctrl = QWidget()
        ctrl.setFixedWidth(240)
        cl = QVBoxLayout(ctrl)
        cl.setContentsMargins(4, 0, 4, 4)

        # Variable selection
        var_group = QGroupBox("Variable")
        var_layout = QVBoxLayout(var_group)

        frow = QHBoxLayout()
        frow.addWidget(QLabel("Filter"))
        self.var_filter = QLineEdit()
        self.var_filter.setPlaceholderText("Type to filter...")
        self.var_filter.textChanged.connect(self._update_filter)
        frow.addWidget(self.var_filter, 1)
        var_layout.addLayout(frow)

        self.var_combo = _QCB()
        self.var_combo.setMaxVisibleItems(20)
        self.var_combo.setSizeAdjustPolicy(_QCB.AdjustToMinimumContentsLengthWithIcon)
        self.var_combo.setMinimumContentsLength(8)
        self.var_combo.setStyleSheet("QComboBox { combobox-popup: 0; }")
        self.var_combo.view().setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        for vn, display in self._all_var_items:
            self.var_combo.addItem(display, vn)
        self.var_combo.currentIndexChanged.connect(self._on_var_changed)
        var_layout.addWidget(self.var_combo)
        cl.addWidget(var_group)

        # Options
        opt_group = QGroupBox("Options")
        opt_layout = QHBoxLayout(opt_group)
        self._fix_axes = False
        self.fix_axes_toggle = ToggleSwitch()
        self.fix_axes_toggle.toggled.connect(self._on_fix_axes)
        opt_layout.addWidget(self.fix_axes_toggle)
        opt_layout.addWidget(QLabel("Fix Axes"))
        opt_layout.addStretch()
        cl.addWidget(opt_group)

        # Runs info
        info_group = QGroupBox("Runs")
        info_layout = QVBoxLayout(info_group)
        import re
        for label in self.run_labels:
            m = re.match(r'^(\d+)(.*)', label)
            text = f"#{m.group(1)} ({m.group(2)})" if m else label
            info_layout.addWidget(QLabel(text))
        cl.addWidget(info_group)

        cl.addStretch()

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        cl.addWidget(close_btn)

        main_layout.addWidget(ctrl)

    def _update_filter(self):
        filt = self.var_filter.text().strip().upper()
        current_var = self.var_combo.currentData()
        self.var_combo.blockSignals(True)
        self.var_combo.clear()
        for vn, display in self._all_var_items:
            if filt and filt not in vn.upper() and filt not in display.upper():
                continue
            self.var_combo.addItem(display, vn)
        if current_var:
            for i in range(self.var_combo.count()):
                if self.var_combo.itemData(i) == current_var:
                    self.var_combo.setCurrentIndex(i)
                    break
        self.var_combo.blockSignals(False)
        if self.var_combo.count() > 0:
            self._on_var_changed()

    def _on_var_changed(self):
        self._plot()

    def _on_fix_axes(self, checked):
        self._fix_axes = checked
        if not checked:
            self._plot()

    def _init_plot(self):
        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.tick_params(labelsize=10)
        ax.grid(ls='--', lw=0.3, color='#444444')
        self._ax = ax
        self.figure.subplots_adjust(
            left=0.12, right=0.97, top=0.90, bottom=0.14)
        apply_dark_figure_style(self.figure)

    def _plot(self):
        import re
        vn = self.var_combo.currentData()
        if vn is None:
            return

        ax = self._ax
        ax.clear()
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.tick_params(labelsize=10)
        ax.grid(ls='--', lw=0.3, color='#444444')

        import matplotlib.pyplot as plt
        cmap = plt.get_cmap('tab10')
        colors = [cmap(i % 10) for i in range(len(self.run_labels))]

        for idx, label in enumerate(self.run_labels):
            cdf = self.cdf_cache.get(label)
            if cdf is None or vn not in cdf['timetraces']:
                continue
            tt = cdf['timetraces'][vn]
            m = re.match(r'^(\d+)(.*)', label)
            legend = f"#{m.group(1)} ({m.group(2)})" if m else label
            valid = np.isfinite(tt['data'])
            if np.any(valid):
                if len(tt['data']) <= 1:
                    ax.plot(tt['time'][valid], tt['data'][valid], 'o',
                            color=colors[idx], markersize=8, label=legend)
                else:
                    ax.plot(tt['time'][valid], tt['data'][valid], '-',
                            color=colors[idx], lw=1.5, marker='.', markersize=5, label=legend)

        units = ''
        for cdf in self.cdf_cache.values():
            if vn in cdf['timetraces']:
                units = cdf['timetraces'][vn].get('units', '')
                break
        ax.set_ylabel(f"{vn} [{units}]" if units else vn, fontsize=12)

        if ax.get_legend_handles_labels()[1]:
            ax.legend(fontsize=8, loc='best', frameon=False)

        self.figure.suptitle(vn, fontsize=12)
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()


# ============================================================
# BiProfilePreviewDialog
# ============================================================

_BI_UNITS = {'Ti': 'keV', 'vT': 'km/s', 'Te': 'keV', 'ne': r'$10^{19}$/m$^3$'}
_BI_LABELS = {'Ti': r'T$_i$', 'vT': r'v$_T$', 'Te': r'T$_e$', 'ne': r'n$_e$'}


class BiProfilePreviewDialog(QDialog):
    """Non-modal dialog for browsing BiProfile data.

    mode='profile': slider browses time, plots profile vs psi_N
    mode='timetrace': slider browses psi_N, plots time trace
    """

    def __init__(self, parent, bi_data, shot_number, params, mode='profile',
                 selected_listbox=None, sdata=None, plot_color='#1f77b4'):
        super().__init__(parent)
        self.bi_data = bi_data
        self.shot_number = shot_number
        self.params = params
        self.mode = mode
        self.selected_listbox = selected_listbox
        self.sdata = sdata
        self.plot_color = plot_color
        self.added_count = 0

        self._artists = []
        self._axes = []
        self._fix_axes = False

        ref = next((p for p in params if p in bi_data), None)
        d = bi_data[ref]
        if mode == 'profile':
            self.browse_arr = d['time']
        else:
            self.browse_arr = d['psin']

        self._play_timer = QTimer()
        self._play_timer.timeout.connect(self._play_step)
        self._playing = False

        mode_str = "Profile" if mode == 'profile' else "TimeTrace"
        self.setWindowTitle(
            f"Browse {mode_str}  #{shot_number} ({', '.join(params)})")
        self.resize(1000, 550)
        self._build_ui()
        self._init_plot()
        self._update_plot(0)

    def _build_ui(self):
        from ui.widgets.toggle_switch import ToggleSwitch

        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(4)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)

        self.figure = Figure((8, 4.5))
        apply_dark_figure_style(self.figure)
        self.canvas = FigureCanvasQTAgg(self.figure)
        toolbar = NavigationToolbar2QT(self.canvas, self)
        left_layout.addWidget(toolbar)
        left_layout.addWidget(self.canvas, stretch=1)

        slider_row = QHBoxLayout()
        prev_btn = QPushButton()
        prev_btn.setIcon(get_icon(QStyle.SP_ArrowBack))
        prev_btn.setFixedWidth(28)
        prev_btn.clicked.connect(lambda: self._step(-1))
        slider_row.addWidget(prev_btn)
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(len(self.browse_arr) - 1)
        self.slider.valueChanged.connect(self._on_slider)
        slider_row.addWidget(self.slider, stretch=1)
        next_btn = QPushButton()
        next_btn.setIcon(get_icon(QStyle.SP_ArrowForward))
        next_btn.setFixedWidth(28)
        next_btn.clicked.connect(lambda: self._step(1))
        slider_row.addWidget(next_btn)
        left_layout.addLayout(slider_row)

        hint = QLabel("Use mouse wheel or arrow keys to navigate.")
        hint.setStyleSheet("color: gray; font-size: 8pt;")
        hint.setAlignment(Qt.AlignCenter)
        left_layout.addWidget(hint)
        main_layout.addWidget(left, stretch=1)

        ctrl = QWidget()
        ctrl.setFixedWidth(240)
        cl = QVBoxLayout(ctrl)
        cl.setContentsMargins(4, 0, 4, 4)

        # View mode (top)
        from PySide6.QtWidgets import QRadioButton, QButtonGroup
        view_group = QGroupBox("View")
        view_layout = QHBoxLayout(view_group)
        self._view_mode = '2D'
        self._view_btn_group = QButtonGroup(self)
        self.view_2d_radio = QRadioButton("2D")
        self.view_3d_radio = QRadioButton("3D")
        self.view_2d_radio.setChecked(True)
        self._view_btn_group.addButton(self.view_2d_radio)
        self._view_btn_group.addButton(self.view_3d_radio)
        view_layout.addWidget(self.view_2d_radio)
        view_layout.addWidget(self.view_3d_radio)
        view_layout.addStretch()
        self.view_2d_radio.toggled.connect(
            lambda c: c and self._on_view_changed("2D"))
        self.view_3d_radio.toggled.connect(
            lambda c: c and self._on_view_changed("3D"))
        cl.addWidget(view_group)

        nav_group = QGroupBox("Navigation")
        nav_layout = QVBoxLayout(nav_group)
        idx_row = QHBoxLayout()
        idx_row.addWidget(QLabel("Index"))
        self.idx_entry = QLineEdit("0")
        self.idx_entry.returnPressed.connect(self._goto_idx)
        idx_row.addWidget(self.idx_entry)
        idx_row.addWidget(QLabel(f"/ {len(self.browse_arr)-1}"))
        go_idx = QPushButton()
        go_idx.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
        go_idx.setFixedSize(24, 24)
        go_idx.clicked.connect(self._goto_idx)
        idx_row.addWidget(go_idx)
        nav_layout.addLayout(idx_row)
        val_row = QHBoxLayout()
        val_label = "Time [s]" if self.mode == 'profile' else "\u03c8\u2099"
        val_row.addWidget(QLabel(val_label))
        self.val_entry = QLineEdit("0.0")
        self.val_entry.returnPressed.connect(self._goto_val)
        val_row.addWidget(self.val_entry)
        go_val = QPushButton()
        go_val.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
        go_val.setFixedSize(24, 24)
        go_val.clicked.connect(self._goto_val)
        val_row.addWidget(go_val)
        nav_layout.addLayout(val_row)
        cl.addWidget(nav_group)

        play_group = QGroupBox("Playback")
        play_layout = QHBoxLayout(play_group)
        self.speed_combo = QComboBox()
        self.speed_combo.addItems(['1 fps', '2 fps', '5 fps', '10 fps', '20 fps', '50 fps'])
        self.speed_combo.setCurrentText('5 fps')
        play_layout.addWidget(self.speed_combo)
        self.play_btn = QPushButton()
        self.play_btn.setIcon(get_icon(QStyle.SP_MediaPlay))
        self.play_btn.clicked.connect(self._toggle_play)
        play_layout.addWidget(self.play_btn)
        stop_btn = QPushButton()
        stop_btn.setIcon(get_icon(QStyle.SP_MediaStop))
        stop_btn.clicked.connect(self._stop_play)
        play_layout.addWidget(stop_btn)
        cl.addWidget(play_group)

        opt_group = QGroupBox("Options")
        opt_layout = QVBoxLayout(opt_group)
        fix_row = QHBoxLayout()
        self.fix_axes_toggle = ToggleSwitch()
        self.fix_axes_toggle.toggled.connect(self._on_fix_axes)
        fix_row.addWidget(self.fix_axes_toggle)
        fix_row.addWidget(QLabel("Fix Axes"))
        fix_row.addStretch()
        opt_layout.addLayout(fix_row)
        self._apply_scale = True
        if self.mode == 'profile' and ('ne' in self.params or 'Te' in self.params):
            scale_row = QHBoxLayout()
            self.scale_toggle = ToggleSwitch()
            self.scale_toggle.toggled.connect(self._on_scale_toggled)
            scale_row.addWidget(self.scale_toggle)
            scale_row.addWidget(QLabel("Unapply TS ne Scale"))
            scale_row.addStretch()
            opt_layout.addLayout(scale_row)
        cl.addWidget(opt_group)

        cl.addStretch()
        self.add_btn = QPushButton("Add to Selected")
        self.add_btn.clicked.connect(self._add_current)
        cl.addWidget(self.add_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self._on_close)
        cl.addWidget(close_btn)
        main_layout.addWidget(ctrl)

    def _init_plot(self):
        self.figure.clear()
        self._axes = []
        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'
        for col, param in enumerate(self.params):
            ax = self.figure.add_subplot(1, 2, col + 1)
            if self.mode == 'profile':
                ax.set_xlabel(r'$\psi_N$', fontsize=12)
                ax.set_xlim(0, 1.1)
                ax.axvline(x=1.0, color=zc, ls='--', lw=0.8, alpha=0.5)
            else:
                ax.set_xlabel('Time [s]', fontsize=12)
            ax.set_ylabel(f'{_BI_LABELS.get(param, param)} [{_BI_UNITS.get(param, "")}]',
                          fontsize=12)
            ax.tick_params(labelsize=10)
            ax.grid(ls='--', lw=0.3, color='#444444')
            if param == 'vT':
                ax.axhline(y=0, color=zc, ls='--', gid='zero_ref')
            self._axes.append(ax)
        self.figure.subplots_adjust(
            left=0.12, right=0.97, top=0.90, bottom=0.14, wspace=0.30)
        apply_dark_figure_style(self.figure)

    def _update_plot(self, idx):
        val = self.browse_arr[idx]
        bi = self.bi_data
        color = self.plot_color

        if self.mode == 'profile':
            ms = int(round(val * 1000))
            title = f"#{self.shot_number}  {ms:06d}ms"
        else:
            title = f"#{self.shot_number}  \u03c8\u2099={val:.4f}"

        self.figure.suptitle(title, fontsize=12)
        self.idx_entry.setText(str(idx))
        self.val_entry.setText(f"{val:.4f}")

        for a in self._artists:
            a.remove()
        self._artists = []

        for col, param in enumerate(self.params):
            ax = self._axes[col]
            if param not in bi:
                continue
            d = bi[param]

            if self.mode == 'profile':
                psin = d['psin']
                t_idx = np.argmin(np.abs(d['time'] - val))
                mean_p = d['mean'][:, t_idx]
                unc_p = d['unc'][:, t_idx]
                valid = ~np.isnan(mean_p)
                if np.any(valid):
                    line, = ax.plot(psin[valid], mean_p[valid], '-', color=color, lw=2)
                    self._artists.append(line)
                    fill = ax.fill_between(psin[valid],
                                           mean_p[valid] - unc_p[valid],
                                           mean_p[valid] + unc_p[valid],
                                           alpha=0.2, color=color)
                    self._artists.append(fill)
                    if not self._fix_axes:
                        y_min = np.nanmin(mean_p[valid] - unc_p[valid])
                        y_max = np.nanmax(mean_p[valid] + unc_p[valid])
                        margin = (y_max - y_min) * 0.1 if y_max > y_min else 0.1
                        if param in _ZERO_BOTTOM:
                            ax.set_ylim(0, y_max + margin)
                        elif param == 'vT':
                            ax.set_ylim(_VT_YMIN, y_max + margin)
                        else:
                            ax.set_ylim(y_min - margin, y_max + margin)
                if self.sdata:
                    self._overlay_raw(ax, param, val, color)
            else:
                time = d['time']
                mean_t = d['mean'][idx, :]
                unc_t = d['unc'][idx, :]
                valid = ~np.isnan(mean_t)
                if np.any(valid):
                    line, = ax.plot(time[valid], mean_t[valid], '-', color=color, lw=2)
                    self._artists.append(line)
                    fill = ax.fill_between(time[valid],
                                           mean_t[valid] - unc_t[valid],
                                           mean_t[valid] + unc_t[valid],
                                           alpha=0.2, color=color)
                    self._artists.append(fill)
                    if not self._fix_axes:
                        y_min = np.nanmin(mean_t[valid] - unc_t[valid])
                        y_max = np.nanmax(mean_t[valid] + unc_t[valid])
                        margin = (y_max - y_min) * 0.1 if y_max > y_min else 0.1
                        if param in _ZERO_BOTTOM:
                            ax.set_ylim(0, y_max + margin)
                        elif param == 'vT':
                            ax.set_ylim(_VT_YMIN, y_max + margin)
                        else:
                            ax.set_ylim(y_min - margin, y_max + margin)
                    ax.set_xlim(time[0], time[-1])

        self.canvas.draw_idle()

    # ---- Navigation ----

    def _on_slider(self, value):
        if getattr(self, '_view_mode', '2D') == '3D':
            self._render_3d_highlight(value)
        else:
            self._update_plot(value)

    def _step(self, delta):
        self.slider.setValue(
            max(0, min(self.slider.maximum(), self.slider.value() + delta)))

    def _goto_idx(self):
        try:
            self.slider.setValue(
                max(0, min(len(self.browse_arr) - 1, int(self.idx_entry.text()))))
        except ValueError:
            pass

    def _goto_val(self):
        try:
            self.slider.setValue(
                int(np.argmin(np.abs(self.browse_arr - float(self.val_entry.text())))))
        except ValueError:
            pass

    def _on_fix_axes(self, checked):
        self._fix_axes = checked
        if not checked:
            self._update_plot(self.slider.value())

    # ---- View mode (2D / 3D) ----

    def _on_view_changed(self, mode):
        self._view_mode = mode
        if hasattr(self, 'fix_axes_toggle'):
            self.fix_axes_toggle.setEnabled(mode == '2D')
        if mode == '3D':
            self._render_3d()
            self._render_3d_highlight(self.slider.value())
        else:
            self._init_plot()
            self._update_plot(self.slider.value())

    def _render_3d(self):
        """Surface plot of param vs (psi_N, time)."""
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        self.figure.clear()
        self._axes = []
        self._3d_highlight = []
        self._3d_data = []
        bi = self.bi_data

        for col, param in enumerate(self.params):
            ax = self.figure.add_subplot(1, len(self.params), col + 1,
                                          projection='3d')
            self._axes.append(ax)
            if param not in bi:
                self._3d_data.append(None)
                self._3d_highlight.append(None)
                continue
            d = bi[param]
            psin = np.asarray(d['psin'])
            time = np.asarray(d['time'])
            Z = np.asarray(d['mean'])  # shape (n_psin, n_time)
            X, Y = np.meshgrid(psin, time, indexing='ij')
            try:
                ax.plot_surface(X, Y, Z, cmap='viridis',
                                edgecolor='none', alpha=0.6)
            except Exception:
                ax.plot_wireframe(X, Y, Z, color='#1f77b4', lw=0.5)
            ax.set_xlabel(r'$\psi_N$', fontsize=10, labelpad=2)
            ax.set_ylabel('Time [s]', fontsize=10, labelpad=2)
            ax.set_zlabel(
                f'{_BI_LABELS.get(param, param)} [{_BI_UNITS.get(param, "")}]',
                fontsize=10, labelpad=2)
            ax.tick_params(labelsize=8, pad=2)
            _shrink_3d_box(ax)
            _layout_3d_axes(ax, col, len(self.params))
            # Robust z-axis bounds (skip outer 10%, percentile-based)
            z_lo, z_hi = _robust_zrange(Z)
            if z_hi is not None and np.isfinite(z_hi):
                if param in _ZERO_BOTTOM:
                    z_lo_use = 0.0
                elif param == 'vT':
                    z_lo_use = _VT_YMIN
                else:
                    z_lo_use = z_lo if z_lo is not None else 0
                ax.set_zlim(z_lo_use, z_hi * 1.05 if z_hi > 0 else z_hi * 0.95)
            self._3d_data.append((psin, time, Z))
            self._3d_highlight.append(None)

        self.figure.suptitle(f"#{self.shot_number}  3D view", fontsize=12)
        self.figure.subplots_adjust(left=0.0, right=1.0, top=0.94, bottom=0.0,
                                    wspace=0.05)
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()

    def _render_3d_highlight(self, idx):
        """Highlight a slice (time or psi_N) on the 3D surface."""
        if not getattr(self, '_3d_data', None):
            return
        if idx >= len(self.browse_arr):
            return
        val = float(self.browse_arr[idx])
        # Sync navigation widgets
        if hasattr(self, 'idx_entry'):
            self.idx_entry.setText(str(idx))
        if hasattr(self, 'val_entry'):
            self.val_entry.setText(f"{val:.4f}")
        if self.mode == 'profile':
            ms = int(round(val * 1000))
            title = f"#{self.shot_number}  {ms:06d}ms  3D view"
        else:
            title = f"#{self.shot_number}  ψₙ={val:.4f}  3D view"

        for col, dataset in enumerate(self._3d_data):
            ax = self._axes[col]
            old = self._3d_highlight[col]
            if old is not None:
                try:
                    old.remove()
                except Exception:
                    pass
            if dataset is None:
                continue
            psin, time, Z = dataset
            if self.mode == 'profile':
                ti = int(np.argmin(np.abs(time - val)))
                line, = ax.plot(psin, np.full_like(psin, float(time[ti])),
                                Z[:, ti], color='red', lw=2.5, zorder=20)
            else:
                pi = int(np.argmin(np.abs(psin - val)))
                line, = ax.plot(np.full_like(time, float(psin[pi])), time,
                                Z[pi, :], color='red', lw=2.5, zorder=20)
            self._3d_highlight[col] = line

        self.figure.suptitle(title, fontsize=12)
        self.canvas.draw_idle()

    # ---- Playback ----

    def _get_fps(self):
        return int(self.speed_combo.currentText().split()[0])

    def _toggle_play(self):
        if self._playing:
            self._stop_play()
        else:
            self._playing = True
            self.play_btn.setIcon(get_icon(QStyle.SP_MediaPause))
            self._play_timer.start(max(1, 1000 // self._get_fps()))

    def _stop_play(self):
        self._playing = False
        self.play_btn.setIcon(get_icon(QStyle.SP_MediaPlay))
        self._play_timer.stop()

    def _play_step(self):
        if self.slider.value() >= self.slider.maximum():
            self._stop_play()
            return
        self.slider.setValue(self.slider.value() + 1)
        self._play_timer.setInterval(max(1, 1000 // self._get_fps()))

    # ---- Raw overlay ----

    def _on_scale_toggled(self, checked):
        self._apply_scale = not checked
        self._update_plot(self.slider.value())

    def _overlay_raw(self, ax, param, t_actual, color):
        from data_loaders.biprofile_loader import map_R_to_psin
        sdata = self.sdata
        efit = sdata.get('efit')
        if efit is None:
            return
        diag = sdata.get('diag')

        if param in ['Ti', 'vT']:
            raw = sdata.get('ces')
            if raw is None:
                return
            t_idx = np.argmin(np.abs(raw['time'] - t_actual))
            R = raw['R']
            psin_mapped = map_R_to_psin(R, efit, t_actual)
            y = raw['Ti'][:, t_idx] if param == 'Ti' else raw['vT'][:, t_idx]
            yerr = raw['Ti_err'][:, t_idx] if param == 'Ti' else raw['vT_err'][:, t_idx]
            use_flags = self._get_use_flags(diag, param, len(R))
            self._plot_raw(ax, psin_mapped, y, yerr, use_flags, color, 'o')
        elif param in ['Te', 'ne']:
            raw = sdata.get('thomson')
            if raw is None:
                return
            t_idx = np.argmin(np.abs(raw['time'] - t_actual))
            R = raw['R']
            psin_mapped = map_R_to_psin(R, efit, t_actual)
            y = raw['Te'][:, t_idx] if param == 'Te' else raw['ne'][:, t_idx]
            yerr = raw['Te_err'][:, t_idx] if param == 'Te' else raw['ne_err'][:, t_idx]
            if param == 'ne' and self._apply_scale and diag:
                scale = self._get_ts_scale(diag, t_actual, len(R))
                y = y * scale
                yerr = yerr * scale
            use_flags = self._get_use_flags(diag, param, len(R))
            self._plot_raw(ax, psin_mapped, y, yerr, use_flags, color, 's')

    @staticmethod
    def _get_use_flags(diag, param, n_ch):
        flags = np.full(n_ch, np.nan)
        if diag is None:
            return flags
        def _pick(arr):
            if arr is None: return np.nan
            if not isinstance(arr, np.ndarray): return float(arr)
            if arr.ndim == 0: return float(arr)
            return float(arr[0])
        if param in ['Ti', 'vT']:
            key = 'ti' if param == 'Ti' else 'vt'
            for i, ch in enumerate(diag.get('ces', {}).get(key, [])):
                if i >= n_ch: break
                flags[i] = _pick(ch.get(f'{key}_use'))
        elif param in ['ne', 'Te']:
            use_key = 'ne_use' if param == 'ne' else 'te_use'
            idx = 0
            for region in ['core', 'edge']:
                for ch in diag.get('thomson', {}).get(region, []):
                    if idx >= n_ch: break
                    flags[idx] = _pick(ch.get(use_key))
                    idx += 1
        return flags

    def _get_ts_scale(self, diag, t_actual, n_ch):
        scales = np.ones(n_ch)
        ts = diag.get('thomson', {})
        bi = self.bi_data
        ref = next((p for p in ('ne', 'Te') if p in bi), None)
        if ref is None:
            return scales
        t_idx = np.argmin(np.abs(bi[ref]['time'] - t_actual))
        idx = 0
        for region in ['core', 'edge']:
            for ch in ts.get(region, []):
                if idx >= n_ch: break
                scale_arr = ch.get('scale')
                if scale_arr is not None:
                    if isinstance(scale_arr, np.ndarray) and len(scale_arr) > 0:
                        scales[idx] = scale_arr[min(t_idx, len(scale_arr) - 1)]
                    elif not isinstance(scale_arr, np.ndarray):
                        scales[idx] = float(scale_arr)
                idx += 1
        return scales

    def _plot_raw(self, ax, psin, y, yerr, use_flags, color, marker):
        yerr = np.where(yerr < 0, 0, yerr)
        valid = (psin >= 0) & (psin <= 1.2) & np.isfinite(y)
        used = valid & (use_flags == 1)
        excluded = valid & (use_flags == 0)
        unknown = valid & np.isnan(use_flags)
        if np.any(used):
            c = ax.errorbar(psin[used], y[used], yerr=yerr[used],
                            fmt=marker, markersize=5, color=color,
                            capsize=5, zorder=10, markeredgewidth=1)
            self._artists.append(c)
        if np.any(excluded):
            c = ax.errorbar(psin[excluded], y[excluded], yerr=yerr[excluded],
                            fmt=marker, markersize=5, color=(0.6, 0.6, 0.6, 0.35),
                            capsize=5, zorder=5, markeredgewidth=1)
            self._artists.append(c)
        if np.any(unknown):
            c = ax.errorbar(psin[unknown], y[unknown], yerr=yerr[unknown],
                            fmt=marker, markersize=5, color=color,
                            capsize=5, zorder=10, markeredgewidth=1)
            self._artists.append(c)

    # ---- Add to Selected ----

    def _add_current(self):
        idx = self.slider.value()
        val = self.browse_arr[idx]
        if self.mode == 'profile':
            ms = int(round(val * 1000))
            entry = f"{self.shot_number:06d}_{val*1e3:06.0f} (Bi)"
            display = f"{ms:06d}ms"
        else:
            entry = f"{self.shot_number:06d}_{idx:03d} (Bi)"
            display = f"idx={idx} (\u03c8\u2099={val:.4f})"
        if self.selected_listbox is not None:
            existing = {self.selected_listbox.item(i).text()
                        for i in range(self.selected_listbox.count())}
            if entry in existing:
                self.add_btn.setText(f"{display} already selected")
                return
            self.selected_listbox.addItem(entry)
        self.added_count += 1
        self.add_btn.setText(f"Added {display}  ({self.added_count} total)")

    def _on_close(self):
        self._stop_play()
        self.accept()


# ============================================================
# BiProfileDerivedPreviewDialog
# ============================================================

class BiProfileDerivedPreviewDialog(_PreviewBase):
    """Browse derived plasma quantities for the currently-loaded shot.

    Variable dropdown lets the user pick any quantity from the derived catalog;
    time slider scans the biprofile time grid. Mirrors TranspProfilePreviewDialog
    UX. 'Add to Selected' pushes `{shot:06d}_{ms:06d}` entries into the parent
    derived tab's selected_listbox.
    """

    _HAS_DALPHA = False   # no Dα data in derived dialog

    def __init__(self, parent, derived_computer, shot, selected_listbox=None):
        from core.derived_quantities import PROFILE_QUANTITIES
        self.derived = derived_computer
        self.shot = shot
        self.time = np.asarray(derived_computer.time)
        self.psin = np.asarray(derived_computer.psin)
        n = len(self.time)
        self._all_var_items = [
            (key, f"{key}  —  {name}  [{unit}]")
            for key, name, _lbl, unit, _req in PROFILE_QUANTITIES
        ]
        self._current_data = None

        super().__init__(parent, n, f"Browse Derived  #{shot}", selected_listbox)
        self._build_ui()
        self._init_plot()
        if self._all_var_items:
            self.var_combo.setCurrentIndex(0)
            self._on_var_changed()

    def _build_ui(self):
        from PySide6.QtWidgets import QComboBox as _QCB

        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(4)

        self._build_left(main_layout,
                         "Use mouse wheel or arrow keys to navigate time points.")

        def _nav(layout):
            var_group = QGroupBox("Variable")
            var_layout = QVBoxLayout(var_group)
            self.var_combo = _QCB()
            self.var_combo.setMaxVisibleItems(20)
            self.var_combo.setSizeAdjustPolicy(_QCB.AdjustToMinimumContentsLengthWithIcon)
            self.var_combo.setMinimumContentsLength(8)
            self.var_combo.setStyleSheet("QComboBox { combobox-popup: 0; }")
            self.var_combo.view().setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
            for vn, display in self._all_var_items:
                self.var_combo.addItem(display, vn)
            self.var_combo.currentIndexChanged.connect(self._on_var_changed)
            var_layout.addWidget(self.var_combo)
            layout.addWidget(var_group)

            nav_group = QGroupBox("Navigation")
            nav_layout = QVBoxLayout(nav_group)
            row1 = QHBoxLayout()
            row1.addWidget(QLabel("Frame"))
            self.frame_entry = QLineEdit("1")
            self.frame_entry.returnPressed.connect(self._goto_frame)
            row1.addWidget(self.frame_entry)
            self.frame_max_label = QLabel(f"/ {len(self.time)}")
            row1.addWidget(self.frame_max_label)
            go = QPushButton()
            go.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
            go.setFixedSize(24, 24)
            go.clicked.connect(self._goto_frame)
            row1.addWidget(go)
            nav_layout.addLayout(row1)

            row2 = QHBoxLayout()
            row2.addWidget(QLabel("Time [s]"))
            self.time_entry = QLineEdit(f"{self.time[0]:.3f}" if len(self.time) else "0.0")
            self.time_entry.returnPressed.connect(self._goto_time)
            row2.addWidget(self.time_entry)
            go2 = QPushButton()
            go2.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
            go2.setFixedSize(24, 24)
            go2.clicked.connect(self._goto_time)
            row2.addWidget(go2)
            nav_layout.addLayout(row2)
            layout.addWidget(nav_group)

        self._build_right(main_layout, _nav)

    def _on_var_changed(self):
        key = self.var_combo.currentData()
        if key is None:
            return
        try:
            d = self.derived.compute(key)
        except Exception as e:
            import traceback; traceback.print_exc()
            d = None
        self._current_data = d
        if d is None:
            self.figure.clear()
            ax = self.figure.add_subplot(1, 1, 1)
            ax.text(0.5, 0.5, f"{key}: not computable\n(check console)",
                    ha='center', va='center', fontsize=12, color='gray',
                    transform=ax.transAxes)
            ax.axis('off')
            apply_dark_figure_style(self.figure)
            self.canvas.draw_idle()
            return
        # Preserve current view mode (2D / 3D) when switching variable
        if self._view_mode == '3D':
            self._render_3d()
            self._render_3d_highlight(self.slider.value())
            self.canvas.draw_idle()
        else:
            self._update_plot(self.slider.value())

    def _init_plot(self):
        self.figure.clear()
        self.canvas.draw_idle()

    def _global_ylim(self):
        """Return (ymin, ymax) over ALL times for the current quantity."""
        d = self._current_data
        if d is None:
            return None
        vals = d['value']
        valid = np.isfinite(vals)
        if not valid.any():
            return None
        ymin = float(np.nanmin(vals[valid]))
        ymax = float(np.nanmax(vals[valid]))
        if ymin == ymax:
            ymin -= 1; ymax += 1
        pad = 0.05 * (ymax - ymin)
        return (ymin - pad, ymax + pad)

    def _update_plot(self, idx):
        d = self._current_data
        if d is None or idx >= len(self.time):
            return
        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)
        psin = d['psin']
        vals = d['value'][:, idx]
        valid = np.isfinite(vals)
        if np.any(valid):
            ax.plot(psin[valid], vals[valid], 'o-',
                    color='#1f77b4', lw=1.6, ms=3)
        ax.set_xlabel(r'$\psi_\mathrm{N}$', fontsize=11)
        # ylabel: include full name
        from core.derived_quantities import get_quantity_meta
        key = self.var_combo.currentData() if hasattr(self, 'var_combo') else ''
        meta = get_quantity_meta(key) if key else None
        if meta:
            qname, qlabel, qunit, _ = meta
            ax.set_ylabel(f'{qname}, {qlabel}  [{qunit}]', fontsize=11)
        else:
            ax.set_ylabel(f"{d['label']} [{d['unit']}]", fontsize=11)
        ax.set_xlim(0, 1.1)
        ax.grid(ls='--', lw=0.3, color='#444444')
        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'
        ax.axvline(x=1.0, color=zc, ls='--', lw=0.8, alpha=0.5)
        if key and key.startswith('nu_star'):
            ax.set_yscale('log')
        # Fix Axes: clamp y to the global range across all times for this var
        if self._fix_axes:
            yl = self._global_ylim()
            if yl is not None:
                ax.set_ylim(*yl)
        t_s = float(self.time[idx])
        ax.set_title(f"#{self.shot}  t = {t_s*1e3:.1f} ms  ({idx+1}/{len(self.time)})",
                     fontsize=10)
        if hasattr(self, 'frame_entry'):
            self.frame_entry.setText(str(idx + 1))
        if hasattr(self, 'time_entry'):
            self.time_entry.setText(f"{t_s:.3f}")
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()

    def _render_3d(self):
        """3D surface: ψ_N × time × quantity."""
        d = self._current_data
        if d is None:
            super()._render_3d()
            return
        try:
            from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
        except Exception:
            super()._render_3d()
            return
        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1, projection='3d')
        psin = np.asarray(d['psin'])
        time = np.asarray(d['time'])
        vals = np.asarray(d['value'])   # (n_psi, n_time)
        # Mask non-finite to avoid surface artifacts
        Z = np.where(np.isfinite(vals), vals, np.nan)
        P, T = np.meshgrid(psin, time, indexing='ij')
        try:
            ax.plot_surface(P, T, Z, cmap='viridis',
                            edgecolor='none', antialiased=False,
                            rstride=1, cstride=1)
        except Exception:
            ax.text2D(0.5, 0.5, "3D render failed", ha='center',
                      va='center', transform=ax.transAxes)
        from core.derived_quantities import get_quantity_meta
        key = self.var_combo.currentData() if hasattr(self, 'var_combo') else ''
        meta = get_quantity_meta(key) if key else None
        if meta:
            qname, qlabel, qunit, _ = meta
            zlabel = f'{qname}, {qlabel}  [{qunit}]'
        else:
            zlabel = f"{d['label']} [{d['unit']}]"
        ax.set_xlabel(r'$\psi_\mathrm{N}$')
        ax.set_ylabel('Time [s]')
        ax.set_zlabel(zlabel)
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()

    def _render_3d_highlight(self, idx):
        """Overlay a highlight curve on the 3D surface for the current time."""
        d = self._current_data
        if d is None or not self.figure.axes:
            return
        ax = self.figure.axes[0]
        # Remove previous highlight artists tagged with gid='3d_hl'
        for artist in list(ax.lines):
            if getattr(artist, 'get_gid', lambda: None)() == '3d_hl':
                artist.remove()
        psin = np.asarray(d['psin'])
        time = np.asarray(d['time'])
        vals = np.asarray(d['value'])
        if idx < 0 or idx >= len(time):
            return
        z = vals[:, idx]
        valid = np.isfinite(z)
        if not np.any(valid):
            return
        t_const = np.full_like(psin[valid], time[idx])
        ax.plot(psin[valid], t_const, z[valid], color='red', lw=2.0,
                gid='3d_hl')
        self.canvas.draw_idle()

    def _goto_frame(self):
        try:
            self.slider.setValue(max(
                0, min(self.slider.maximum(), int(self.frame_entry.text()) - 1)))
        except ValueError:
            pass

    def _goto_time(self):
        try:
            t = float(self.time_entry.text())
            self.slider.setValue(int(np.argmin(np.abs(self.time - t))))
        except ValueError:
            pass

    def _add_current(self):
        idx = self.slider.value()
        if idx >= len(self.time):
            return
        t_s = float(self.time[idx])
        entry = f"{self.shot:06d}_{int(round(t_s*1e3)):06d}"
        display = f"{int(round(t_s*1e3)):06d}ms"
        if self.selected_listbox is None:
            self.add_btn.setText(f"{display}: no parent list")
            return
        existing = {self.selected_listbox.item(i).text()
                    for i in range(self.selected_listbox.count())}
        if entry in existing:
            self.add_btn.setText(f"{display} already selected")
            return
        self.selected_listbox.addItem(entry)
        self.added_count += 1
        self.add_btn.setText(f"Added {display}  ({self.added_count} total)")
