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

_ZERO_BOTTOM = {'Ti', 'Te', 'ne'}


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

        # Options
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
        self._update_plot(value)

    def _step(self, delta):
        self.slider.setValue(
            max(0, min(self.slider.maximum(), self.slider.value() + delta)))

    def _on_fix_axes(self, checked):
        self._fix_axes = checked
        if not checked:
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
        self.time_arr = getattr(data, 'time_prof', data.time)

        super().__init__(parent, len(self.time_arr),
                         f"Browse  #{shot_number} ({source})",
                         selected_listbox)
        self._build_ui()
        self._init_plot()
        self._update_plot(0)

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
            self.slider.setValue(
                max(0, min(len(self.time_arr) - 1, int(self.frame_entry.text()) - 1)))
        except ValueError:
            pass

    def _goto_time(self):
        try:
            self.slider.setValue(
                int(np.argmin(np.abs(self.time_arr - float(self.time_entry.text())))))
        except ValueError:
            pass

    # ---- Plot ----

    def _init_plot(self):
        self.figure.clear()
        self._axes = []
        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'
        for col, (pname, plabel) in enumerate([
            (self.param1_name, self.param1_label),
            (self.param2_name, self.param2_label),
        ]):
            ax = self.figure.add_subplot(1, 2, col + 1)
            ax.set_xlabel('R [m]', fontsize=12)
            ax.set_ylabel(plabel, fontsize=12)
            ax.tick_params(labelsize=10)
            ax.grid(ls='--', lw=0.3, color='#444444')
            if 'v' in pname.lower():
                ax.axhline(y=0, color=zc, ls='--', gid='zero_ref')
            self._axes.append(ax)
        self.figure.subplots_adjust(
            left=0.12, right=0.97, top=0.90, bottom=0.14, wspace=0.30)
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
                else:
                    ax.set_ylim(y_min - margin, y_max + margin)

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
                    ax2.set_ylim(np.nanmin(y) - margin, np.nanmax(y) + margin)

        if not self._fix_axes:
            ax1.set_xlim(time[0], time[-1])

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
            self.var_combo.setSizeAdjustPolicy(_QCB.AdjustToMinimumContentsLengthWithIcon)
            self.var_combo.setMinimumContentsLength(8)
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
        self.var_combo.setSizeAdjustPolicy(_QCB.AdjustToMinimumContentsLengthWithIcon)
        self.var_combo.setMinimumContentsLength(8)
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
                        else:
                            ax.set_ylim(y_min - margin, y_max + margin)
                    ax.set_xlim(time[0], time[-1])

        self.canvas.draw_idle()

    # ---- Navigation ----

    def _on_slider(self, value):
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
