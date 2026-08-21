"""
TRANSP UFILE — single tab combining Time Trace, Profile, 2D views.

Pick one TRANSP UFILE run directory, choose any UFILE from the dropdown.
Plot dispatches by classification:
    time_trace        → 1D line vs time
    time_trace_multi  → NBI/EC sources (NBI gets 3-panel Pinj/Vacc/fractions)
    profile_time      → profile lines for each time (legend)
    geometry          → LIM (R, Z) parametric / generic 2D imshow / 3D MMX
                        flux surfaces; LIM auto-overlays on every plot
"""

import os
import numpy as np

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QScrollArea, QLayout,
    QGroupBox, QLabel, QPushButton, QFileDialog, QApplication,
    QMessageBox, QSpinBox, QDoubleSpinBox, QSlider, QDialog,
    QDialogButtonBox, QComboBox, QRadioButton, QButtonGroup,
)
from PySide6.QtCore import Qt

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style


# ============================================================
# Helpers
# ============================================================

def _channel_count_from_scalars(uf):
    for s in uf.get('SCALARS', []) or []:
        nm = (s.get('name') or '').upper()
        if any(nm.startswith(p) for p in
               ('NBEAM', 'NGYRO', 'NSOURCE', 'NSYS')):
            try:
                return int(round(float(s.get('data'))))
            except Exception:
                pass
    return None


def _detect_nbi_packing(uf):
    F = uf.get('F', {}).get('data')
    if F is None or F.ndim != 2:
        return None
    n = F.shape[1]
    if n < 4 or n % 4 != 0:
        return None
    quarter = n // 4
    last_max = np.nanmax(np.abs(F[:, -quarter:])) if F.size else 0.0
    if 0 < last_max <= 1.5:
        return quarter
    return None


def _resolve_n_sources(ext, uf):
    n = _channel_count_from_scalars(uf)
    if n is not None:
        return n
    if ext.upper() in ('NBI', 'NB2'):
        n = _detect_nbi_packing(uf)
        if n is not None:
            return n
    F = uf.get('F', {}).get('data')
    n_total = F.shape[1] if (F is not None and F.ndim == 2) else 0
    from data_loaders.transp_loader import get_catalog
    n_quants = int(get_catalog().get(ext.upper(), {}).get('n_quants', 1) or 1)
    if n_quants > 0 and n_total % n_quants == 0:
        return n_total // n_quants
    return n_total


# ============================================================
# Tab
# ============================================================

class TranspUFileTab:
    """Single-tab UFILE viewer covering all classes."""

    def __init__(self, main_window):
        self.main = main_window
        self._run = None
        self._last_dir = None
        self._cax = None

        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10
        self._color_mode = 'Gradient(viridis)'
        # Profile view mode: '2D' (lines + colorbar + legend) or '3D' (surface)
        self._view_3d = False
        # Max legend entries before truncating with ellipsis
        self._legend_max_items = 20
        self._settings_key = "transp_ufile"

        self.frame = QWidget()
        self.canvas = None
        # Load persisted settings BEFORE widgets are built so radio buttons,
        # combos, etc. reflect the restored state (was: created with defaults
        # then settings overwrote the model — UI fell out of sync).
        self._restore_settings()
        self._create_widgets()

    # ============================================================
    # Widgets
    # ============================================================

    def _create_widgets(self):
        self.figure = Figure((10, 6))
        apply_dark_figure_style(self.figure)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.draw()

        canvas_widget = QWidget()
        cl = QVBoxLayout(canvas_widget)
        cl.setContentsMargins(0, 0, 0, 0)
        cl.addWidget(self.canvas)

        # Re-snap colorbar to compressed-box right edge on resize
        def _on_resize(_event):
            ax = self.figure.axes[0] if self.figure.axes else None
            if ax is not None and self._cax is not None:
                self._reposition_colorbar(ax)
        self.figure.canvas.mpl_connect('resize_event', _on_resize)

        scroll = QScrollArea()
        scroll.setFixedWidth(CONTROL_PANEL_WIDTH)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        control = QWidget()
        control_layout = QVBoxLayout(control)
        control_layout.setContentsMargins(9, 9, 9, 9)
        control_layout.setSizeConstraint(QLayout.SetMinimumSize)
        scroll.setWidget(control)
        scroll.viewport().setAutoFillBackground(False)
        control.setAutoFillBackground(False)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(canvas_widget)
        splitter.addWidget(scroll)

        main_layout = QVBoxLayout(self.frame)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(splitter)

        self._create_load_section(control)
        self._create_var_section(control)
        self._create_plot_controls(control)
        self._create_legend_section(control)
        self._create_time_section(control)
        control_layout.addStretch()

    def _create_load_section(self, parent):
        group = QGroupBox("1. Load UFILE Run")
        layout = QVBoxLayout(group)
        btn = QPushButton("Open Run Directory...")
        btn.clicked.connect(self._open_run_dir)
        layout.addWidget(btn)
        parent.layout().addWidget(group)

    def _create_var_section(self, parent):
        group = QGroupBox("2. Select UFILE Variable")
        group.setEnabled(False)
        self._var_group = group
        layout = QVBoxLayout(group)

        self.var_combo = QComboBox()
        self.var_combo.setMaxVisibleItems(20)
        self.var_combo.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.var_combo.setMinimumContentsLength(10)
        self.var_combo.setMaximumWidth(CONTROL_PANEL_WIDTH - 30)
        self.var_combo.setStyleSheet("QComboBox { combobox-popup: 0; }")
        self.var_combo.view().setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.var_combo.currentIndexChanged.connect(self._on_var_changed)
        layout.addWidget(self.var_combo)

        parent.layout().addWidget(group)

    def _create_plot_controls(self, parent):
        group = QGroupBox("3. Plot")
        layout = QVBoxLayout(group)
        btn_row = QHBoxLayout()
        plot_btn = QPushButton("Plot")
        plot_btn.clicked.connect(self._plot)
        btn_row.addWidget(plot_btn, 3)
        opt_btn = QPushButton("Style")
        opt_btn.clicked.connect(self._show_style_dialog)
        btn_row.addWidget(opt_btn, 1)
        layout.addLayout(btn_row)
        parent.layout().addWidget(group)

    def _create_legend_section(self, parent):
        group = QGroupBox("4. View")
        group.setStyleSheet(
            "QGroupBox:disabled { color: #888; }"
            "QGroupBox:disabled::title { color: #888; }"
        )
        group.setEnabled(False)
        self._legend_group = group
        layout = QVBoxLayout(group)
        layout.setContentsMargins(9, 2, 9, 9)

        hint = QLabel("active when a profile UFILE is selected")
        hint.setStyleSheet("color: #888; font-size: 9pt;")
        layout.addWidget(hint)

        radio_row = QHBoxLayout()
        self.view_btn_group = QButtonGroup(group)
        self.view_2d_radio = QRadioButton("2D")
        self.view_3d_radio = QRadioButton("3D")
        self.view_btn_group.addButton(self.view_2d_radio, 0)
        self.view_btn_group.addButton(self.view_3d_radio, 1)
        if self._view_3d:
            self.view_3d_radio.setChecked(True)
        else:
            self.view_2d_radio.setChecked(True)
        self.view_btn_group.buttonClicked.connect(self._on_view_changed)
        radio_row.addWidget(self.view_2d_radio)
        radio_row.addWidget(self.view_3d_radio)
        radio_row.addStretch()
        layout.addLayout(radio_row)

        parent.layout().addWidget(group)

    def _on_view_changed(self, btn):
        self._view_3d = (btn is self.view_3d_radio)
        if self.figure.axes and self._current_class() == 'profile_time':
            self._plot()

    def _on_legend_mode_changed(self, btn):
        self._profile_legend_mode = (
            'colorbar' if btn is self.colorbar_radio else 'legend')
        if self.figure.axes and self._current_class() == 'profile_time':
            self._plot()

    def _create_time_section(self, parent):
        group = QGroupBox("5. Time Slice")
        group.setStyleSheet(
            "QGroupBox:disabled { color: #888; }"
            "QGroupBox:disabled::title { color: #888; }"
        )
        group.setEnabled(False)
        self._time_group = group
        layout = QVBoxLayout(group)
        layout.setContentsMargins(9, 2, 9, 9)

        hint = QLabel("active when a 3D UFILE (e.g., MMX) is selected")
        hint.setStyleSheet("color: #888; font-size: 9pt;")
        layout.addWidget(hint)

        row = QHBoxLayout()
        row.addWidget(QLabel("t [s]"))
        self.time_spin = QDoubleSpinBox()
        self.time_spin.setDecimals(3)
        self.time_spin.setRange(0.0, 1000.0)
        self.time_spin.setSingleStep(0.05)
        self.time_spin.setFixedWidth(90)
        self.time_spin.valueChanged.connect(self._on_time_changed)
        row.addWidget(self.time_spin)
        layout.addLayout(row)

        self.time_slider = QSlider(Qt.Horizontal)
        self.time_slider.setRange(0, 1000)
        self.time_slider.valueChanged.connect(self._on_slider_changed)
        layout.addWidget(self.time_slider)

        self.time_status = QLabel("")
        self.time_status.setStyleSheet("color: #888; font-size: 9pt;")
        self.time_status.setWordWrap(True)
        layout.addWidget(self.time_status)
        parent.layout().addWidget(group)

    # ============================================================
    # Load / status
    # ============================================================

    def _open_run_dir(self):
        start = self._last_dir or os.path.expanduser("~")
        dlg = QFileDialog(self.frame, "Select TRANSP UFILE directory", start)
        dlg.setFileMode(QFileDialog.Directory)
        dlg.setOption(QFileDialog.ShowDirsOnly, True)
        dlg.setOption(QFileDialog.DontUseNativeDialog, True)
        dlg.setWindowModality(Qt.NonModal)
        dlg.fileSelected.connect(self._on_run_selected)
        dlg.show()

    def _on_run_selected(self, path):
        from data_loaders.transp_loader import scan_run_directory
        self._last_dir = path
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            run = scan_run_directory(path)
        finally:
            QApplication.restoreOverrideCursor()
        if run is None:
            QMessageBox.information(
                self.frame, "No UFILEs",
                "No `<prefix><shot>.<ext>` files found in:\n" + path)
            return
        self._run = run
        self._var_group.setEnabled(True)
        self._update_var_combo()
        self._on_var_changed()
        self._update_status()

    def _update_status(self):
        from ui.ui_constants import show_status
        if self._run is None:
            show_status(self.frame, 'UFILE', "No run loaded", 'gray')
        else:
            r = self._run
            show_status(self.frame, 'UFILE',
                        f"#{r['shot']}  {r['prefix']} · {r['label']}  —  "
                        f"{len(r['ufiles'])} files", 'green')

    # ============================================================
    # Variable combo — all classes, sorted by group
    # ============================================================

    # Physical-category sort order: 평형 → profile → 가열 → 그 외(history)
    _CAT_ORDER = {
        'equilibrium': 0,
        'profile':     1,
        'heating':     2,
        'history':     3,
    }

    def _update_var_combo(self):
        from data_loaders.transp_loader import pretty_label, get_catalog
        catalog = get_catalog()
        current = self.var_combo.currentData()
        items = []
        if self._run:
            for ext, _cls in self._run['classes'].items():
                # LIM auto-overlays on geometry plots — don't list separately
                if ext.upper() == 'LIM':
                    continue
                cat = catalog.get(ext.upper(), {}).get('category', 'unknown')
                cat_rank = self._CAT_ORDER.get(cat, 9)
                label = pretty_label(ext, self._run['ufiles'][ext])
                items.append((cat_rank, ext, label))

        # Sort by physical category, then by extension alphabetically
        items.sort(key=lambda it: (it[0], it[1]))

        self.var_combo.blockSignals(True)
        self.var_combo.clear()
        for _rank, ext, label in items:
            self.var_combo.addItem(label, ext)
        restore = current or getattr(self, '_pending_variable', None)
        if restore:
            for i in range(self.var_combo.count()):
                if self.var_combo.itemData(i) == restore:
                    self.var_combo.setCurrentIndex(i)
                    break
            self._pending_variable = None
        self.var_combo.blockSignals(False)

    def _current_uf(self):
        if self._run is None:
            return None
        ext = self.var_combo.currentData()
        if not ext or ext not in self._run['ufiles']:
            return None
        return self._run['ufiles'][ext]

    def _current_class(self):
        if self._run is None:
            return None
        ext = self.var_combo.currentData()
        return self._run['classes'].get(ext)

    # ============================================================
    # Var change → toggle time slider (3D only) + auto re-plot
    # ============================================================

    @staticmethod
    def _extract_time(uf):
        for k in ('X0', 'X1', 'X2'):
            ax = uf.get(k, {})
            name = (ax.get('name') or '').upper()
            units = (ax.get('units') or '').upper()
            if 'TIME' in name or units in ('S', 'SEC', 'SECS', 'SECOND',
                                           'SECONDS', 'MS', 'MSEC'):
                return ax.get('data')
        return uf.get('X0', {}).get('data')

    def _on_var_changed(self):
        uf = self._current_uf()
        cls = self._current_class()

        # Time Labels group active only for profile_time UFILEs
        is_profile = (cls == 'profile_time')
        self._legend_group.setEnabled(bool(is_profile))

        is_3d = (uf is not None and uf.get('NDIM', 0) >= 3
                 and uf.get('F', {}).get('data') is not None
                 and uf['F']['data'].ndim >= 3)
        self._time_group.setEnabled(bool(is_3d))
        if is_3d:
            t = self._extract_time(uf)
            if t is not None and len(t) > 0:
                tmin, tmax = float(t[0]), float(t[-1])
                self.time_spin.blockSignals(True)
                self.time_slider.blockSignals(True)
                self.time_spin.setRange(tmin, tmax)
                if not (tmin <= self.time_spin.value() <= tmax):
                    self.time_spin.setValue(tmax)
                rng = max(tmax - tmin, 1e-9)
                v = (self.time_spin.value() - tmin) / rng
                self.time_slider.setValue(int(round(v * 1000)))
                self.time_spin.blockSignals(False)
                self.time_slider.blockSignals(False)
                self.time_status.setText(
                    f"range: [{tmin:.3f}, {tmax:.3f}] s  ({len(t)} frames)")
        else:
            self.time_status.setText("")

        # Auto-replot if there's already a plot on the canvas
        if self.figure.axes and self._run is not None and uf is not None:
            self._plot()

    def _on_time_changed(self, _v):
        tmin = self.time_spin.minimum()
        tmax = self.time_spin.maximum()
        rng = max(tmax - tmin, 1e-9)
        self.time_slider.blockSignals(True)
        v = (self.time_spin.value() - tmin) / rng
        self.time_slider.setValue(int(round(v * 1000)))
        self.time_slider.blockSignals(False)
        if self._current_class() == 'geometry':
            uf = self._current_uf()
            if uf is not None and uf.get('NDIM', 0) >= 3:
                self._plot()

    def _on_slider_changed(self, ival):
        tmin = self.time_spin.minimum()
        tmax = self.time_spin.maximum()
        rng = max(tmax - tmin, 1e-9)
        self.time_spin.blockSignals(True)
        self.time_spin.setValue(tmin + (ival / 1000.0) * rng)
        self.time_spin.blockSignals(False)
        if self._current_class() == 'geometry':
            uf = self._current_uf()
            if uf is not None and uf.get('NDIM', 0) >= 3:
                self._plot()

    # ============================================================
    # Color
    # ============================================================

    def _cmap_name(self):
        mode = self._color_mode
        s, e = mode.find('('), mode.find(')')
        return mode[s + 1:e] if s != -1 and e != -1 else 'viridis'

    def _get_plot_colors(self, n):
        import matplotlib.pyplot as plt
        mode = self._color_mode
        s, e = mode.find('('), mode.find(')')
        cmap_name = mode[s + 1:e] if s != -1 and e != -1 else 'viridis'
        cmap = plt.get_cmap(cmap_name)
        if mode.startswith('Fixed'):
            colors = (getattr(cmap, 'colors', None)
                      or [cmap(i / max(n - 1, 1)) for i in range(n)])
            return [colors[i % len(colors)] for i in range(n)]
        if n == 1:
            return [cmap(0.5)]
        return [cmap(i / (n - 1)) for i in range(n)]

    # ============================================================
    # Colorbar (explicit cax, repositioned to compressed box)
    # ============================================================

    def _make_cbar_axes(self):
        cax = self.figure.add_axes([0.93, 0.25, 0.018, 0.50])
        self._cax = cax
        return cax

    def _reposition_colorbar(self, ax):
        cax = getattr(self, '_cax', None)
        if cax is None:
            return
        # Only reposition for aspect='equal' plots (LIM/MMX) where the
        # axes box compresses. For regular axes (profile), leave the
        # colorbar at the explicit position so its tick labels don't run
        # past the figure edge.
        try:
            aspect = ax.get_aspect()
        except Exception:
            aspect = None
        if aspect not in ('equal', 1.0, 1):
            return
        self.canvas.draw()
        pos = ax.get_position()
        cax_w = 0.018
        cax_pad = 0.012
        cax_h = pos.height * 0.75
        cax_left = pos.x1 + cax_pad
        cax_bottom = pos.y0 + (pos.height - cax_h) / 2.0
        cax.set_position([cax_left, cax_bottom, cax_w, cax_h])

    # ============================================================
    # Plot dispatch
    # ============================================================

    def _plot(self):
        uf = self._current_uf()
        if uf is None:
            return
        cls = self._current_class()

        self.figure.clear()
        self._cax = None

        if cls == 'time_trace':
            self._plot_time_trace(uf)
        elif cls == 'time_trace_multi':
            self._plot_multi(uf)
        elif cls == 'profile_time':
            self._plot_profile(uf)
        elif cls == 'geometry':
            self._plot_geometry(uf)
        else:
            ax = self.figure.add_subplot(1, 1, 1)
            ax.text(0.5, 0.5, "(unsupported class)",
                    transform=ax.transAxes, ha='center', va='center',
                    color='#888')

        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()

    # ----- time_trace -----

    def _plot_time_trace(self, uf):
        ax = self.figure.add_subplot(1, 1, 1)
        F = uf['F']['data']
        t_idx = uf.get('_time_axis', 0)
        t = uf.get(f'X{t_idx}', {}).get('data')
        units = (uf['F'].get('units') or '').strip()
        F_name = uf['F'].get('name', '')
        ax.plot(t, F, '.-', lw=1.5, color='#1f77b4', markersize=3)
        ax.set_xlabel('Time [s]', fontsize=self.label_fontsize)
        ax.set_ylabel(f"{F_name}  [{units}]" if units else F_name,
                      fontsize=self.label_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(ls='--', lw=0.3, color='#444444')
        self.figure.subplots_adjust(left=0.12, right=0.97, top=0.95,
                                    bottom=0.10)

    # ----- time_trace_multi -----

    def _plot_multi(self, uf):
        ext = self.var_combo.currentData()
        F = uf['F']['data']
        t_idx = uf.get('_time_axis', 0)
        t = uf.get(f'X{t_idx}', {}).get('data')
        if F.ndim != 2:
            return
        if t_idx == 1:
            F = F.T

        is_nbi_4q = (ext.upper() in ('NBI', 'NB2')
                     and _detect_nbi_packing(uf) is not None)

        if is_nbi_4q:
            n_sources = _detect_nbi_packing(uf)
            axes = self.figure.subplots(2, 1, sharex=True, squeeze=False)[:, 0]
            colors = self._get_plot_colors(max(n_sources, 1))

            # Mean current fractions per source (over time) — shown in
            # the Pinj legend so users see (full, half, third) at a glance.
            cfull_all  = F[:, 2 * n_sources:3 * n_sources]
            chalf_all  = F[:, 3 * n_sources:4 * n_sources]
            with np.errstate(invalid='ignore'):
                cfull_mean  = np.nanmean(cfull_all, axis=0)
                chalf_mean  = np.nanmean(chalf_all, axis=0)
            cthird_mean = np.clip(1.0 - cfull_mean - chalf_mean, 0.0, 1.0)

            # Pinj
            for s in range(n_sources):
                lbl = (f"source {s+1}  "
                       f"({cfull_mean[s]:.2f}, "
                       f"{chalf_mean[s]:.2f}, "
                       f"{cthird_mean[s]:.2f})")
                axes[0].plot(t, F[:, s], '.-', color=colors[s], lw=1.5,
                             markersize=3, label=lbl)
            axes[0].set_ylabel('Pinj  [W]', fontsize=self.label_fontsize)
            axes[0].legend(fontsize=self.legend_fontsize, loc='best',
                           frameon=False, ncol=min(n_sources, 2))

            # Vacc
            for s in range(n_sources):
                axes[1].plot(t, F[:, n_sources + s], '.-',
                             color=colors[s], lw=1.5, markersize=3)
            axes[1].set_ylabel('Vacc  [V]', fontsize=self.label_fontsize)

            for ax_q in axes:
                ax_q.grid(ls='--', lw=0.3, color='#444444')
                # Force tick labels on every panel even with sharex=True
                ax_q.tick_params(labelsize=self.tick_fontsize,
                                 labelbottom=True)
            axes[-1].set_xlabel('Time [s]', fontsize=self.label_fontsize)
            self.figure.subplots_adjust(left=0.12, right=0.97, top=0.97,
                                        bottom=0.08, hspace=0.20)
            return

        # Generic multi-channel (ECP/ECA/ECB): one axes, one line per source
        ax = self.figure.add_subplot(1, 1, 1)
        n_total = F.shape[1]
        n_active = max(1, min(int(_resolve_n_sources(ext, uf)), n_total))
        colors = self._get_plot_colors(max(n_active, 1))
        for ci in range(n_active):
            ax.plot(t, F[:, ci], '.-', color=colors[ci], lw=1.5,
                    markersize=3, label=f"source {ci+1}")
        units = (uf['F'].get('units') or '').strip()
        F_name = uf['F'].get('name', '')
        ax.set_xlabel('Time [s]', fontsize=self.label_fontsize)
        ax.set_ylabel(f"{F_name}  [{units}]" if units else F_name,
                      fontsize=self.label_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(ls='--', lw=0.3, color='#444444')
        if ax.get_legend_handles_labels()[1]:
            ax.legend(fontsize=self.legend_fontsize, loc='best',
                      frameon=False, ncol=2)
        self.figure.subplots_adjust(left=0.12, right=0.97, top=0.95,
                                    bottom=0.10)

    # ----- profile_time -----

    def _plot_profile(self, uf):
        F = uf['F']['data']
        t_idx = uf.get('_time_axis')
        if F is None or F.ndim != 2 or t_idx is None:
            return
        sp_idx = 1 - t_idx
        t = uf.get(f'X{t_idx}', {}).get('data')
        x = uf.get(f'X{sp_idx}', {}).get('data')
        if t is None or x is None:
            return
        F2 = F if F.shape == (len(t), len(x)) else (
            F.T if F.shape == (len(x), len(t)) else None)
        if F2 is None:
            return

        units = (uf['F'].get('units') or '').strip()
        F_name = uf['F'].get('name', '')
        sp_meta = uf.get(f'X{sp_idx}', {})
        x_label = sp_meta.get('name', '')
        if sp_meta.get('units'):
            x_label = f"{x_label} [{sp_meta['units']}]"
        z_label = f"{F_name}  [{units}]" if units else F_name

        # 3D surface view (value vs spatial × time)
        if self._view_3d:
            from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
            ax = self.figure.add_subplot(1, 1, 1, projection='3d')
            X_grid, T_grid = np.meshgrid(np.asarray(x, dtype=float),
                                          np.asarray(t, dtype=float),
                                          indexing='xy')
            try:
                ax.plot_surface(X_grid, T_grid, F2, cmap=self._cmap_name(),
                                edgecolor='none', alpha=0.85)
            except Exception:
                ax.plot_wireframe(X_grid, T_grid, F2, color='#1f77b4',
                                  lw=0.5)
            ax.set_xlabel(x_label or 'X', fontsize=self.label_fontsize,
                          labelpad=4)
            ax.set_ylabel('Time [s]', fontsize=self.label_fontsize,
                          labelpad=4)
            ax.set_zlabel(z_label, fontsize=self.label_fontsize, labelpad=4)
            ax.tick_params(labelsize=self.tick_fontsize, pad=2)
            self.figure.subplots_adjust(left=0.02, right=0.96,
                                        top=0.98, bottom=0.02)
            return

        ax = self.figure.add_subplot(1, 1, 1)
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        from matplotlib.lines import Line2D
        cmap = plt.get_cmap(self._cmap_name())
        n_t = len(t)

        for i in range(n_t):
            color = cmap(i / max(n_t - 1, 1))
            ax.plot(x, F2[i, :], '.-', color=color, lw=1.5, markersize=3,
                    label=f"{float(t[i])*1e3:.0f} ms")

        ax.set_xlabel(x_label or 'X', fontsize=self.label_fontsize)
        ax.set_ylabel(f"{F_name}  [{units}]" if units else F_name,
                      fontsize=self.label_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(ls='--', lw=0.3, color='#444444')

        # Reserve right margin for the always-drawn colorbar
        self.figure.subplots_adjust(left=0.10, right=0.84, top=0.97,
                                    bottom=0.10)

        # Always draw the time colorbar
        norm = Normalize(vmin=float(t[0]) * 1e3,
                         vmax=float(t[-1]) * 1e3)
        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cax = self.figure.add_axes([0.87, 0.18, 0.018, 0.70])
        self._cax = cax
        cb = self.figure.colorbar(sm, cax=cax)
        cb.set_label('Time [ms]', fontsize=self.label_fontsize)
        cb.ax.tick_params(labelsize=self.tick_fontsize)

        # Legend (always shown in 2D mode, with ellipsis for many items)
        handles, labels = ax.get_legend_handles_labels()
        max_n = int(self._legend_max_items)
        if len(labels) > max_n:
            idxs = sorted(set(np.linspace(0, len(labels) - 1,
                                          max_n - 1).astype(int)))
            shown_h = [handles[i] for i in idxs]
            shown_l = [labels[i] for i in idxs]
            shown_h.append(Line2D([0], [0], color='none'))
            shown_l.append(f"…  (+{len(labels) - len(idxs)} more)")
            ax.legend(shown_h, shown_l, fontsize=self.legend_fontsize,
                      loc='best', frameon=False,
                      ncol=2 if len(shown_l) > 12 else 1)
        else:
            ax.legend(handles, labels, fontsize=self.legend_fontsize,
                      loc='best', frameon=False,
                      ncol=2 if n_t > 12 else 1)

    # ----- geometry -----

    def _plot_geometry(self, uf):
        ax = self.figure.add_subplot(1, 1, 1)
        F = uf.get('F', {}).get('data')
        F_name = uf.get('F', {}).get('name', '')
        F_units = (uf.get('F', {}).get('units') or '').strip()
        ndim = uf.get('NDIM', 0)

        if F is None or F.size == 0:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha='center', va='center', color='#888')
        elif ndim == 1 or F.ndim == 1:
            self._render_geo_1d(ax, uf, F, F_name, F_units)
        elif ndim == 2 or F.ndim == 2:
            self._render_geo_2d(ax, uf, F, F_name, F_units)
        elif ndim >= 3 or F.ndim >= 3:
            self._render_geo_3d(ax, uf, F, F_name, F_units)

        # LIM auto-overlay (if present)
        self._overlay_lim(ax)

        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(ls='--', lw=0.3, color='#444444')
        if ax.get_legend_handles_labels()[1]:
            ax.legend(fontsize=self.legend_fontsize, loc='best',
                      frameon=False)

        self.figure.subplots_adjust(left=0.12, right=0.97, top=0.90,
                                    bottom=0.10)
        self._reposition_colorbar(ax)

    def _render_geo_1d(self, ax, uf, F, F_name, F_units):
        x0 = uf.get('X0', {}).get('data')
        x0_name = uf.get('X0', {}).get('name', '')
        x0_units = uf.get('X0', {}).get('units', '')
        xs = (x0 if (x0 is not None and len(x0) == len(F))
              else np.arange(len(F)))
        x0u = (x0_units or '').strip().upper()
        fu = F_units.strip().upper()
        is_rz = (x0u in {'M', 'CM', 'MM'} and fu == x0u
                 and x0 is not None and len(x0) == len(F))
        ax.plot(xs, F, '-', lw=1.4, label=F_name)
        if is_rz:
            ax.set_xlim(1.0, 2.5)
            ax.set_ylim(-1.5, 1.5)
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlabel('R [m]', fontsize=self.label_fontsize)
            ax.set_ylabel('Z [m]', fontsize=self.label_fontsize)
        else:
            ax.set_xlabel(
                f"{x0_name} [{x0_units}]" if x0_units else x0_name,
                fontsize=self.label_fontsize)
            ax.set_ylabel(
                f"{F_name}  [{F_units}]" if F_units else F_name,
                fontsize=self.label_fontsize)

    def _render_geo_2d(self, ax, uf, F, F_name, F_units):
        x0 = uf.get('X0', {}).get('data')
        x0_name = uf.get('X0', {}).get('name', '')
        x0_units = uf.get('X0', {}).get('units', '')
        x1 = uf.get('X1', {}).get('data')
        x1_name = uf.get('X1', {}).get('name', '')
        x1_units = uf.get('X1', {}).get('units', '')
        extent = None
        F2 = F
        if (x0 is not None and len(x0) == F.shape[0]
                and x1 is not None and len(x1) == F.shape[1]):
            extent = [float(x1[0]), float(x1[-1]),
                      float(x0[0]), float(x0[-1])]
        elif (x0 is not None and len(x0) == F.shape[1]
              and x1 is not None and len(x1) == F.shape[0]):
            F2 = F.T
            extent = [float(x0[0]), float(x0[-1]),
                      float(x1[0]), float(x1[-1])]
            x0_name, x1_name = x1_name, x0_name
            x0_units, x1_units = x1_units, x0_units
        u0 = (x0_units or '').strip().upper()
        u1 = (x1_units or '').strip().upper()
        aspect = ('equal' if u0 in {'M', 'CM', 'MM'} and u1 == u0
                  else 'auto')
        im = ax.imshow(F2, aspect=aspect, origin='lower',
                       extent=extent, interpolation='nearest')
        cax = self._make_cbar_axes()
        cb = self.figure.colorbar(im, cax=cax)
        cb.set_label(f"{F_name}  [{F_units}]" if F_units else F_name,
                     fontsize=self.label_fontsize)
        ax.set_xlabel(f"{x1_name} [{x1_units}]" if x1_units else x1_name,
                      fontsize=self.label_fontsize)
        ax.set_ylabel(f"{x0_name} [{x0_units}]" if x0_units else x0_name,
                      fontsize=self.label_fontsize)

    def _render_geo_3d(self, ax, uf, F, F_name, F_units):
        time_axis = None
        for k_idx, k in enumerate(('X0', 'X1', 'X2')):
            ax_meta = uf.get(k, {})
            nm = (ax_meta.get('name') or '').upper()
            un = (ax_meta.get('units') or '').upper()
            if 'TIME' in nm or un in ('S', 'SEC', 'SECS', 'SECOND',
                                      'SECONDS', 'MS', 'MSEC'):
                time_axis = k_idx
                break
        if time_axis is None:
            time_axis = 0
        t = uf.get(f'X{time_axis}', {}).get('data')
        if t is None or len(t) == 0:
            return
        t_target = float(self.time_spin.value())
        ti = int(np.argmin(np.abs(t - t_target)))
        actual_t = float(t[ti])

        slicer = [slice(None), slice(None), slice(None)]
        slicer[time_axis] = ti
        F2 = F[tuple(slicer)]

        rem = [i for i in range(3) if i != time_axis]
        ax_a = uf.get(f'X{rem[0]}', {})
        ax_b = uf.get(f'X{rem[1]}', {})
        name_a = (ax_a.get('name') or '').upper()
        name_b = (ax_b.get('name') or '').upper()
        if ('MOMENT' in name_a or 'INDEX' in name_a) and 'MOMENT' not in name_b:
            F2 = F2.T
            ax_a, ax_b = ax_b, ax_a

        x_radial = ax_a.get('data')
        n_radial = F2.shape[0]
        if x_radial is None or len(x_radial) != n_radial:
            x_radial = np.arange(n_radial)

        rz = self._mmx_to_rz(F2, F_units, x_radial=x_radial)
        if rz is not None:
            R, Z = rz
            self._draw_flux_surfaces(ax, R, Z, x_radial, actual_t)
            return

        # Fallback: raw moment overlay
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        cmap = plt.get_cmap(self._cmap_name())
        n_moments = F2.shape[1]
        for j in range(n_moments):
            color = cmap(j / max(n_moments - 1, 1))
            ax.plot(x_radial, F2[:, j], '-', color=color, lw=1.0)
        sm = ScalarMappable(
            norm=Normalize(vmin=0, vmax=max(n_moments - 1, 1)), cmap=cmap)
        sm.set_array([])
        cax = self._make_cbar_axes()
        cb = self.figure.colorbar(sm, cax=cax)
        cb.set_label(ax_b.get('name', 'moment index'),
                     fontsize=self.label_fontsize)
        ax.set_xlabel(ax_a.get('name', ''), fontsize=self.label_fontsize)
        ax.set_ylabel(f"{F_name}  [{F_units}]" if F_units else F_name,
                      fontsize=self.label_fontsize)
        ax.set_title(f"t = {actual_t*1e3:.0f} ms (raw moments)",
                     fontsize=self.label_fontsize)

    @staticmethod
    def _mmx_to_rz(F2, F_units, x_radial=None, n_theta=256):
        if F2.ndim != 2:
            return None
        n_x, n_mom = F2.shape
        if n_mom < 4 or n_mom % 4 != 0:
            return None
        n_modes = n_mom // 4
        Rcm = F2[:, 0:n_modes].copy()
        Rsm = F2[:, n_modes:2 * n_modes].copy()
        Zcm = F2[:, 2 * n_modes:3 * n_modes].copy()
        Zsm = F2[:, 3 * n_modes:].copy()
        if x_radial is not None and len(x_radial) == n_x:
            x = np.asarray(x_radial, dtype=float).reshape(n_x, 1)
            m_idx = np.arange(n_modes).reshape(1, n_modes)
            scale = np.where(m_idx >= 1, np.power(x, m_idx), 1.0)
            Rcm *= scale; Rsm *= scale; Zcm *= scale; Zsm *= scale
        theta = np.linspace(0.0, 2.0 * np.pi, n_theta)
        m = np.arange(n_modes)
        cos_m = np.cos(np.outer(m, theta))
        sin_m = np.sin(np.outer(m, theta))
        R = Rcm @ cos_m + Rsm @ sin_m
        Z = Zcm @ cos_m + Zsm @ sin_m
        if (F_units or '').strip().upper() == 'CM':
            R = R / 100.0
            Z = Z / 100.0
        if not (np.all(np.isfinite(R)) and np.all(np.isfinite(Z))):
            return None
        return R, Z

    def _draw_flux_surfaces(self, ax, R, Z, x_radial, actual_t):
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        cmap = plt.get_cmap(self._cmap_name())
        n_x = R.shape[0]
        x_arr = (np.asarray(x_radial, dtype=float)
                 if x_radial is not None and len(x_radial) == n_x
                 else np.linspace(0.0, 1.0, n_x))
        idxs = np.unique(np.linspace(0, n_x - 1,
                                     min(20, n_x)).astype(int))
        x_min = float(np.nanmin(x_arr))
        x_max = float(np.nanmax(x_arr))
        norm = Normalize(vmin=x_min, vmax=x_max)
        for ix in idxs:
            color = cmap(norm(float(x_arr[ix])))
            ax.plot(R[ix, :], Z[ix, :], '-', color=color, lw=0.9, alpha=0.9)
        bdy_idx = int(np.argmin(np.abs(x_arr - 1.0)))
        ax.plot(R[bdy_idx, :], Z[bdy_idx, :], 'k-', lw=1.6,
                label='boundary')

        ax.set_xlabel('R [m]', fontsize=self.label_fontsize)
        ax.set_ylabel('Z [m]', fontsize=self.label_fontsize)
        ax.set_xlim(1.0, 2.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_aspect('equal', adjustable='box')
        ax.set_title(f"t = {actual_t*1e3:.0f} ms",
                     fontsize=self.label_fontsize)

        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cax = self._make_cbar_axes()
        cb = self.figure.colorbar(sm, cax=cax)
        cb.set_label(r'$x = \sqrt{\phi / \phi_{\mathrm{lim}}}$',
                     fontsize=self.label_fontsize)
        cb.ax.tick_params(labelsize=self.tick_fontsize)

    def _overlay_lim(self, ax):
        if self._run is None or 'LIM' not in self._run.get('ufiles', {}):
            return
        uf = self._run['ufiles']['LIM']
        F = uf.get('F', {}).get('data')
        x0 = uf.get('X0', {}).get('data')
        if F is None or x0 is None or F.ndim != 1 or len(x0) != len(F):
            return
        x0u = (uf.get('X0', {}).get('units') or '').strip().upper()
        fu = (uf.get('F', {}).get('units') or '').strip().upper()
        if x0u not in {'M', 'CM', 'MM'} or fu != x0u:
            return
        R = np.append(np.asarray(x0, dtype=float), float(x0[0]))
        Z = np.append(np.asarray(F, dtype=float), float(F[0]))
        ax.plot(R, Z, '-', color='gray', lw=1.0, alpha=0.6,
                label='LIM', zorder=1)

    # ============================================================
    # Options dialog
    # ============================================================

    def _show_style_dialog(self):
        W = 150
        dlg = QDialog(self.frame)
        dlg.setWindowTitle("Plot Options")
        dlg.setMinimumWidth(300)
        dl = QVBoxLayout(dlg)

        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Color"))
        color_combo = QComboBox()
        color_combo.setFixedWidth(W)
        color_combo.addItems([
            "Gradient(viridis)", "Gradient(hot)", "Gradient(jet)",
            "Gradient(coolwarm)", "Gradient(plasma)", "Gradient(inferno)",
            "Fixed(tab10)", "Fixed(tab20)", "Fixed(Set1)", "Fixed(Set2)",
        ])
        color_combo.setCurrentText(self._color_mode)
        color_row.addWidget(color_combo)
        dl.addLayout(color_row)

        # Max legend items before ellipsis
        max_row = QHBoxLayout()
        max_row.addWidget(QLabel("Max legend items"))
        max_spin = QSpinBox(); max_spin.setFixedWidth(W)
        max_spin.setRange(4, 200)
        max_spin.setValue(int(self._legend_max_items))
        max_spin.setProperty('attr', '_legend_max_items')
        max_spin.setProperty('default', 20)
        max_row.addWidget(max_spin)
        dl.addLayout(max_row)

        for name, attr, lo, hi, default in [
            ("Label font size", 'label_fontsize', 6, 24, 12),
            ("Legend font size", 'legend_fontsize', 4, 20, 8),
            ("Tick font size", 'tick_fontsize', 6, 20, 10),
        ]:
            row = QHBoxLayout()
            row.addWidget(QLabel(name))
            spin = QSpinBox(); spin.setFixedWidth(W)
            spin.setRange(lo, hi); spin.setValue(getattr(self, attr))
            spin.setProperty('attr', attr); spin.setProperty('default', default)
            row.addWidget(spin); dl.addLayout(row)

        btns = QDialogButtonBox(
            QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Ok
            | QDialogButtonBox.Cancel)

        def _apply():
            self._color_mode = color_combo.currentText()
            for spin in dlg.findChildren(QSpinBox):
                setattr(self, spin.property('attr'), spin.value())
            self._plot()

        def _reset():
            color_combo.setCurrentText("Gradient(viridis)")
            for spin in dlg.findChildren(QSpinBox):
                spin.setValue(spin.property('default'))

        btns.accepted.connect(lambda: (_apply(), dlg.accept()))
        btns.rejected.connect(dlg.reject)
        btns.button(QDialogButtonBox.RestoreDefaults).clicked.connect(_reset)
        dl.addWidget(btns)
        dlg.show()

    # ============================================================
    # Settings persistence
    # ============================================================

    def _restore_settings(self):
        from config.user_settings import get_tab_settings
        s = get_tab_settings(self._settings_key)
        self._color_mode = s.get("color_mode", "Gradient(viridis)")
        self.label_fontsize = s.get("label_fontsize", 12)
        self.legend_fontsize = s.get("legend_fontsize", 8)
        self.tick_fontsize = s.get("tick_fontsize", 10)
        self._view_3d = bool(s.get("view_3d", False))
        self._legend_max_items = int(s.get("legend_max_items", 20))
        self._last_dir = s.get("last_dir", None)
        saved_var = s.get("variable", "")
        if saved_var:
            self._pending_variable = saved_var

    def save_settings(self):
        from config.user_settings import get_tab_settings, set_tab_settings
        s = get_tab_settings(self._settings_key)
        s["color_mode"] = self._color_mode
        s["label_fontsize"] = self.label_fontsize
        s["legend_fontsize"] = self.legend_fontsize
        s["tick_fontsize"] = self.tick_fontsize
        s["view_3d"] = bool(self._view_3d)
        s["legend_max_items"] = int(self._legend_max_items)
        if self._last_dir:
            s["last_dir"] = self._last_dir
        s["variable"] = self.var_combo.currentData() or ""
        set_tab_settings(self._settings_key, s)
