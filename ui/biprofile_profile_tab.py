"""
BiProfile Profile tab - browse Bayesian inference fitted profiles with time slider.
Overlays raw CES/Thomson data with USE-based coloring. Matches PRISM v2.3.4 styling.
"""

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QScrollArea, QLayout,
    QGroupBox, QGridLayout, QLabel, QLineEdit, QPushButton, QFrame,
    QSlider, QApplication, QMessageBox, QStyle, QSpinBox,
    QDialog, QDialogButtonBox, QColorDialog,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style
from ui.theme import ThemeManager

UNITS = {'Ti': 'keV', 'vT': 'km/s', 'Te': 'keV', 'ne': r'$10^{19}$/m$^3$'}
UNITS_QT = {'Ti': 'keV', 'vT': 'km/s', 'Te': 'keV', 'ne': '1e19/m3'}
LABELS = {'Ti': r'T$_i$', 'vT': r'v$_T$', 'Te': r'T$_e$', 'ne': r'n$_e$'}
_ZERO_BOTTOM = {'Ti', 'Te', 'ne'}

# Default color per tab: same color for both params in a tab
_DEFAULT_COLOR = '#1f77b4'


class BiProfileTab:
    """Profile tab with PRISM-matching styling. One color per tab."""

    def __init__(self, main_window, params):
        self.main = main_window
        self.params = params
        self._current_shot = None
        self._current_data = None
        self._axes = []
        self._artists = []
        self._legend_entries = {}  # {col: [(handle, label), ...]}

        self.plot_color = _DEFAULT_COLOR
        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10
        self._fixed_ylim = {}  # {col: (ymin, ymax)} computed on fetch
        self._fix_ylim = False  # toggle: True=fixed, False=dynamic

        self.frame = QWidget()
        self.canvas = None
        self._create_widgets()

    def _create_widgets(self):
        self.figure = Figure((10, 6))
        apply_dark_figure_style(self.figure)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.draw()

        canvas_widget = QWidget()
        cl = QVBoxLayout(canvas_widget)
        cl.setContentsMargins(0, 0, 0, 0)
        cl.addWidget(self.canvas)

        scroll_area = QScrollArea()
        scroll_area.setFixedWidth(CONTROL_PANEL_WIDTH)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        control_frame = QWidget()
        control_layout = QVBoxLayout(control_frame)
        control_layout.setContentsMargins(9, 9, 9, 9)
        control_layout.setSizeConstraint(QLayout.SetMinimumSize)
        scroll_area.setWidget(control_frame)
        scroll_area.viewport().setAutoFillBackground(False)
        control_frame.setAutoFillBackground(False)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(canvas_widget)
        splitter.addWidget(scroll_area)

        main_layout = QVBoxLayout(self.frame)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(splitter)

        self._create_shot_input(control_frame)
        self._create_plot_controls(control_frame)
        self._create_time_browser(control_frame)
        control_layout.addStretch()

    def _create_shot_input(self, parent):
        group = QGroupBox(f"1. Load {', '.join(self.params)} Data")
        grid = QGridLayout(group)
        grid.setColumnStretch(1, 1)
        grid.addWidget(QLabel("Shot"), 0, 0)
        self.shot_entry = QLineEdit()
        self.shot_entry.returnPressed.connect(self._fetch)
        grid.addWidget(self.shot_entry, 0, 1)
        btn_updown = QWidget()
        bl = QVBoxLayout(btn_updown)
        bl.setContentsMargins(0, 0, 0, 0); bl.setSpacing(0)
        sty = "padding: 0px; border-radius: 2px;"
        up = QPushButton()
        up.setIcon(QApplication.style().standardIcon(QStyle.SP_ArrowUp))
        up.setFixedSize(24, 15); up.setStyleSheet(sty)
        up.clicked.connect(lambda: self._adjust_shot(1)); bl.addWidget(up)
        dn = QPushButton()
        dn.setIcon(QApplication.style().standardIcon(QStyle.SP_ArrowDown))
        dn.setFixedSize(24, 15); dn.setStyleSheet(sty)
        dn.clicked.connect(lambda: self._adjust_shot(-1)); bl.addWidget(dn)
        grid.addWidget(btn_updown, 0, 2)
        self.fetch_button = QPushButton("Fetch")
        self.fetch_button.setFixedWidth(70)
        self.fetch_button.clicked.connect(self._fetch)
        grid.addWidget(self.fetch_button, 0, 3)
        parent.layout().addWidget(group)

    def _create_plot_controls(self, parent):
        from ui.widgets.toggle_switch import ToggleSwitch

        group = QGroupBox("2. Plot")
        layout = QVBoxLayout(group)

        # Status label (above Plot button)
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        # Plot + Option buttons
        btn_row = QHBoxLayout()
        plot_btn = QPushButton("Plot")
        plot_btn.clicked.connect(self._force_replot)
        btn_row.addWidget(plot_btn, 3)
        opt_btn = QPushButton("Option")
        opt_btn.clicked.connect(self._show_style_dialog)
        btn_row.addWidget(opt_btn, 1)
        layout.addLayout(btn_row)

        sep = QFrame(); sep.setFrameShape(QFrame.HLine); sep.setFrameShadow(QFrame.Sunken)
        layout.addWidget(sep)

        # Show Nodes toggle
        nodes_row = QHBoxLayout()
        self.show_nodes_toggle = ToggleSwitch()
        self.show_nodes_toggle.toggled.connect(self._on_show_nodes_toggled)
        nodes_row.addWidget(self.show_nodes_toggle)
        nodes_row.addWidget(QLabel("Show Nodes"))
        nodes_row.addStretch()
        layout.addLayout(nodes_row)

        # Apply Scale toggle (ne, Te tab only — BITS SCALE for Thomson ne)
        self._apply_scale = False
        if 'ne' in self.params or 'Te' in self.params:
            scale_row = QHBoxLayout()
            self.scale_toggle = ToggleSwitch()
            self.scale_toggle.toggled.connect(self._on_scale_toggled)
            scale_row.addWidget(self.scale_toggle)
            scale_row.addWidget(QLabel("Apply TS Scale"))
            scale_row.addStretch()
            layout.addLayout(scale_row)

        sep2 = QFrame(); sep2.setFrameShape(QFrame.HLine); sep2.setFrameShadow(QFrame.Sunken)
        layout.addWidget(sep2)

        self.info_label = QLabel("")
        self.info_label.setWordWrap(True)
        self.info_label.setStyleSheet("color: #888; font-size: 11px;")
        layout.addWidget(self.info_label)
        parent.layout().addWidget(group)

    def _create_time_browser(self, parent):
        from ui.widgets.toggle_switch import ToggleSwitch

        group = QGroupBox("3. Browse Data")
        layout = QVBoxLayout(group)
        self.time_label = QLabel("Time: --")
        self.time_label.setStyleSheet("font-weight: bold; font-size: 13px;")
        layout.addWidget(self.time_label)
        self.time_slider = QSlider(Qt.Horizontal)
        self.time_slider.setEnabled(False)
        self.time_slider.valueChanged.connect(self._on_slider)
        self.time_slider.sliderReleased.connect(self._on_slider_released)
        layout.addWidget(self.time_slider)
        self.range_label = QLabel("")
        self.range_label.setStyleSheet("color: #888; font-size: 11px;")
        layout.addWidget(self.range_label)

        # Fix Y-axis toggle
        fix_row = QHBoxLayout()
        self.fix_ylim_toggle = ToggleSwitch(checked=self._fix_ylim)
        self.fix_ylim_toggle.toggled.connect(self._on_fix_ylim_toggled)
        fix_row.addWidget(self.fix_ylim_toggle)
        fix_row.addWidget(QLabel("Fix Y-axis"))
        fix_row.addStretch()
        layout.addLayout(fix_row)

        # ymin/ymax grid: label | entry | ~ | entry
        ylim_grid = QGridLayout()
        ylim_grid.setSpacing(4)
        p1, p2 = self.params[0], self.params[1]
        u1, u2 = UNITS_QT[p1], UNITS_QT[p2]

        ylim_grid.addWidget(QLabel(f"{p1} [{u1}]"), 0, 0)
        self.ymin_entry_0 = QLineEdit()
        self.ymin_entry_0.returnPressed.connect(self._apply_custom_ylim)
        ylim_grid.addWidget(self.ymin_entry_0, 0, 1)
        tilde0 = QLabel("~"); tilde0.setAlignment(Qt.AlignCenter); tilde0.setFixedWidth(12)
        ylim_grid.addWidget(tilde0, 0, 2)
        self.ymax_entry_0 = QLineEdit()
        self.ymax_entry_0.returnPressed.connect(self._apply_custom_ylim)
        ylim_grid.addWidget(self.ymax_entry_0, 0, 3)

        ylim_grid.addWidget(QLabel(f"{p2} [{u2}]"), 1, 0)
        self.ymin_entry_1 = QLineEdit()
        self.ymin_entry_1.returnPressed.connect(self._apply_custom_ylim)
        ylim_grid.addWidget(self.ymin_entry_1, 1, 1)
        tilde1 = QLabel("~"); tilde1.setAlignment(Qt.AlignCenter); tilde1.setFixedWidth(12)
        ylim_grid.addWidget(tilde1, 1, 2)
        self.ymax_entry_1 = QLineEdit()
        self.ymax_entry_1.returnPressed.connect(self._apply_custom_ylim)
        ylim_grid.addWidget(self.ymax_entry_1, 1, 3)

        ylim_grid.setColumnStretch(1, 1)
        ylim_grid.setColumnStretch(3, 1)

        layout.addLayout(ylim_grid)

        self.apply_ylim_btn = QPushButton("Apply Y-axis")
        self.apply_ylim_btn.clicked.connect(self._apply_custom_ylim)
        layout.addWidget(self.apply_ylim_btn)

        self._ylim_entries = [
            (self.ymin_entry_0, self.ymax_entry_0),
            (self.ymin_entry_1, self.ymax_entry_1),
        ]
        self._set_ylim_entries_enabled(self._fix_ylim)

        self.browse_status = QLabel("")
        self.browse_status.setWordWrap(True)
        self.browse_status.setStyleSheet("color: #888; font-size: 11px;")
        layout.addWidget(self.browse_status)

        self._browse_group = group
        self._browse_group.setEnabled(False)
        parent.layout().addWidget(group)

    def _set_ylim_entries_enabled(self, enabled):
        for ymin_e, ymax_e in self._ylim_entries:
            ymin_e.setEnabled(enabled)
            ymax_e.setEnabled(enabled)
        if hasattr(self, 'apply_ylim_btn'):
            self.apply_ylim_btn.setEnabled(enabled)

    def _on_fix_ylim_toggled(self, checked):
        self._fix_ylim = checked
        self._set_ylim_entries_enabled(checked)
        if self._current_data:
            if checked:
                # Auto-fill if empty, and sync entry values to _fixed_ylim
                for col, (ymin_e, ymax_e) in enumerate(self._ylim_entries):
                    if not ymin_e.text() and not ymax_e.text() and col in self._fixed_ylim:
                        lo, hi = self._fixed_ylim[col]
                        ymin_e.setText(f"{lo:.4g}")
                        ymax_e.setText(f"{hi:.4g}")
                self._apply_custom_ylim()
            else:
                self._do_slider_update()

    def _apply_custom_ylim(self):
        """Apply user-entered ymin/ymax per param"""
        for col, (ymin_e, ymax_e) in enumerate(self._ylim_entries):
            try:
                lo = float(ymin_e.text())
                hi = float(ymax_e.text())
            except ValueError:
                continue
            if lo < hi:
                self._fixed_ylim[col] = (lo, hi)
        if self._current_data:
            self._do_slider_update()

    def _on_show_nodes_toggled(self, checked):
        if self._current_data is None or not self._plotted:
            return
        self._do_slider_update()

    def _on_scale_toggled(self, checked):
        self._apply_scale = checked
        if self._current_data is None or not self._plotted:
            return
        self._do_slider_update()

    # ---- Style Dialog ----

    def _show_style_dialog(self):
        W = 150
        dlg = QDialog(self.frame)
        dlg.setWindowTitle("Plot Options")
        dlg.setMinimumWidth(300)
        dl = QVBoxLayout(dlg)

        # Color picker
        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Color"))
        self._color_btn = QPushButton()
        self._color_btn.setFixedSize(W, 24)
        self._dialog_color = self.plot_color
        self._color_btn.setStyleSheet(
            f"background-color: {self._dialog_color}; border: 1px solid #555; border-radius: 4px;")
        self._color_btn.clicked.connect(self._pick_color)
        color_row.addWidget(self._color_btn)
        dl.addLayout(color_row)

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
            QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept); btns.rejected.connect(dlg.reject)
        def _reset():
            self._dialog_color = _DEFAULT_COLOR
            self._color_btn.setStyleSheet(
                f"background-color: {self._dialog_color}; border: 1px solid #555; border-radius: 4px;")
            for spin in dlg.findChildren(QSpinBox):
                spin.setValue(spin.property('default'))
        btns.button(QDialogButtonBox.RestoreDefaults).clicked.connect(_reset)
        dl.addWidget(btns)

        if dlg.exec() == QDialog.Accepted:
            self.plot_color = self._dialog_color
            for spin in dlg.findChildren(QSpinBox):
                setattr(self, spin.property('attr'), spin.value())
            if self._current_data:
                self._force_replot()

    def _pick_color(self):
        c = QColorDialog.getColor(QColor(self._dialog_color), self.frame, "Plot Color")
        if c.isValid():
            self._dialog_color = c.name()
            self._color_btn.setStyleSheet(
                f"background-color: {self._dialog_color}; border: 1px solid #555; border-radius: 4px;")

    # ---- Actions ----

    def _update_status(self, message, color='blue'):
        self.status_label.setStyleSheet(f"color: {color}; font-weight: bold; font-size: 9pt;")
        self.status_label.setText(message)
        QApplication.processEvents()

    def _adjust_shot(self, delta):
        try:
            self.shot_entry.setText(str(max(1, int(self.shot_entry.text()) + delta)))
        except ValueError:
            pass

    def _fetch(self):
        shot_text = self.shot_entry.text().strip()
        if not shot_text:
            return
        try:
            shot = int(shot_text)
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
            return
        self._update_status(f"Loading #{shot}...", color='blue')
        self.fetch_button.setText("Loading..."); self.fetch_button.setEnabled(False)
        QApplication.setOverrideCursor(Qt.WaitCursor); QApplication.processEvents()
        try:
            sdata = self.main.fetch_biprofile_shot(shot)
        finally:
            self.fetch_button.setText("Fetch"); self.fetch_button.setEnabled(True)
            QApplication.restoreOverrideCursor()
        if not sdata:
            self._update_status(f"#{shot}: No data", color='red')
            return
        bi = sdata['bi']
        ref = next((p for p in self.params if p in bi), None)
        if ref is None:
            self._update_status(f"#{shot}: No data for {', '.join(self.params)}", color='red')
            return

        self._current_shot = shot
        self._current_data = sdata

        # Precompute fixed ylim from all time slices
        self._fixed_ylim = {}
        for col, param in enumerate(self.params):
            if param not in bi:
                continue
            d = bi[param]
            all_mean = d['mean']
            all_unc = d['unc']
            lo = np.nanmin(all_mean - all_unc)
            hi = np.nanmax(all_mean + all_unc)
            margin = (hi - lo) * 0.1 if hi > lo else 0.1
            if param in _ZERO_BOTTOM:
                self._fixed_ylim[col] = (0, hi + margin)
            else:
                self._fixed_ylim[col] = (lo - margin, hi + margin)

        time_arr = bi[ref]['time']
        self.time_slider.blockSignals(True)
        self.time_slider.setMinimum(0); self.time_slider.setMaximum(len(time_arr) - 1)
        self.time_slider.setValue(0); self.time_slider.setEnabled(True)
        self.time_slider.blockSignals(False)
        self._plotted = False
        self._browse_group.setEnabled(False)
        loaded = [p for p in self.params if p in bi]
        self._update_status(f"#{shot}: {', '.join(loaded)} loaded  \u2192 Press Plot", color='green')
        self.info_label.setText(
            f"EFIT: {bi[ref].get('efit_used','?')}\n"
            f"Fit: {bi[ref].get('fit_func','?')}\n"
            f"{len(time_arr)} time pts, {len(bi[ref]['psin'])} \u03c8\u2099 pts")
        self.range_label.setText(f"{time_arr[0]:.4f} ~ {time_arr[-1]:.4f} s")

    def _force_replot(self):
        if self._current_data:
            self._init_axes()
            self._update_plot(self.time_slider.value())
            self._plotted = True
            self._browse_group.setEnabled(True)
            self._set_ylim_entries_enabled(self._fix_ylim)
            self._update_status(
                f"#{self._current_shot}: {', '.join(p for p in self.params if p in self._current_data['bi'])} plotted",
                color='green')
            self.browse_status.setText("Drag slider to browse time points")

    def _on_slider(self, value):
        """During drag: lightweight update (fit line only, no raw overlay)"""
        if self._current_data is None or not self._plotted:
            return
        self._update_plot(value, lightweight=True)
        bi = self._current_data['bi']
        ref = next((p for p in self.params if p in bi), None)
        t = bi[ref]['time'][value]
        ms = int(round(t * 1000))
        self.browse_status.setText(f"Dragging... t={t:.4f}s ({ms:06d}ms)")

    def _on_slider_released(self):
        """On release: full update with raw overlay"""
        if self._current_data is None or not self._plotted:
            return
        self._update_plot(self.time_slider.value(), lightweight=False)
        bi = self._current_data['bi']
        ref = next((p for p in self.params if p in bi), None)
        t = bi[ref]['time'][self.time_slider.value()]
        ms = int(round(t * 1000))
        self.browse_status.setText(f"t={t:.4f}s ({ms:06d}ms) — raw overlay updated")

    def _do_slider_update(self):
        """Called by toggles/apply to refresh at current slider position"""
        if not self._plotted:
            return
        self._update_plot(self.time_slider.value(), lightweight=False)

    # ---- Plotting ----

    def _init_axes(self):
        self.figure.clear()
        self._axes = []
        self._artists = []
        self._legend_entries = {}

        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'

        for col, param in enumerate(self.params):
            ax = self.figure.add_subplot(1, 2, col + 1)
            ax.set_xlabel(r'$\psi_N$', fontsize=self.label_fontsize)
            ax.set_ylabel(f'{LABELS.get(param, param)} [{UNITS.get(param, "")}]', fontsize=self.label_fontsize)
            ax.tick_params(labelsize=self.tick_fontsize)
            ax.set_xlim(0, 1.1)
            ax.grid(ls='--', lw=0.3, color='#444444')
            ax.axvline(x=1.0, color=zc, ls='--', lw=0.8, alpha=0.5, gid='zero_ref')
            self._axes.append(ax)
            self._legend_entries[col] = []

        self.figure.subplots_adjust(left=0.10, right=0.97, top=0.92, bottom=0.10, wspace=0.20)
        apply_dark_figure_style(self.figure)

    def _update_plot(self, t_idx, lightweight=False):
        """Update plot. lightweight=True skips raw overlay for fast slider dragging."""
        sdata = self._current_data
        shot = self._current_shot
        bi = sdata['bi']
        ref = next((p for p in self.params if p in bi), None)
        t_actual = bi[ref]['time'][t_idx]
        ms = int(round(t_actual * 1000))

        self.time_label.setText(f"Time: {t_actual:.4f} s  ({t_idx}/{len(bi[ref]['time'])-1})")

        for a in self._artists:
            a.remove()
        self._artists = []

        for ax in self._axes:
            leg = ax.get_legend()
            if leg:
                leg.remove()

        color = self.plot_color
        self.figure.suptitle(f'#{shot}  {ms:06d}ms', fontsize=self.label_fontsize)

        for col, param in enumerate(self.params):
            ax = self._axes[col]
            self._legend_entries[col] = []

            if param not in bi:
                continue

            d = bi[param]
            psin, mean_p, unc_p = d['psin'], d['mean'][:, t_idx], d['unc'][:, t_idx]
            valid = ~np.isnan(mean_p)

            if np.any(valid):
                line, = ax.plot(psin[valid], mean_p[valid], '-', color=color, lw=2)
                self._artists.append(line)
                self._legend_entries[col].append((line, 'fit'))

                fill = ax.fill_between(psin[valid],
                                       mean_p[valid] - unc_p[valid],
                                       mean_p[valid] + unc_p[valid],
                                       alpha=0.2, color=color)
                self._artists.append(fill)

                if self._fix_ylim and col in self._fixed_ylim:
                    ax.set_ylim(self._fixed_ylim[col])
                else:
                    y_min = np.nanmin(mean_p[valid] - unc_p[valid])
                    y_max = np.nanmax(mean_p[valid] + unc_p[valid])
                    margin = (y_max - y_min) * 0.1 if y_max > y_min else 0.1
                    ax.set_ylim(0 if param in _ZERO_BOTTOM else y_min - margin,
                                y_max + margin)

            # Raw overlay only on full update (not during drag)
            if not lightweight:
                self._overlay_raw(ax, col, param, t_actual, color)

            entries = self._legend_entries[col]
            if entries:
                ax.legend([h for h, _ in entries], [l for _, l in entries],
                          fontsize=self.legend_fontsize, loc='best', frameon=False)

        self.canvas.draw_idle()

    # ---- Raw overlay ----

    def _overlay_raw(self, ax, col, param, t_actual, color):
        from data_loaders.biprofile_loader import map_R_to_psin
        sdata = self._current_data
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
            ch_labels = [f'CES{i+1:02d}' for i in range(len(R))]
            self._plot_raw_with_use(ax, col, psin_mapped, y, yerr, use_flags,
                                    color, 'o', 'CES', ch_labels)
        elif param in ['Te', 'ne']:
            raw = sdata.get('thomson')
            if raw is None:
                return
            t_idx = np.argmin(np.abs(raw['time'] - t_actual))
            R = raw['R']
            psin_mapped = map_R_to_psin(R, efit, t_actual)
            y = raw['Te'][:, t_idx] if param == 'Te' else raw['ne'][:, t_idx]
            yerr = raw['Te_err'][:, t_idx] if param == 'Te' else raw['ne_err'][:, t_idx]

            # Apply BITS SCALE to ne if toggle is on
            if param == 'ne' and self._apply_scale and diag:
                scale = self._get_ts_scale(diag, t_actual, len(R))
                y = y * scale
                yerr = yerr * scale

            use_flags = self._get_use_flags(diag, param, len(R))
            n_core = sum(1 for ch in (diag or {}).get('thomson', {}).get('core', [])
                         if ch.get('ne_use') is not None or ch.get('te_use') is not None)
            ch_labels = []
            for i in range(len(R)):
                if i < n_core:
                    ch_labels.append(f'TS_C{i+1}')
                else:
                    ch_labels.append(f'TS_E{i-n_core+1}')
            label = 'Thomson (scaled)' if (param == 'ne' and self._apply_scale) else 'Thomson'
            self._plot_raw_with_use(ax, col, psin_mapped, y, yerr, use_flags,
                                    color, 's', label, ch_labels)

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
        """Get per-channel BITS SCALE values for Thomson ne at given time.

        SCALE is time-dependent (shape=(252,) per channel).
        Thomson loader sorts channels by R (core first, then edge).
        """
        scales = np.ones(n_ch)
        ts = diag.get('thomson', {})

        # Time index in BINE grid
        bi = self._current_data['bi']
        ref = next((p for p in ('ne', 'Te') if p in bi), None)
        if ref is None:
            return scales
        t_idx = np.argmin(np.abs(bi[ref]['time'] - t_actual))

        idx = 0
        for region in ['core', 'edge']:
            for ch in ts.get(region, []):
                if idx >= n_ch:
                    break
                scale_arr = ch.get('scale')
                if scale_arr is not None:
                    if isinstance(scale_arr, np.ndarray) and len(scale_arr) > 0:
                        scales[idx] = scale_arr[min(t_idx, len(scale_arr) - 1)]
                    elif not isinstance(scale_arr, np.ndarray):
                        scales[idx] = float(scale_arr)
                idx += 1
        return scales

    def _plot_raw_with_use(self, ax, col, psin, y, yerr, use_flags, color,
                          marker, diag_name, ch_labels=None):
        yerr = np.where(yerr < 0, 0, yerr)
        valid = (psin >= 0) & (psin <= 1.2) & np.isfinite(y)
        used = valid & (use_flags == 1)
        excluded = valid & (use_flags == 0)
        unknown = valid & np.isnan(use_flags)

        labeled = False
        if np.any(used):
            c = ax.errorbar(psin[used], y[used], yerr=yerr[used],
                            fmt=marker, markersize=5, color=color,
                            capsize=5, zorder=10, markeredgewidth=1)
            self._artists.append(c)
            self._legend_entries[col].append((c, diag_name))
            labeled = True
        if np.any(excluded):
            c = ax.errorbar(psin[excluded], y[excluded], yerr=yerr[excluded],
                            fmt=marker, markersize=5,
                            color=(0.6, 0.6, 0.6, 0.35),
                            capsize=5, zorder=5, markeredgewidth=1)
            self._artists.append(c)
        if np.any(unknown):
            c = ax.errorbar(psin[unknown], y[unknown], yerr=yerr[unknown],
                            fmt=marker, markersize=5, color=color,
                            capsize=5, zorder=10, markeredgewidth=1)
            self._artists.append(c)
            if not labeled:
                self._legend_entries[col].append((c, diag_name))

        # Show Nodes: annotate channel labels above data points
        if ch_labels and self.show_nodes_toggle.isChecked():
            show_mask = valid  # annotate all valid points (used + excluded + unknown)
            indices = np.where(show_mask)[0]
            for idx in indices:
                ann = ax.annotate(
                    ch_labels[idx], (psin[idx], y[idx]),
                    textcoords='offset points', xytext=(0, 5),
                    ha='center', fontsize=7, alpha=0.8,
                    clip_on=True, annotation_clip=True)
                self._artists.append(ann)
