"""
BiProfile Time Trace tab - all traces in gray (LineCollection), selected psi_N highlighted.
Matches PRISM v2.3.4 styling. One color per tab, configurable via Option dialog.
"""

import numpy as np
from matplotlib.collections import LineCollection
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
_DEFAULT_COLOR = '#1f77b4'


class BiTimeTraceTab:
    """Time trace tab with PRISM-matching styling. One color per tab."""

    def __init__(self, main_window, params):
        self.main = main_window
        self.params = params
        self._current_shot = None
        self._current_data = None
        self._gray_drawn = False
        self._highlight_artists = {}
        self._axes = []

        self.plot_color = _DEFAULT_COLOR
        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10

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
        self._create_psin_browser(control_frame)
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
        group = QGroupBox("2. Plot")
        layout = QVBoxLayout(group)

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

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
        self.info_label = QLabel("")
        self.info_label.setWordWrap(True)
        self.info_label.setStyleSheet("color: #888; font-size: 11px;")
        layout.addWidget(self.info_label)
        parent.layout().addWidget(group)

    def _create_psin_browser(self, parent):
        from ui.widgets.toggle_switch import ToggleSwitch

        group = QGroupBox("3. Browse Data")
        layout = QVBoxLayout(group)
        self.psin_label = QLabel(u"\u03c8\u2099: --")
        self.psin_label.setStyleSheet("font-weight: bold; font-size: 13px;")
        layout.addWidget(self.psin_label)
        self.psin_slider = QSlider(Qt.Horizontal)
        self.psin_slider.setEnabled(False)
        self.psin_slider.valueChanged.connect(self._on_slider)
        layout.addWidget(self.psin_slider)
        self.range_label = QLabel("")
        self.range_label.setStyleSheet("color: #888; font-size: 11px;")
        layout.addWidget(self.range_label)

        # Fix Y-axis toggle
        fix_row = QHBoxLayout()
        self.fix_ylim_toggle = ToggleSwitch(checked=False)
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
        self._fix_ylim = False
        self._fixed_ylim = {}
        self._set_ylim_entries_enabled(False)

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
        if self._current_data and self._gray_drawn:
            if checked:
                for col, (ymin_e, ymax_e) in enumerate(self._ylim_entries):
                    if not ymin_e.text() and not ymax_e.text() and col in self._fixed_ylim:
                        lo, hi = self._fixed_ylim[col]
                        ymin_e.setText(f"{lo:.4g}")
                        ymax_e.setText(f"{hi:.4g}")
                self._apply_custom_ylim()
            else:
                self._force_replot()

    def _apply_custom_ylim(self):
        for col, (ymin_e, ymax_e) in enumerate(self._ylim_entries):
            try:
                lo = float(ymin_e.text())
                hi = float(ymax_e.text())
            except ValueError:
                continue
            if lo < hi:
                self._fixed_ylim[col] = (lo, hi)
        if self._current_data and self._gray_drawn:
            self._force_replot()

    # ---- Style Dialog ----

    def _show_style_dialog(self):
        W = 150
        dlg = QDialog(self.frame)
        dlg.setWindowTitle("Plot Options")
        dlg.setMinimumWidth(300)
        dl = QVBoxLayout(dlg)

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
        self._gray_drawn = False
        self._fixed_ylim = {}  # reset for new shot
        psin = bi[ref]['psin']
        n_psi = len(psin)
        initial_idx = n_psi // 2
        self.psin_slider.blockSignals(True)
        self.psin_slider.setMinimum(0); self.psin_slider.setMaximum(n_psi - 1)
        self.psin_slider.setValue(initial_idx); self.psin_slider.setEnabled(True)
        self.psin_slider.blockSignals(False)
        loaded = [p for p in self.params if p in bi]
        n_time = len(bi[ref]['time'])
        self._browse_group.setEnabled(False)
        self._update_status(f"#{shot}: {', '.join(loaded)} loaded  \u2192 Press Plot", color='green')
        self.info_label.setText(
            f"{n_psi} \u03c8\u2099 points, {n_time} time points\n"
            f"Time: {bi[ref]['time'][0]:.3f} ~ {bi[ref]['time'][-1]:.3f} s")
        self.range_label.setText(f"\u03c8\u2099: {psin[0]:.4f} ~ {psin[-1]:.4f}")

    def _force_replot(self):
        if self._current_data:
            self._draw_all(self.psin_slider.value())
            self._browse_group.setEnabled(True)
            self._set_ylim_entries_enabled(self._fix_ylim)
            self._update_status(
                f"#{self._current_shot}: {', '.join(p for p in self.params if p in self._current_data['bi'])} plotted",
                color='green')
            self.browse_status.setText(u"Drag slider to browse \u03c8\u2099 positions")

    def _on_slider(self, value):
        if self._current_data is None or not self._gray_drawn:
            return
        self._update_psin_label(value)

        bi = self._current_data['bi']
        ref = next((p for p in self.params if p in bi), None)
        psin_val = bi[ref]['psin'][value]

        self.browse_status.setText(f"\u03c8\u2099 = {psin_val:.4f}  (idx {value})")

        self.figure.suptitle(
            f'#{self._current_shot}  \u03c8\u2099={psin_val:.3f}',
            fontsize=self.label_fontsize)

        # Swap highlight artists directly (no debounce — fast enough)
        for col, param in enumerate(self.params):
            if param not in bi or col >= len(self._axes):
                continue
            ax = self._axes[col]
            for a in self._highlight_artists.get(col, []):
                a.remove()
            self._highlight_artists[col] = []
            self._draw_highlight(ax, col, param, bi[param], value)
        self.canvas.draw_idle()

    # ---- Plotting ----

    def _draw_all(self, psi_idx):
        bi = self._current_data['bi']
        shot = self._current_shot

        self.figure.clear()
        self._axes = []
        self._highlight_artists = {}

        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'

        for col, param in enumerate(self.params):
            ax = self.figure.add_subplot(1, 2, col + 1)
            self._axes.append(ax)
            self._highlight_artists[col] = []

            if param not in bi:
                ax.set_xlabel('Time [s]', fontsize=self.label_fontsize)
                ax.set_ylabel(f'{LABELS.get(param, param)} [{UNITS.get(param, "")}]', fontsize=self.label_fontsize)
                continue

            d = bi[param]
            time, mean = d['time'], d['mean']

            segments = []
            for i in range(mean.shape[0]):
                valid = ~np.isnan(mean[i, :])
                if np.any(valid):
                    segments.append(np.column_stack([time[valid], mean[i, valid]]))
            if segments:
                ax.add_collection(LineCollection(
                    segments, colors='#555555', linewidths=0.4, alpha=0.5))

            ax.set_ylabel(f'{LABELS.get(param, param)} [{UNITS.get(param, "")}]', fontsize=self.label_fontsize)
            ax.set_xlabel('Time [s]', fontsize=self.label_fontsize)
            ax.tick_params(labelsize=self.tick_fontsize)
            ax.grid(ls='--', lw=0.3, color='#444444')

            if param == 'vT':
                ax.axhline(y=0, color=zc, ls='--', gid='zero_ref')

            # Compute default ylim from core psi_N range (skip boundary outliers)
            psin = d['psin']
            core_mask = (psin >= 0.05) & (psin <= 0.95)
            core_mean = mean[core_mask, :]
            core_valid = ~np.isnan(core_mean)
            if np.any(core_valid):
                vals = core_mean[core_valid]
                y_max = np.nanmax(vals)
                margin = y_max * 0.15 if y_max > 0 else 0.1
                if param in _ZERO_BOTTOM:
                    auto_ylim = (0, y_max + margin)
                else:
                    y_min = np.nanmin(vals)
                    margin_abs = (y_max - y_min) * 0.15 if y_max > y_min else 0.1
                    auto_ylim = (y_min - margin_abs, y_max + margin_abs)

                # Precompute for Fix Y-axis (first plot only)
                if col not in self._fixed_ylim:
                    self._fixed_ylim[col] = auto_ylim

                if self._fix_ylim and col in self._fixed_ylim:
                    ax.set_ylim(self._fixed_ylim[col])
                else:
                    ax.set_ylim(auto_ylim)
                ax.set_xlim(time[0], time[-1])

            self._draw_highlight(ax, col, param, d, psi_idx)

        psin_val = bi[next(p for p in self.params if p in bi)]['psin'][psi_idx]
        self.figure.suptitle(
            f'#{shot}  \u03c8\u2099={psin_val:.3f}', fontsize=self.label_fontsize)
        self.figure.subplots_adjust(left=0.10, right=0.97, top=0.92, bottom=0.10, wspace=0.20)
        apply_dark_figure_style(self.figure)
        self._update_psin_label(psi_idx)
        self.canvas.draw()
        self._gray_drawn = True

    def _draw_highlight(self, ax, col, param, d, psi_idx):
        time = d['time']
        mean, unc = d['mean'][psi_idx, :], d['unc'][psi_idx, :]
        valid = ~np.isnan(mean)
        color = self.plot_color

        artists = []
        if np.any(valid):
            # fill first (behind), then line on top — both above gray (zorder=10)
            fill = ax.fill_between(time[valid],
                                   mean[valid] - unc[valid],
                                   mean[valid] + unc[valid],
                                   alpha=0.25, color=color, zorder=10)
            artists.append(fill)
            line, = ax.plot(time[valid], mean[valid], '-', color=color, lw=2, zorder=11)
            artists.append(line)
        self._highlight_artists[col] = artists

    def _update_psin_label(self, psi_idx):
        bi = self._current_data['bi']
        ref = next((p for p in self.params if p in bi), None)
        psin = bi[ref]['psin']
        self.psin_label.setText(
            f"\u03c8\u2099: {psin[psi_idx]:.4f}  ({psi_idx}/{len(psin)-1})")
