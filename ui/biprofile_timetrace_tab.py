"""
BiProfile Time Trace tab - Available/Selected workflow with Preview dialog.
Plots multiple time traces (one per selected psi_N) with uncertainty bands.
"""

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QScrollArea, QLayout,
    QGroupBox, QGridLayout, QLabel, QLineEdit, QPushButton, QFrame,
    QListWidget, QAbstractItemView, QFileDialog,
    QApplication, QMessageBox, QStyle, QSpinBox,
    QDialog, QDialogButtonBox,
    QTableWidget, QTableWidgetItem, QHeaderView,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QShortcut, QKeySequence

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style, get_icon
from ui.theme import ThemeManager

UNITS = {'Ti': 'keV', 'vT': 'km/s', 'Te': 'keV', 'ne': r'$10^{19}$/m$^3$'}
UNITS_QT = {'Ti': 'keV', 'vT': 'km/s', 'Te': 'keV', 'ne': '1e19/m3'}
LABELS = {'Ti': r'T$_i$', 'vT': r'v$_T$', 'Te': r'T$_e$', 'ne': r'n$_e$'}
_ZERO_BOTTOM = {'Ti', 'Te', 'ne'}
_DEFAULT_COLOR = '#1f77b4'

COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']


class BiTimeTraceTab:
    """BiProfile TimeTrace tab: Available/Selected + Preview + multi-trace Plot."""

    def __init__(self, main_window, params):
        self.main = main_window
        self.params = params
        self._current_shot = None
        self._current_data = None

        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10

        param_key = ''.join(p.lower() for p in params)
        self._settings_key = f"bi_{param_key}_timetrace"

        self.frame = QWidget()
        self.canvas = None
        self._create_widgets()
        self._restore_settings()

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
        self._create_select_data(control_frame)
        self._create_plot_controls(control_frame)
        self._create_save_controls(control_frame)
        control_layout.addStretch()

    # ---- 1. Load Data ----

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

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        grid.addWidget(self.status_label, 1, 0, 1, 4)

        parent.layout().addWidget(group)

    # ---- 2. Select Data ----

    def _create_select_data(self, parent):
        group = QGroupBox("2. Select Data")
        layout = QVBoxLayout(group)

        self.browse_button = QPushButton("Fetch a shot to browse")
        self.browse_button.setEnabled(False)
        self.browse_button.clicked.connect(self._open_browse)
        layout.addWidget(self.browse_button)

        lists_row = QHBoxLayout()

        avail_col = QVBoxLayout()
        avail_col.addWidget(QLabel(u"\u03c8\u2099"))
        self.available_listbox = QListWidget()
        self.available_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.available_listbox.setFixedHeight(150)
        avail_col.addWidget(self.available_listbox)
        lists_row.addLayout(avail_col, stretch=1)

        btn_col = QVBoxLayout()
        btn_col.addStretch()
        add_btn = QPushButton()
        add_btn.setIcon(get_icon(QStyle.SP_ArrowForward))
        add_btn.setFixedWidth(30)
        add_btn.clicked.connect(self._add_items)
        btn_col.addWidget(add_btn)
        rm_btn = QPushButton()
        rm_btn.setIcon(get_icon(QStyle.SP_ArrowBack))
        rm_btn.setFixedWidth(30)
        rm_btn.clicked.connect(self._remove_items)
        btn_col.addWidget(rm_btn)
        btn_col.addStretch()
        lists_row.addLayout(btn_col)

        sel_col = QVBoxLayout()
        sel_col.addWidget(QLabel("Selected"))
        self.selected_listbox = QListWidget()
        self.selected_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.selected_listbox.setFixedHeight(150)
        sel_col.addWidget(self.selected_listbox)
        lists_row.addLayout(sel_col, stretch=1)

        layout.addLayout(lists_row)

        QShortcut(QKeySequence(Qt.Key_Delete), self.selected_listbox
                  ).activated.connect(self._remove_items)

        parent.layout().addWidget(group)

    # ---- 3. Plot ----

    def _create_plot_controls(self, parent):
        group = QGroupBox("3. Plot")
        layout = QVBoxLayout(group)

        btn_row = QHBoxLayout()
        plot_btn = QPushButton("Plot")
        plot_btn.clicked.connect(self._plot)
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

    # ---- Settings ----

    def _restore_settings(self):
        from config.user_settings import get_tab_settings
        s = get_tab_settings(self._settings_key)
        if s.get("shot"):
            self.shot_entry.setText(str(s["shot"]))
        self._color_mode = s.get("color_mode", "Gradient(viridis)")
        self.label_fontsize = s.get("label_fontsize", 12)
        self.legend_fontsize = s.get("legend_fontsize", 8)
        self.tick_fontsize = s.get("tick_fontsize", 10)
        if s.get("show_nodes") and hasattr(self, 'show_nodes_toggle'):
            self.show_nodes_toggle.setChecked(s["show_nodes"], animate=False)

    def save_settings(self):
        from config.user_settings import get_tab_settings, set_tab_settings
        s = get_tab_settings(self._settings_key)
        s["shot"] = self.shot_entry.text()
        s["color_mode"] = getattr(self, '_color_mode', 'Gradient(viridis)')
        s["label_fontsize"] = self.label_fontsize
        s["legend_fontsize"] = self.legend_fontsize
        s["tick_fontsize"] = self.tick_fontsize
        if hasattr(self, 'show_nodes_toggle'):
            s["show_nodes"] = self.show_nodes_toggle.isChecked()
        set_tab_settings(self._settings_key, s)

    # ---- Actions ----

    def _adjust_shot(self, delta):
        try:
            self.shot_entry.setText(str(max(1, int(self.shot_entry.text()) + delta)))
        except ValueError:
            pass

    def _update_status(self, message, color='blue'):
        self.status_label.setStyleSheet(f"color: {color}; font-weight: bold; font-size: 9pt;")
        self.status_label.setText(message)
        QApplication.processEvents()

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
            sdata = self.main.fetch_biprofile_shot(shot, params=self.params)
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

        # Populate available list: {shot}_{idx} (Bi) — idx is psi_N index
        self.available_listbox.clear()
        for i, psi in enumerate(bi[ref]['psin']):
            self.available_listbox.addItem(f"{shot:06d}_{i:03d} (Bi)")

        loaded = [p for p in self.params if p in bi]
        self._update_status(f"#{shot}: {', '.join(loaded)} loaded", color='green')
        self.info_label.setText(
            f"{len(bi[ref]['psin'])} \u03c8\u2099 pts, {len(bi[ref]['time'])} time pts\n"
            f"Time: {bi[ref]['time'][0]:.3f} ~ {bi[ref]['time'][-1]:.3f} s")

        self.browse_button.setEnabled(True)
        self.browse_button.setText(f"Browse #{shot}")

    def _open_browse(self):
        if not self._current_data:
            return
        from ui.widgets.preview_dialog import BiProfilePreviewDialog
        dlg = BiProfilePreviewDialog(
            self.frame, self._current_data['bi'], self._current_shot,
            self.params, mode='timetrace',
            selected_listbox=self.selected_listbox)
        dlg.show()

    def _add_items(self):
        existing = {self.selected_listbox.item(i).text()
                    for i in range(self.selected_listbox.count())}
        for item in self.available_listbox.selectedItems():
            if item.text() not in existing:
                self.selected_listbox.addItem(item.text())

    def _remove_items(self):
        for item in reversed(self.selected_listbox.selectedItems()):
            self.selected_listbox.takeItem(self.selected_listbox.row(item))

    # ---- Style Dialog ----

    # ---- 4. Save Data ----

    def _create_save_controls(self, parent):
        group = QGroupBox("4. Save Data")
        layout = QVBoxLayout(group)
        btn = QPushButton("Preview && Save")
        btn.clicked.connect(self._preview_and_save)
        layout.addWidget(btn)
        parent.layout().addWidget(group)

    def _preview_and_save(self):
        lines = self._build_csv_lines()
        if not lines:
            QMessageBox.warning(self.frame, "Warning", "No data to preview")
            return
        self._show_data_preview(lines)

    def _build_csv_lines(self):
        entries = [self.selected_listbox.item(i).text()
                   for i in range(self.selected_listbox.count())]
        if not entries or not self._current_data:
            return None

        bi = self._current_data['bi']
        shot = self._current_shot
        lines = []
        lines.append(f"# BiProfile Time Trace\n")
        header_parts = ["Shot", "psi_N", "Time[s]"]
        for p in self.params:
            units = UNITS.get(p, '')
            header_parts.append(f"{p}[{units}]")
            header_parts.append(f"{p}_unc")
        lines.append("#" + ','.join(f"{h:>12s}" for h in header_parts) + "\n")

        for entry in entries:
            s, psi_idx = self._parse_entry(entry)
            ref = next((p for p in self.params if p in bi), None)
            if ref is None:
                continue
            d_ref = bi[ref]
            psi_actual = d_ref['psin'][psi_idx] if psi_idx < len(d_ref['psin']) else 0
            time = d_ref['time']

            for t_idx, t in enumerate(time):
                row = f" {s:>12d},{psi_actual:>12.6f},{t:>12.4f}"
                for p in self.params:
                    if p in bi:
                        m = bi[p]['mean'][psi_idx, t_idx]
                        u = bi[p]['unc'][psi_idx, t_idx]
                        row += f",{m:>12.6e},{u:>12.6e}"
                    else:
                        row += f",{'':>12s},{'':>12s}"
                lines.append(row + "\n")

        return lines

    def _show_data_preview(self, lines):
        from PySide6.QtGui import QFont, QKeySequence

        headers = []
        data_rows = []
        title_line = ""
        last_comment_with_comma = None

        for line in lines:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                content = stripped.lstrip('#').strip()
                if ',' in content:
                    last_comment_with_comma = content
                else:
                    title_line = content
            else:
                if not headers and last_comment_with_comma:
                    headers = [h.strip() for h in last_comment_with_comma.split(',')]
                data_rows.append([c.strip() for c in stripped.split(',')])

        if not data_rows:
            QMessageBox.information(self.frame, "Preview", "No data to display")
            return

        n_cols = len(headers) if headers else len(data_rows[0])
        n_rows = len(data_rows)

        dlg = QDialog(self.frame)
        dlg.setWindowTitle(f"Data Preview — {title_line}" if title_line else "Data Preview")
        dlg.resize(min(120 * n_cols, 1200), min(30 * n_rows + 100, 700))
        dl = QVBoxLayout(dlg)

        info = QLabel(f"{n_rows} rows \u00d7 {n_cols} columns")
        info.setStyleSheet("color: gray;")
        dl.addWidget(info)

        table = QTableWidget(n_rows, n_cols)
        table.setFont(QFont('Courier', 9))
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        if headers:
            table.setHorizontalHeaderLabels(headers)
        table.horizontalHeader().setStretchLastSection(True)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        table.verticalHeader().setDefaultSectionSize(24)

        for r, row in enumerate(data_rows):
            for c, val in enumerate(row):
                if c < n_cols:
                    table.setItem(r, c, QTableWidgetItem(val))

        table.resizeColumnsToContents()
        for c in range(n_cols):
            if table.columnWidth(c) > 150:
                table.setColumnWidth(c, 150)

        dl.addWidget(table)

        def _copy_selection():
            sel = table.selectedRanges()
            if not sel:
                return
            rows = range(sel[0].topRow(), sel[0].bottomRow() + 1)
            cols = range(sel[0].leftColumn(), sel[0].rightColumn() + 1)
            text = '\n'.join(
                '\t'.join((table.item(r, c).text() if table.item(r, c) else '')
                          for c in cols) for r in rows)
            QApplication.clipboard().setText(text)
        QShortcut(QKeySequence.Copy, table).activated.connect(_copy_selection)

        btn_row = QHBoxLayout()
        btn_row.addStretch()

        save_btn = QPushButton("Save as .csv")
        def _save():
            import os
            path, _ = QFileDialog.getSaveFileName(
                dlg, "Save Data", os.path.expanduser("~"),
                "CSV files (*.csv);;All files (*)")
            if not path:
                return
            if not os.path.splitext(path)[1]:
                path += '.csv'
            with open(path, 'w') as f:
                f.writelines(lines)
            QMessageBox.information(dlg, "Saved", f"Data saved to {path}")
        save_btn.clicked.connect(_save)
        btn_row.addWidget(save_btn)

        copy_btn = QPushButton("Copy to Clipboard")
        def _copy_all():
            text = '\t'.join(headers) + '\n' if headers else ''
            text += '\n'.join('\t'.join(row) for row in data_rows)
            QApplication.clipboard().setText(text)
            copy_btn.setText("Copied!")
            from PySide6.QtCore import QTimer
            QTimer.singleShot(1500, lambda: copy_btn.setText("Copy to Clipboard"))
        copy_btn.clicked.connect(_copy_all)
        btn_row.addWidget(copy_btn)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dlg.close)
        btn_row.addWidget(close_btn)

        dl.addLayout(btn_row)
        dlg.show()

    # ---- Style Dialog ----

    def _show_style_dialog(self):
        W = 150
        dlg = QDialog(self.frame)
        dlg.setWindowTitle("Plot Options")
        dlg.setMinimumWidth(300)
        dl = QVBoxLayout(dlg)

        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Color"))
        from PySide6.QtWidgets import QComboBox as _QCB
        color_combo = _QCB()
        color_combo.setFixedWidth(W)
        color_combo.addItems([
            "Gradient(viridis)", "Gradient(hot)", "Gradient(jet)", "Gradient(coolwarm)",
            "Fixed(tab10)", "Fixed(tab20)", "Fixed(Set1)", "Fixed(Set2)", "Fixed(Set3)",
        ])
        color_combo.setCurrentText(getattr(self, '_color_mode', "Gradient(viridis)"))
        color_row.addWidget(color_combo)
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
        self._style_dialog = dlg
        dlg.show()

    def _get_plot_colors(self, n):
        import matplotlib.pyplot as plt
        mode = getattr(self, '_color_mode', 'Gradient(viridis)')
        start = mode.find('('); end = mode.find(')')
        cmap_name = mode[start+1:end] if start != -1 and end != -1 else 'viridis'
        cmap = plt.get_cmap(cmap_name)
        if mode.startswith('Fixed'):
            return [cmap.colors[i % len(cmap.colors)] for i in range(n)]
        else:  # Gradient
            if n == 1: return [cmap(0.5)]
            return [cmap(i / (n - 1)) for i in range(n)]

    # ---- Plotting ----

    def _parse_entry(self, entry):
        """Parse '039551_055 (Bi)' -> (shot, psi_idx)"""
        main = entry.split('(')[0].strip()
        parts = main.split('_')
        shot = int(parts[0])
        psi_idx = int(parts[1])
        return shot, psi_idx

    def _plot(self):
        entries = [self.selected_listbox.item(i).text()
                   for i in range(self.selected_listbox.count())]
        if not entries:
            return
        if not self._current_data:
            return

        bi = self._current_data['bi']
        shot = self._current_shot

        self.figure.clear()
        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'

        axes = []
        for col, param in enumerate(self.params):
            ax = self.figure.add_subplot(1, 2, col + 1)
            ax.set_xlabel('Time [s]', fontsize=self.label_fontsize)
            ax.set_ylabel(f'{LABELS.get(param, param)} [{UNITS.get(param, "")}]',
                          fontsize=self.label_fontsize)
            ax.tick_params(labelsize=self.tick_fontsize)
            ax.grid(ls='--', lw=0.3, color='#444444')
            if param == 'vT':
                ax.axhline(y=0, color=zc, ls='--', gid='zero_ref')
            axes.append(ax)

        colors = self._get_plot_colors(len(entries))
        for idx, entry in enumerate(entries):
            s, psi_idx = self._parse_entry(entry)
            color = colors[idx]

            for col, param in enumerate(self.params):
                ax = axes[col]
                if param not in bi:
                    continue

                d = bi[param]
                psin = d['psin']
                psi_actual = psin[psi_idx] if psi_idx < len(psin) else psin[-1]
                time = d['time']
                mean_t = d['mean'][psi_idx, :]
                unc_t = d['unc'][psi_idx, :]
                valid = ~np.isnan(mean_t)

                if np.any(valid):
                    ax.plot(time[valid], mean_t[valid], '-', color=color, lw=1.5,
                            label=f'\u03c8\u2099={psi_actual:.3f}')
                    ax.fill_between(time[valid],
                                    mean_t[valid] - unc_t[valid],
                                    mean_t[valid] + unc_t[valid],
                                    alpha=0.15, color=color)

        for col, param in enumerate(self.params):
            ax = axes[col]
            if param in _ZERO_BOTTOM:
                yl = ax.get_ylim()
                ax.set_ylim(0, yl[1])
            if ax.get_legend_handles_labels()[1]:
                ax.legend(fontsize=self.legend_fontsize, loc='best', frameon=False)

        self.figure.subplots_adjust(
            left=0.10, right=0.97, top=0.95, bottom=0.10, wspace=0.20)
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()
