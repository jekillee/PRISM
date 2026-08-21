"""
BiProfile Profile tab - Available/Selected workflow with Preview dialog.
Plots multiple fitted profiles (one per selected time) with raw CES/Thomson overlay.
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

from ui.ui_constants import save_file_async
from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, apply_dark_figure_style, apply_listbox_arrow_icons,
    apply_shot_arrow_icons,
)
from ui.theme import ThemeManager

UNITS = {'Ti': 'keV', 'vT': 'km/s', 'Te': 'keV', 'ne': r'$10^{19}$/m$^3$'}
UNITS_QT = {'Ti': 'keV', 'vT': 'km/s', 'Te': 'keV', 'ne': '1e19/m3'}
LABELS = {'Ti': r'T$_i$', 'vT': r'v$_T$', 'Te': r'T$_e$', 'ne': r'n$_e$'}
_ZERO_BOTTOM = {'Ti', 'Te', 'ne'}
_DEFAULT_COLOR = '#1f77b4'

COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']


class BiProfileTab:
    """BiProfile Profile tab: Available/Selected + Preview + multi-profile Plot."""

    def __init__(self, main_window, params):
        self.main = main_window
        self.params = params
        self._current_shot = None
        self._current_data = None

        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10
        self._apply_scale = True

        param_key = ''.join(p.lower() for p in params)
        self._settings_key = f"bi_{param_key}_profile"

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
        up.setFixedSize(24, 15); up.setStyleSheet(sty)
        up.clicked.connect(lambda: self._adjust_shot(1)); bl.addWidget(up)
        dn = QPushButton()
        dn.setFixedSize(24, 15); dn.setStyleSheet(sty)
        dn.clicked.connect(lambda: self._adjust_shot(-1)); bl.addWidget(dn)
        apply_shot_arrow_icons(up, dn)
        grid.addWidget(btn_updown, 0, 2)
        self.fetch_button = QPushButton("Fetch")
        self.fetch_button.setFixedWidth(70)
        self.fetch_button.clicked.connect(self._fetch)
        grid.addWidget(self.fetch_button, 0, 3)

        parent.layout().addWidget(group)

    # ---- 2. Select Data ----

    def _create_select_data(self, parent):
        group = QGroupBox("2. Select Data")
        layout = QVBoxLayout(group)

        # Preview button
        self.browse_button = QPushButton("Fetch a shot to browse")
        self.browse_button.setEnabled(False)
        self.browse_button.clicked.connect(self._open_browse)
        layout.addWidget(self.browse_button)

        # Available / Selected lists
        lists_row = QHBoxLayout()

        avail_col = QVBoxLayout()
        avail_col.addWidget(QLabel("Time [ms]"))
        self.available_listbox = QListWidget()
        self.available_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.available_listbox.setFixedHeight(150)
        avail_col.addWidget(self.available_listbox)
        lists_row.addLayout(avail_col, stretch=1)

        btn_col = QVBoxLayout()
        btn_col.addStretch()
        add_btn = QPushButton()
        add_btn.setFixedWidth(30)
        add_btn.clicked.connect(self._add_items)
        btn_col.addWidget(add_btn)
        rm_btn = QPushButton()
        rm_btn.setFixedWidth(30)
        rm_btn.clicked.connect(self._remove_items)
        btn_col.addWidget(rm_btn)
        btn_col.addStretch()
        apply_listbox_arrow_icons(add_btn, rm_btn)
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
        from ui.widgets.toggle_switch import ToggleSwitch

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

        sep = QFrame(); sep.setFrameShape(QFrame.HLine); sep.setFrameShadow(QFrame.Sunken)
        layout.addWidget(sep)

        # Show Nodes toggle
        nodes_row = QHBoxLayout()
        self.show_nodes_toggle = ToggleSwitch()
        nodes_row.addWidget(self.show_nodes_toggle)
        nodes_row.addWidget(QLabel("Show Nodes"))
        nodes_row.addStretch()
        layout.addLayout(nodes_row)

        # Deapply TS Scale toggle (ne,Te only) — scale is applied by default
        self._apply_scale = True
        if 'ne' in self.params or 'Te' in self.params:
            scale_row = QHBoxLayout()
            self.scale_toggle = ToggleSwitch()
            self.scale_toggle.toggled.connect(lambda c: setattr(self, '_apply_scale', not c))
            scale_row.addWidget(self.scale_toggle)
            scale_row.addWidget(QLabel("Unapply TS ne Scale"))
            scale_row.addStretch()
            layout.addLayout(scale_row)

        sep2 = QFrame(); sep2.setFrameShape(QFrame.HLine); sep2.setFrameShadow(QFrame.Sunken)
        layout.addWidget(sep2)

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
        if "apply_scale" in s and hasattr(self, 'scale_toggle'):
            self._apply_scale = s["apply_scale"]
            self.scale_toggle.setChecked(not self._apply_scale, animate=False)

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
        s["apply_scale"] = self._apply_scale
        set_tab_settings(self._settings_key, s)

    # ---- Actions ----

    def _adjust_shot(self, delta):
        try:
            self.shot_entry.setText(str(max(1, int(self.shot_entry.text()) + delta)))
        except ValueError:
            pass

    def _update_status(self, message, color='blue'):
        from ui.ui_constants import show_status
        show_status(self.frame, "BiProfile Profile", message, color)

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

        # Populate available list
        self.available_listbox.clear()
        for t in bi[ref]['time']:
            self.available_listbox.addItem(f"{shot:06d}_{t*1e3:06.0f} (Bi)")

        loaded = [p for p in self.params if p in bi]
        self._update_status(f"#{shot}: {', '.join(loaded)} loaded", color='green')
        self.info_label.setText(
            f"EFIT: {bi[ref].get('efit_used','?')}\n"
            f"Fit: {bi[ref].get('fit_func','?')}\n"
            f"{len(bi[ref]['time'])} time pts, {len(bi[ref]['psin'])} \u03c8\u2099 pts")

        self.browse_button.setEnabled(True)
        self.browse_button.setText(f"Browse #{shot}")

    def _open_browse(self):
        if not self._current_data:
            return
        from ui.widgets.preview_dialog import BiProfilePreviewDialog
        dlg = BiProfilePreviewDialog(
            self.frame, self._current_data['bi'], self._current_shot,
            self.params, mode='profile',
            selected_listbox=self.selected_listbox,
            sdata=self._current_data)
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
        lines = []
        lines.append(f"# BiProfile Data\n")
        header_parts = ["Shot", "Time[s]", "psi_N"]
        for p in self.params:
            units = UNITS.get(p, '')
            header_parts.append(f"{p}[{units}]")
            header_parts.append(f"{p}_unc")
        lines.append("#" + ','.join(f"{h:>12s}" for h in header_parts) + "\n")

        for entry in entries:
            s, t_target = self._parse_entry(entry)
            for col, param in enumerate(self.params):
                if param not in bi:
                    continue
                d = bi[param]
                t_idx = np.argmin(np.abs(d['time'] - t_target))
                t_actual = d['time'][t_idx]
                psin = d['psin']
                mean_p = d['mean'][:, t_idx]
                unc_p = d['unc'][:, t_idx]

                if col == 0:
                    # First param: store psin and values
                    n_psi = len(psin)
                    shot_col = [str(s)] * n_psi
                    time_col = [f"{t_actual:.4f}"] * n_psi
                    psin_col = psin
                    data_cols = {param: (mean_p, unc_p)}
                else:
                    data_cols[param] = (mean_p, unc_p)

            for i in range(n_psi):
                row = f" {shot_col[i]:>12s},{time_col[i]:>12s},{psin_col[i]:>12.6f}"
                for p in self.params:
                    if p in data_cols:
                        m, u = data_cols[p]
                        row += f",{m[i]:>12.6e},{u[i]:>12.6e}"
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
            def _write(path):
                if not path:
                    return
                if not os.path.splitext(path)[1]:
                    path += '.csv'
                with open(path, 'w') as f:
                    f.writelines(lines)
                QMessageBox.information(dlg, "Saved", f"Data saved to {path}")
            # Non-modal save dialog keeps the main window movable while it's open
            save_file_async(dlg, "Save Data", os.path.expanduser("~"),
                            "CSV files (*.csv);;All files (*)", _write)
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
        """Get n colors from current color mode"""
        import matplotlib.pyplot as plt
        mode = getattr(self, '_color_mode', 'Gradient(viridis)')
        start = mode.find('('); end = mode.find(')')
        cmap_name = mode[start+1:end] if start != -1 and end != -1 else 'viridis'
        cmap = plt.get_cmap(cmap_name)
        if mode.startswith('Fixed'):
            return [cmap.colors[i % len(cmap.colors)] for i in range(n)]
        else:  # Gradient
            if n == 1:
                return [cmap(0.5)]
            return [cmap(i / (n - 1)) for i in range(n)]

    # ---- Plotting ----

    def _parse_entry(self, entry):
        """Parse '039551_005020 (Bi)' → (shot, time_s)"""
        main = entry.split('(')[0].strip()
        parts = main.split('_')
        shot = int(parts[0])
        time_s = float(parts[1]) / 1e3
        return shot, time_s

    def _plot(self):
        entries = [self.selected_listbox.item(i).text()
                   for i in range(self.selected_listbox.count())]
        if not entries:
            return
        if not self._current_data:
            return

        bi = self._current_data['bi']
        sdata = self._current_data
        shot = self._current_shot

        self.figure.clear()
        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'

        axes = []
        for col, param in enumerate(self.params):
            ax = self.figure.add_subplot(1, 2, col + 1)
            ax.set_xlabel(r'$\psi_N$', fontsize=self.label_fontsize)
            ax.set_ylabel(f'{LABELS.get(param, param)} [{UNITS.get(param, "")}]',
                          fontsize=self.label_fontsize)
            ax.tick_params(labelsize=self.tick_fontsize)
            ax.set_xlim(0, 1.1)
            ax.grid(ls='--', lw=0.3, color='#444444')
            ax.axvline(x=1.0, color=zc, ls='--', lw=0.8, alpha=0.5)
            axes.append(ax)

        colors = self._get_plot_colors(len(entries))
        for idx, entry in enumerate(entries):
            s, t_target = self._parse_entry(entry)
            color = colors[idx]
            ms = int(round(t_target * 1000))

            for col, param in enumerate(self.params):
                ax = axes[col]
                if param not in bi:
                    continue

                d = bi[param]
                t_idx = np.argmin(np.abs(d['time'] - t_target))
                t_actual = d['time'][t_idx]
                psin = d['psin']
                mean_p = d['mean'][:, t_idx]
                unc_p = d['unc'][:, t_idx]
                valid = ~np.isnan(mean_p)

                if np.any(valid):
                    ax.plot(psin[valid], mean_p[valid], '-', color=color, lw=2,
                            label=f'#{s} {ms:06d}ms (Bi)')
                    ax.fill_between(psin[valid],
                                    mean_p[valid] - unc_p[valid],
                                    mean_p[valid] + unc_p[valid],
                                    alpha=0.15, color=color)

                # Raw overlay
                self._overlay_raw(ax, param, t_actual, color)

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

    # ---- Raw overlay ----

    def _overlay_raw(self, ax, param, t_actual, color):
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
        bi = self._current_data['bi']
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

    @staticmethod
    def _plot_raw(ax, psin, y, yerr, use_flags, color, marker):
        yerr = np.where(yerr < 0, 0, yerr)
        valid = (psin >= 0) & (psin <= 1.2) & np.isfinite(y)
        used = valid & (use_flags == 1)
        excluded = valid & (use_flags == 0)
        unknown = valid & np.isnan(use_flags)

        if np.any(used):
            ax.errorbar(psin[used], y[used], yerr=yerr[used],
                        fmt=marker, markersize=5, color=color,
                        capsize=5, zorder=10, markeredgewidth=1)
        if np.any(excluded):
            ax.errorbar(psin[excluded], y[excluded], yerr=yerr[excluded],
                        fmt=marker, markersize=5, color=(0.6, 0.6, 0.6, 0.35),
                        capsize=5, zorder=5, markeredgewidth=1)
        if np.any(unknown):
            ax.errorbar(psin[unknown], y[unknown], yerr=yerr[unknown],
                        fmt=marker, markersize=5, color=color,
                        capsize=5, zorder=10, markeredgewidth=1)
