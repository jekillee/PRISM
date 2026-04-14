"""
TRANSP CDF Time Trace tab - load CDF files and plot time trace variables.
Supports multiple CDFs for comparison.
"""

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QScrollArea, QLayout,
    QGroupBox, QLabel, QLineEdit, QPushButton, QFrame,
    QListWidget, QAbstractItemView, QComboBox, QFileDialog,
    QApplication, QMessageBox, QStyle, QSpinBox,
    QDialog, QDialogButtonBox,
    QTableWidget, QTableWidgetItem, QHeaderView,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QShortcut, QKeySequence

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style, get_icon
from ui.theme import ThemeManager


class TranspTimeTraceTab:
    """TRANSP CDF Time Trace tab: load CDF → select variables → plot."""

    def __init__(self, main_window):
        self.main = main_window
        self._cdf_cache = {}  # {label: cdf_data}

        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10
        self._settings_key = "transp_timetrace"

        self.frame = QWidget()
        self.canvas = None
        self._create_widgets()
        self._restore_settings()

    # ================================================================
    # Widget creation
    # ================================================================

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

        self._create_load_section(control_frame)
        self._create_select_section(control_frame)
        self._create_var_section(control_frame)
        self._create_plot_controls(control_frame)
        self._create_save_section(control_frame)
        control_layout.addStretch()

    # ---- 1. Load CDF ----

    def _create_load_section(self, parent):
        group = QGroupBox("1. Load CDF")
        layout = QVBoxLayout(group)

        open_btn = QPushButton("Open CDF...")
        open_btn.clicked.connect(self._open_cdf)
        layout.addWidget(open_btn)

        self.cdf_status = QLabel("")
        self.cdf_status.setWordWrap(True)
        layout.addWidget(self.cdf_status)

        parent.layout().addWidget(group)

    # ---- 2. Select Data ----

    def _create_select_section(self, parent):
        group = QGroupBox("2. Select Data")
        group.setEnabled(False)
        self._select_group = group
        layout = QVBoxLayout(group)

        # Preview button
        self.preview_button = QPushButton("Load a CDF to browse")
        self.preview_button.setEnabled(False)
        self.preview_button.clicked.connect(self._open_preview)
        layout.addWidget(self.preview_button)

        # Selected run IDs (populated on CDF load)
        layout.addWidget(QLabel("Selected Runs"))
        self.selected_listbox = QListWidget()
        self.selected_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.selected_listbox.setFixedHeight(100)
        layout.addWidget(self.selected_listbox)

        rm_btn = QPushButton("Remove Selected Run")
        rm_btn.clicked.connect(self._remove_items)
        layout.addWidget(rm_btn)

        QShortcut(QKeySequence(Qt.Key_Delete), self.selected_listbox
                  ).activated.connect(self._remove_items)

        parent.layout().addWidget(group)

    def _create_var_section(self, parent):
        group = QGroupBox("3. Select TRANSP Variable")
        group.setEnabled(False)
        self._var_group = group
        layout = QVBoxLayout(group)

        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Filter"))
        self.var_filter = QLineEdit()
        self.var_filter.setPlaceholderText("Type to filter...")
        self.var_filter.textChanged.connect(self._update_var_combo)
        filter_row.addWidget(self.var_filter, 1)
        layout.addLayout(filter_row)

        self.var_combo = QComboBox()
        self.var_combo.setMaxVisibleItems(20)
        self.var_combo.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.var_combo.setMinimumContentsLength(10)
        self.var_combo.setStyleSheet("QComboBox { combobox-popup: 0; }")
        self.var_combo.view().setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        layout.addWidget(self.var_combo)

        parent.layout().addWidget(group)

    # ---- 4. Plot ----

    def _create_plot_controls(self, parent):
        group = QGroupBox("4. Plot")
        layout = QVBoxLayout(group)

        btn_row = QHBoxLayout()
        plot_btn = QPushButton("Plot")
        plot_btn.clicked.connect(self._plot)
        btn_row.addWidget(plot_btn, 3)
        opt_btn = QPushButton("Option")
        opt_btn.clicked.connect(self._show_style_dialog)
        btn_row.addWidget(opt_btn, 1)
        layout.addLayout(btn_row)

        parent.layout().addWidget(group)

    # ---- 5. Save Data ----

    def _create_save_section(self, parent):
        group = QGroupBox("5. Save Data")
        layout = QVBoxLayout(group)
        btn = QPushButton("Preview && Save")
        btn.clicked.connect(self._preview_and_save)
        layout.addWidget(btn)
        parent.layout().addWidget(group)

    # ================================================================
    # Actions
    # ================================================================

    @staticmethod
    def _default_cdf_dir():
        import os, socket
        host = socket.gethostname()
        user = os.environ.get('USER', os.environ.get('USERNAME', ''))
        if host.startswith('nkstar'):
            return f"/home/users/{user}"
        elif host.startswith('ukstar'):
            transp_dir = f"/UKSTAR_HOME/data/transp/{user}"
            if os.path.isdir(transp_dir):
                return transp_dir
            return f"/UKSTAR_HOME/{user}"
        return os.path.expanduser("~")

    @staticmethod
    def _parse_run_label(label):
        """Split run label into (shot_str, run_suffix). e.g. '39551W09' → ('39551','W09')"""
        import re
        m = re.match(r'^(\d+)(.*)', label)
        if m:
            return m.group(1), m.group(2)
        return label, ''

    def _open_cdf(self):
        start_dir = getattr(self, '_last_cdf_dir', self._default_cdf_dir())
        dlg = QFileDialog(self.frame, "Open TRANSP CDF", start_dir,
                          "TRANSP CDF (*.CDF);;All files (*)")
        dlg.setFileMode(QFileDialog.ExistingFile)
        dlg.setWindowModality(Qt.NonModal)
        dlg.fileSelected.connect(self._on_cdf_selected)
        dlg.show()

    def _on_cdf_selected(self, path):
        import os
        self._last_cdf_dir = os.path.dirname(path)

        # Skip if same file already loaded
        basename = os.path.splitext(os.path.basename(path))[0]
        if basename in self._cdf_cache:
            return

        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            from data_loaders.transp_cdf_loader import load_transp_cdf
            cdf = load_transp_cdf(path)
            self._cdf_cache[cdf['label']] = cdf
        finally:
            QApplication.restoreOverrideCursor()

        self._update_var_combo()
        self._update_selected_runs()
        self._update_cdf_status()
        self._select_group.setEnabled(bool(self._cdf_cache))
        self._var_group.setEnabled(bool(self._cdf_cache))

        if self._cdf_cache:
            label = list(self._cdf_cache.keys())[-1]
            self._preview_label = label
            shot, run = self._parse_run_label(label)
            self.preview_button.setEnabled(True)
            self.preview_button.setText(f"Browse #{shot} ({run})")

    def _open_preview(self):
        label = getattr(self, '_preview_label', None)
        if not label or label not in self._cdf_cache:
            return

        from ui.widgets.preview_dialog import TranspTimeTracePreviewDialog
        dlg = TranspTimeTracePreviewDialog(
            self.frame, {label: self._cdf_cache[label]}, [label])
        dlg.show()

    def _update_cdf_status(self, color='green'):
        if not self._cdf_cache:
            self.cdf_status.setStyleSheet("color: #888; font-weight: bold; font-size: 9pt;")
            self.cdf_status.setText("No CDF loaded")
            return
        parts = []
        for label in sorted(self._cdf_cache):
            shot, run = self._parse_run_label(label)
            parts.append(f"#{shot} ({run})")
        n_t = sum(len(c['timetraces']) for c in self._cdf_cache.values())
        self.cdf_status.setStyleSheet(f"color: {color}; font-weight: bold; font-size: 9pt;")
        self.cdf_status.setText(f"{', '.join(parts)}  —  {n_t} time traces")

    def _update_var_combo(self):
        """Refresh variable combo from all loaded CDFs, applying filter."""
        filt = self.var_filter.text().strip().upper()
        current_var = self.var_combo.currentData()

        all_vars = {}
        for cdf in self._cdf_cache.values():
            for vn, vd in cdf['timetraces'].items():
                if vn not in all_vars:
                    parts = [vn]
                    if vd['long_name'] and vd['long_name'] != vn:
                        parts.append(f"- {vd['long_name']}")
                    if vd['units']:
                        parts.append(f"[{vd['units']}]")
                    all_vars[vn] = ' '.join(parts)

        filtered = sorted(all_vars.items())
        if filt:
            filtered = [(k, v) for k, v in filtered
                        if filt in k.upper() or filt in v.upper()]

        self.var_combo.blockSignals(True)
        self.var_combo.clear()
        for vn, display in filtered:
            self.var_combo.addItem(display, vn)
        restore_var = current_var or getattr(self, '_pending_variable', None)
        if restore_var:
            for i in range(self.var_combo.count()):
                if self.var_combo.itemData(i) == restore_var:
                    self.var_combo.setCurrentIndex(i)
                    break
            self._pending_variable = None
        self.var_combo.blockSignals(False)

    def _update_selected_runs(self):
        """Add newly loaded CDF run IDs to selected list."""
        existing = {self.selected_listbox.item(i).data(Qt.UserRole)
                    for i in range(self.selected_listbox.count())}
        for label in sorted(self._cdf_cache):
            if label not in existing:
                shot, run = self._parse_run_label(label)
                from PySide6.QtWidgets import QListWidgetItem
                item = QListWidgetItem(f"#{shot} ({run})")
                item.setData(Qt.UserRole, label)
                self.selected_listbox.addItem(item)

    def _remove_items(self):
        for item in reversed(self.selected_listbox.selectedItems()):
            label = item.data(Qt.UserRole)
            self.selected_listbox.takeItem(self.selected_listbox.row(item))
            self._cdf_cache.pop(label, None)
        self._update_var_combo()
        self._update_cdf_status()
        if not self._cdf_cache:
            self._select_group.setEnabled(False)
            self._var_group.setEnabled(False)

    # ================================================================
    # Plotting
    # ================================================================

    def _get_plot_colors(self, n):
        import matplotlib.pyplot as plt
        mode = getattr(self, '_color_mode', 'Gradient(viridis)')
        start = mode.find('('); end = mode.find(')')
        cmap_name = mode[start+1:end] if start != -1 and end != -1 else 'viridis'
        cmap = plt.get_cmap(cmap_name)
        if mode.startswith('Fixed'):
            return [cmap.colors[i % len(cmap.colors)] for i in range(n)]
        else:
            if n == 1:
                return [cmap(0.5)]
            return [cmap(i / (n - 1)) for i in range(n)]

    def _plot(self):
        var_name = self.var_combo.currentData()
        if not var_name:
            return
        run_labels = [self.selected_listbox.item(i).data(Qt.UserRole)
                      for i in range(self.selected_listbox.count())]
        if not run_labels:
            return

        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)
        ax.set_xlabel('Time [s]', fontsize=self.label_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(ls='--', lw=0.3, color='#444444')

        # Filter runs that have this variable
        plot_items = []
        for label in run_labels:
            cdf = self._cdf_cache.get(label)
            if cdf and var_name in cdf['timetraces']:
                plot_items.append((label, cdf['timetraces'][var_name]))

        if not plot_items:
            return

        colors = self._get_plot_colors(len(plot_items))
        units_seen = set()

        for idx, (label, tt) in enumerate(plot_items):
            try:
                units_seen.add(tt['units'])
                shot, run = self._parse_run_label(label)
                ax.plot(tt['time'], tt['data'], '-', color=colors[idx],
                        lw=1.5, marker='.', markersize=3, label=f'#{shot} ({run})')
            except Exception as e:
                print(f"[TRANSP] Error plotting {var_name} ({label}): {e}")

        # Y-axis: variable name + units
        display_text = self.var_combo.currentText()
        ax.set_ylabel(display_text, fontsize=self.label_fontsize)

        if ax.get_legend_handles_labels()[1]:
            ax.legend(fontsize=self.legend_fontsize, loc='best', frameon=False)

        self.figure.subplots_adjust(
            left=0.12, right=0.97, top=0.95, bottom=0.10)
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()

    # ================================================================
    # Style dialog
    # ================================================================

    # ================================================================
    # Preview & Save
    # ================================================================

    def _build_csv_lines(self):
        """Build CSV lines in PRISM export format."""
        import re
        var_name = self.var_combo.currentData()
        run_labels = [self.selected_listbox.item(i).data(Qt.UserRole)
                      for i in range(self.selected_listbox.count())]
        if not var_name or not run_labels:
            return None

        units = ''
        for cdf in self._cdf_cache.values():
            if var_name in cdf['timetraces']:
                units = cdf['timetraces'][var_name].get('units', '')
                break

        lines = []
        lines.append(f"# TRANSP Time Trace: {var_name}\n")
        lines.append(f"#%10s,%10s,%10s,%14s\n" % (
            "Shot", "Run", "Time[s]", f"{var_name}[{units}]"))

        for label in run_labels:
            cdf = self._cdf_cache.get(label)
            if cdf is None or var_name not in cdf['timetraces']:
                continue
            tt = cdf['timetraces'][var_name]
            m = re.match(r'^(\d+)(.*)', label)
            shot = m.group(1) if m else label
            run = m.group(2) if m else ''

            for t, v in zip(tt['time'], tt['data']):
                lines.append(" %10s,%10s,%10.4f,%14.6e\n" % (shot, run, t, v))

        return lines

    def _preview_and_save(self):
        lines = self._build_csv_lines()
        if not lines:
            QMessageBox.warning(self.frame, "Warning", "No data to preview")
            return
        self._show_data_preview(lines)

    def _show_data_preview(self, lines):
        """Show data in PRISM-style spreadsheet dialog."""
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

    # ================================================================
    # Style dialog
    # ================================================================

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
            "Gradient(coolwarm)",
            "Fixed(tab10)", "Fixed(tab20)", "Fixed(Set1)", "Fixed(Set2)",
            "Fixed(Set3)",
        ])
        color_combo.setCurrentText(getattr(self, '_color_mode',
                                           "Gradient(viridis)"))
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

    # ================================================================
    # Settings
    # ================================================================

    def _restore_settings(self):
        from config.user_settings import get_tab_settings
        s = get_tab_settings(self._settings_key)
        self._color_mode = s.get("color_mode", "Gradient(viridis)")
        self.label_fontsize = s.get("label_fontsize", 12)
        self.legend_fontsize = s.get("legend_fontsize", 8)
        self.tick_fontsize = s.get("tick_fontsize", 10)
        saved_var = s.get("variable", "")
        if saved_var:
            self._pending_variable = saved_var

    def save_settings(self):
        from config.user_settings import get_tab_settings, set_tab_settings
        s = get_tab_settings(self._settings_key)
        s["color_mode"] = getattr(self, '_color_mode', 'Gradient(viridis)')
        s["label_fontsize"] = self.label_fontsize
        s["legend_fontsize"] = self.legend_fontsize
        s["tick_fontsize"] = self.tick_fontsize
        s["variable"] = self.var_combo.currentData() or ""
        set_tab_settings(self._settings_key, s)
