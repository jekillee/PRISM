"""
EFIT Scalar tab - load EFIT AEQDSK scalar time traces from MDS+ or a-files
and plot selected variables. Supports multiple shots/trees for comparison.
"""

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QSplitter, QScrollArea, QLayout,
    QGroupBox, QLabel, QLineEdit, QPushButton, QFrame, QComboBox,
    QListWidget, QAbstractItemView, QFileDialog,
    QApplication, QMessageBox, QStyle, QSpinBox,
    QDialog, QDialogButtonBox,
    QTableWidget, QTableWidgetItem, QHeaderView,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QShortcut, QKeySequence

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style, get_icon
from ui.theme import ThemeManager


class EfitScalarTab:
    """EFIT Scalar tab: load AEQDSK data -> select variables -> plot."""

    def __init__(self, main_window):
        self.main = main_window
        self._data_cache = {}  # {cache_key: efit_data}

        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10
        self._settings_key = "efit_scalar"

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

    # ---- 1. Load Data ----

    def _create_load_section(self, parent):
        group = QGroupBox("1. Load Data")
        layout = QVBoxLayout(group)

        grid = QGridLayout()
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
        up.setIcon(get_icon(QStyle.SP_ArrowUp))
        up.setFixedSize(24, 15); up.setStyleSheet(sty)
        up.clicked.connect(lambda: self._adjust_shot(1)); bl.addWidget(up)
        dn = QPushButton()
        dn.setIcon(get_icon(QStyle.SP_ArrowDown))
        dn.setFixedSize(24, 15); dn.setStyleSheet(sty)
        dn.clicked.connect(lambda: self._adjust_shot(-1)); bl.addWidget(dn)
        grid.addWidget(btn_updown, 0, 2)

        self.tree_combo = QComboBox()
        self.tree_combo.addItems([
            "efitrt1", "efitrt2", "efit01", "efit02", "efit04",
        ])
        self.tree_combo.setCurrentText("efitrt1")
        grid.addWidget(self.tree_combo, 0, 3)

        fetch_btn = QPushButton("Fetch")
        fetch_btn.setFixedWidth(70)
        fetch_btn.clicked.connect(self._fetch)
        grid.addWidget(fetch_btn, 0, 4)
        layout.addLayout(grid)

        # Open a-File
        open_btn = QPushButton("Open a-File...")
        open_btn.clicked.connect(self._open_file)
        layout.addWidget(open_btn)

        self.load_status = QLabel("")
        self.load_status.setWordWrap(True)
        layout.addWidget(self.load_status)

        parent.layout().addWidget(group)

    # ---- 2. Select Data ----

    def _create_select_section(self, parent):
        group = QGroupBox("2. Select Data")
        group.setEnabled(False)
        self._select_group = group
        layout = QVBoxLayout(group)

        # Browse button
        self.preview_button = QPushButton("Load data to browse")
        self.preview_button.setEnabled(False)
        self.preview_button.clicked.connect(self._open_preview)
        layout.addWidget(self.preview_button)

        # Selected runs list
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

    # ---- 3. Select Variable ----

    def _create_var_section(self, parent):
        group = QGroupBox("3. Select Variable")
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

    def _make_cache_key(self, shot, tree):
        return f"{shot}_{tree}"

    def _parse_cache_key(self, key):
        """Split cache key into (shot, tree).
        MDS+: '39551_efit01' -> ('39551', 'efit01')
        a-file: 'a023879.011700_kin_0' -> ('23879', 'kin_0 @11700ms')
        """
        import re
        # Try a/g-file pattern first
        m = re.match(r'^[ag](\d+)\.(\d+)_(.+)$', key)
        if m:
            return m.group(1), f"{m.group(3)} @{m.group(2)}ms"
        m = re.match(r'^[ag](\d+)\.(\d+)$', key)
        if m:
            return m.group(1), f"@{m.group(2)}ms"
        # MDS+ pattern
        parts = key.split('_', 1)
        if len(parts) == 2:
            return parts[0], parts[1]
        return key, ''

    def _adjust_shot(self, delta):
        try:
            self.shot_entry.setText(str(max(1, int(self.shot_entry.text()) + delta)))
        except ValueError:
            pass

    def _fetch(self):
        shot_text = self.shot_entry.text().strip()
        if not shot_text:
            QMessageBox.warning(self.frame, "Warning", "Enter a shot number.")
            return
        try:
            shot_number = int(shot_text)
        except ValueError:
            QMessageBox.warning(self.frame, "Warning", "Shot number must be an integer.")
            return

        tree = self.tree_combo.currentText()
        cache_key = self._make_cache_key(shot_number, tree)

        if cache_key in self._data_cache:
            return

        self.load_status.setStyleSheet("color: blue; font-weight: bold; font-size: 9pt;")
        self.load_status.setText(f"Loading #{shot_number} ({tree})...")
        QApplication.setOverrideCursor(Qt.WaitCursor)
        QApplication.processEvents()
        try:
            from data_loaders.efit_viewer_loader import load_efit_mds
            efit_data = load_efit_mds(shot_number, tree, load_profiles=False, load_2d=False)
            self._data_cache[cache_key] = efit_data
        except Exception as e:
            QApplication.restoreOverrideCursor()
            self.load_status.setStyleSheet("color: red; font-weight: bold; font-size: 9pt;")
            self.load_status.setText(f"Failed: {e}")
            return
        finally:
            QApplication.restoreOverrideCursor()

        self._update_var_combo()
        self._update_selected_runs()
        self._update_load_status()
        self._select_group.setEnabled(bool(self._data_cache))
        self._var_group.setEnabled(bool(self._data_cache))

        if self._data_cache:
            self._preview_label = cache_key
            shot, tree_name = self._parse_cache_key(cache_key)
            self.preview_button.setEnabled(True)
            self.preview_button.setText(f"Browse #{shot} ({tree_name})")

    def _open_file(self):
        start_dir = getattr(self, '_last_file_dir', None)
        if start_dir is None:
            import os
            start_dir = os.path.expanduser("~")

        dlg = QFileDialog(self.frame, "Open EFIT a-File", start_dir,
                          "EFIT a-files (a*);;All files (*)")
        dlg.setFileMode(QFileDialog.ExistingFiles)
        dlg.setOption(QFileDialog.DontUseNativeDialog, True)
        dlg.setWindowModality(Qt.NonModal)
        dlg.filesSelected.connect(self._on_files_selected)
        dlg.show()

    def _on_files_selected(self, paths):
        import os
        if not paths:
            return
        self._last_file_dir = os.path.dirname(paths[0])

        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            from data_loaders.efit_viewer_loader import load_afile
            for path in paths:
                basename = os.path.basename(path)
                if basename in self._data_cache:
                    continue
                efit_data = load_afile(path)
                self._data_cache[basename] = efit_data
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load file:\n{e}")
        finally:
            QApplication.restoreOverrideCursor()

        self._update_var_combo()
        self._update_selected_runs()
        self._update_load_status()
        self._select_group.setEnabled(bool(self._data_cache))
        self._var_group.setEnabled(bool(self._data_cache))

        if self._data_cache:
            self._preview_label = list(self._data_cache.keys())[-1]
            self.preview_button.setEnabled(True)
            self.preview_button.setText(f"Browse ({len(self._data_cache)} loaded)")

    def _open_preview(self):
        """Browse EFIT scalars — all runs, variable selector."""
        if not self._data_cache:
            return

        # Build pseudo-cache for TranspTimeTracePreviewDialog — all runs, sorted
        pseudo_cache = {}
        run_labels = []
        for key, efit in sorted(self._data_cache.items()):
            display_label = key.replace('_', '', 1)
            pseudo_cache[display_label] = {'timetraces': {}}
            run_labels.append(display_label)
            for vn, vd in efit.get('scalars', {}).items():
                pseudo_cache[display_label]['timetraces'][vn] = {
                    'time': vd['time'],
                    'data': vd['data'],
                    'long_name': vd.get('label', vn),
                'units': '',
            }

        from ui.widgets.preview_dialog import TranspTimeTracePreviewDialog
        dlg = TranspTimeTracePreviewDialog(self.frame, pseudo_cache, run_labels)
        dlg.show()

    def _update_load_status(self, color='green'):
        if not self._data_cache:
            self.load_status.setStyleSheet(
                "color: #888; font-weight: bold; font-size: 9pt;")
            self.load_status.setText("No data loaded")
            return
        n_runs = len(self._data_cache)
        n_v = max(len(d.get('scalars', {})) for d in self._data_cache.values())
        self.load_status.setStyleSheet(
            f"color: {color}; font-weight: bold; font-size: 9pt;")
        self.load_status.setText(f"{n_runs} run(s) loaded  —  {n_v} variables")

    def _update_var_combo(self):
        """Refresh variable combo from intersection of all loaded data."""
        filt = self.var_filter.text().strip().upper()
        current_var = self.var_combo.currentData()

        # Intersection of variable names (skip empty scalar dicts)
        var_sets = [set(efit.get('scalars', {}).keys())
                    for efit in self._data_cache.values()
                    if efit.get('scalars')]
        common_vars = set.intersection(*var_sets) if var_sets else set()

        all_vars = {}
        first_efit = next(iter(self._data_cache.values()), {})
        for vn in common_vars:
            vd = first_efit.get('scalars', {}).get(vn, {})
            parts = [vn]
            label = vd.get('label', '')
            if label and label != vn:
                parts.append(f"- {label}")
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
        """Add newly loaded run entries to selected list."""
        existing = {self.selected_listbox.item(i).data(Qt.UserRole)
                    for i in range(self.selected_listbox.count())}
        for key in sorted(self._data_cache):
            if key not in existing:
                shot, tree = self._parse_cache_key(key)
                from PySide6.QtWidgets import QListWidgetItem
                item = QListWidgetItem(f"#{shot} ({tree})")
                item.setData(Qt.UserRole, key)
                self.selected_listbox.addItem(item)

    def _remove_items(self):
        for item in reversed(self.selected_listbox.selectedItems()):
            key = item.data(Qt.UserRole)
            self.selected_listbox.takeItem(self.selected_listbox.row(item))
            self._data_cache.pop(key, None)
        self._update_var_combo()
        self._update_load_status()
        if not self._data_cache:
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
        run_keys = [self.selected_listbox.item(i).data(Qt.UserRole)
                    for i in range(self.selected_listbox.count())]
        if not run_keys:
            return

        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)
        ax.set_xlabel('Time [s]', fontsize=self.label_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(ls='--', lw=0.3, color='#444444')

        # Filter runs that have this variable
        plot_items = []
        for key in run_keys:
            efit = self._data_cache.get(key)
            if efit and var_name in efit.get('scalars', {}):
                plot_items.append((key, efit['scalars'][var_name]))

        if not plot_items:
            return

        colors = self._get_plot_colors(len(plot_items))
        for idx, (key, sd) in enumerate(plot_items):
            try:
                shot, tree = self._parse_cache_key(key)
                if len(sd['data']) == 1:
                    # Single point (a-file) — marker only
                    ax.plot(sd['time'], sd['data'], 'o', color=colors[idx],
                            markersize=6, label=f'#{shot} ({tree})')
                else:
                    ax.plot(sd['time'], sd['data'], '-', color=colors[idx],
                            lw=1.5, marker='.', markersize=3,
                            label=f'#{shot} ({tree})')
            except Exception as e:
                print(f"[EFIT] Error plotting {var_name} ({key}): {e}")

        # Y-axis label from combo display text
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
    # Preview & Save
    # ================================================================

    def _build_csv_lines(self):
        """Build CSV lines in PRISM export format."""
        var_name = self.var_combo.currentData()
        run_keys = [self.selected_listbox.item(i).data(Qt.UserRole)
                    for i in range(self.selected_listbox.count())]
        if not var_name or not run_keys:
            return None

        units = ''
        for efit in self._data_cache.values():
            if var_name in efit.get('scalars', {}):
                units = efit['scalars'][var_name].get('units', '')
                break

        lines = []
        lines.append(f"# EFIT Scalar: {var_name}\n")
        lines.append(f"#%10s,%10s,%10s,%14s\n" % (
            "Shot", "Tree", "Time[s]", f"{var_name}[{units}]"))

        for key in run_keys:
            efit = self._data_cache.get(key)
            if efit is None or var_name not in efit.get('scalars', {}):
                continue
            sd = efit['scalars'][var_name]
            shot, tree = self._parse_cache_key(key)

            for t, v in zip(sd['time'], sd['data']):
                lines.append(" %10s,%10s,%10.4f,%14.6e\n" % (shot, tree, t, v))

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
