"""
EFIT p-File (PEQDSK) Profile tab - load p-files and plot kinetic profiles.
Plots ne, te, ni, ti, nb, pb, ptot, omeg, etc. at selected time points.
Supports multiple p-files for comparison.
"""

import os
import re
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QScrollArea, QLayout,
    QGroupBox, QGridLayout, QLabel, QLineEdit, QPushButton, QFrame,
    QListWidget, QAbstractItemView, QComboBox, QFileDialog,
    QApplication, QMessageBox, QStyle, QSpinBox,
    QDialog, QDialogButtonBox,
    QTableWidget, QTableWidgetItem, QHeaderView,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QShortcut, QKeySequence

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style, get_icon
from ui.theme import ThemeManager


def _load_pfile(filepath):
    """Parse a p-file (PEQDSK kinetic profiles).

    Returns dict with:
        'shot': str, 'tree': str (suffix), 'time': ndarray (single element, seconds),
        'psin': ndarray, 'profiles': {varname: ndarray(1, npts)},
        'units': {varname: str}, 'derivatives': {varname: ndarray(1, npts)}
    """
    basename = os.path.basename(filepath)
    m = re.match(r'p(\d+)\.(\d+)(?:_(.+))?', basename)
    if m:
        shot_str = str(int(m.group(1)))
        time_ms = int(m.group(2))
        time_s = time_ms / 1000.0
        suffix = m.group(3) if m.group(3) else 'pfile'
    else:
        shot_str = basename
        time_s = 0.0
        time_ms = 0
        suffix = 'pfile'

    with open(filepath, 'r') as f:
        lines = f.readlines()

    profiles = {}
    units = {}
    derivatives = {}
    psin = None

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        # Try to match block header:
        #   {nrows} psinorm {varname}({units}) d{varname}/dpsiN
        header_pat = re.match(
            r'^\s*(\d+)\s+psinorm\s+(\w+)\(([^)]*)\)\s+d\w+/dpsiN\s*$', line)
        if not header_pat:
            i += 1
            continue

        nrows = int(header_pat.group(1))
        varname = header_pat.group(2)
        unit_str = header_pat.group(3)
        i += 1

        psi_col = []
        val_col = []
        der_col = []

        for _ in range(nrows):
            if i >= len(lines):
                break
            tokens = lines[i].split()
            i += 1
            if len(tokens) >= 3:
                psi_col.append(float(tokens[0]))
                val_col.append(float(tokens[1]))
                der_col.append(float(tokens[2]))

        psi_arr = np.array(psi_col)
        val_arr = np.array(val_col)
        der_arr = np.array(der_col)

        # Use psinorm from the first block as the canonical grid
        if psin is None:
            psin = psi_arr

        profiles[varname] = val_arr.reshape(1, -1)
        units[varname] = unit_str
        derivatives[varname] = der_arr.reshape(1, -1)

    if psin is None:
        psin = np.array([0.0])

    return {
        'shot': shot_str,
        'tree': suffix,
        'time': np.array([time_s]),
        'psin': psin,
        'profiles': profiles,
        'units': units,
        'derivatives': derivatives,
    }


class EfitPfileTab:
    """EFIT p-File tab: load p-files -> select entries -> plot kinetic profiles."""

    def __init__(self, main_window):
        self.main = main_window
        self._pfile_cache = {}          # {cache_key: pfile_data}
        self._current_key = None        # cache key of last loaded dataset
        self._all_profile_vars = {}     # {var_name: display_text}
        self._all_units = {}            # {var_name: unit_str}

        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10
        self._settings_key = "efit_pfile"

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

    # ---- 1. Load p-File ----

    def _create_load_section(self, parent):
        group = QGroupBox("1. Load p-File")
        layout = QVBoxLayout(group)

        open_btn = QPushButton("Open p-File...")
        open_btn.clicked.connect(self._open_pfile)
        layout.addWidget(open_btn)

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        parent.layout().addWidget(group)

    # ---- 2. Select Data ----

    def _create_select_section(self, parent):
        group = QGroupBox("2. Select Data")
        group.setEnabled(False)
        self._select_group = group
        layout = QVBoxLayout(group)

        # Browse button
        self.preview_button = QPushButton("Load a p-file to browse")
        self.preview_button.setEnabled(False)
        self.preview_button.clicked.connect(self._open_preview)
        layout.addWidget(self.preview_button)

        # Selected p-files (auto-populated on load)
        layout.addWidget(QLabel("Loaded p-Files"))
        self.selected_listbox = QListWidget()
        self.selected_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.selected_listbox.setFixedHeight(120)
        layout.addWidget(self.selected_listbox)

        rm_btn = QPushButton("Remove Selected")
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
        from PySide6.QtWidgets import QRadioButton, QButtonGroup

        group = QGroupBox("4. Plot")
        layout = QVBoxLayout(group)

        # X-axis radios (only psiN available for p-files)
        xaxis_row = QHBoxLayout()
        xaxis_row.addWidget(QLabel("X-axis"))
        xaxis_row.addStretch(2)
        self._xaxis_group = QButtonGroup()
        for i, (key, label) in enumerate([
            ('psiN', '\u03c8\u2099'),                    # ψₙ
            ('rho_pol', '\u03c1\u209a\u2092\u2097'),     # ρₚₒₗ
            ('rho_tor', '\u03c1\u209c\u2092\u1d63'),     # ρₜₒᵣ
        ]):
            rb = QRadioButton(label)
            rb.setProperty('xaxis_key', key)
            self._xaxis_group.addButton(rb, i)
            xaxis_row.addWidget(rb, 1)
            if key == 'psiN':
                rb.setChecked(True)
            else:
                rb.setEnabled(False)
        layout.addLayout(xaxis_row)

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
    def _cache_key(shot, suffix, time_s=None):
        if time_s is not None:
            ms = int(round(time_s * 1000))
            return f"{shot}_{ms:06d}_{suffix}"
        return f"{shot}_{suffix}"

    def _open_pfile(self):
        start_dir = getattr(self, '_last_pfile_dir', '')
        if not start_dir:
            start_dir = os.path.expanduser("~")
        dlg = QFileDialog(self.frame, "Open p-File", start_dir,
                          "p-files (p*);;All files (*)")
        dlg.setFileMode(QFileDialog.ExistingFiles)
        dlg.setOption(QFileDialog.DontUseNativeDialog, True)
        dlg.setWindowModality(Qt.NonModal)
        dlg.filesSelected.connect(self._on_pfiles_selected)
        dlg.show()

    def _on_pfiles_selected(self, paths):
        if not paths:
            return
        self._last_pfile_dir = os.path.dirname(paths[0])

        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            for path in paths:
                pdata = _load_pfile(path)
                key = self._cache_key(pdata['shot'], pdata['tree'], pdata['time'][0])
                if key in self._pfile_cache:
                    continue
                self._pfile_cache[key] = pdata
                self._current_key = key
        except Exception as e:
            self.status_label.setStyleSheet(
                "color: red; font-weight: bold; font-size: 9pt;")
            self.status_label.setText(f"p-file error: {e}")
        finally:
            QApplication.restoreOverrideCursor()

        self._rebuild_var_index()
        self._update_var_combo()
        self._update_status('green')
        self._populate_selected()
        self._select_group.setEnabled(bool(self._pfile_cache))
        self._var_group.setEnabled(bool(self._pfile_cache))

        if self._pfile_cache:
            self.preview_button.setEnabled(True)
            self.preview_button.setText(
                f"Browse ({len(self._pfile_cache)} loaded)")

    def _update_status(self, color='green'):
        if not self._pfile_cache:
            self.status_label.setStyleSheet(
                "color: #888; font-weight: bold; font-size: 9pt;")
            self.status_label.setText("No p-file loaded")
            return
        n_runs = len(self._pfile_cache)
        n_v = max(len(p['profiles']) for p in self._pfile_cache.values())
        self.status_label.setStyleSheet(
            f"color: {color}; font-weight: bold; font-size: 9pt;")
        self.status_label.setText(
            f"{n_runs} p-file(s) loaded  \u2014  {n_v} variables")

    def _rebuild_var_index(self):
        """Rebuild intersection of profile variable names across all cached datasets."""
        self._all_profile_vars.clear()
        self._all_units.clear()
        var_sets = [set(pd['profiles'].keys()) for pd in self._pfile_cache.values()]
        common = set.intersection(*var_sets) if var_sets else set()
        first = next(iter(self._pfile_cache.values()), {})
        for vn in common:
            unit = first.get('units', {}).get(vn, '')
            display = f"{vn} [{unit}]" if unit else vn
            self._all_profile_vars[vn] = display
            self._all_units[vn] = unit

    def _update_var_combo(self):
        """Refresh combo items based on filter text."""
        filt = self.var_filter.text().strip().lower()
        current_var = self.var_combo.currentData()

        self.var_combo.blockSignals(True)
        self.var_combo.clear()
        for vn in sorted(self._all_profile_vars):
            display = self._all_profile_vars[vn]
            if filt and filt not in vn.lower() and filt not in display.lower():
                continue
            self.var_combo.addItem(display, vn)

        # Restore previous selection or pending from settings
        restore_var = current_var or getattr(self, '_pending_variable', None)
        if restore_var:
            for i in range(self.var_combo.count()):
                if self.var_combo.itemData(i) == restore_var:
                    self.var_combo.setCurrentIndex(i)
                    break
            self._pending_variable = None
        self.var_combo.blockSignals(False)

    def _populate_selected(self):
        """Auto-add newly loaded p-files to selected list."""
        existing = {self.selected_listbox.item(i).data(Qt.UserRole)
                    for i in range(self.selected_listbox.count())}
        for key, pdata in self._pfile_cache.items():
            if key not in existing:
                shot = pdata['shot']
                suffix = pdata['tree']
                t_ms = int(round(pdata['time'][0] * 1000))
                from PySide6.QtWidgets import QListWidgetItem
                item = QListWidgetItem(f"#{shot} {t_ms:06d}ms ({suffix})")
                item.setData(Qt.UserRole, key)
                self.selected_listbox.addItem(item)

    def _open_preview(self):
        """Browse p-file profiles — merge all cached data."""
        if not self._pfile_cache:
            return

        all_times = []
        all_shots = []
        all_profiles = {}
        psin = None
        all_units = {}
        for pdata in sorted(self._pfile_cache.values(), key=lambda p: (str(p['shot']), p['time'][0])):
            if psin is None:
                psin = pdata.get('psin', np.linspace(0, 1, 100))
                all_units = pdata.get('units', {})
            all_times.append(pdata['time'])
            all_shots.extend([str(pdata['shot'])] * len(pdata['time']))
            for vn, data_2d in pdata.get('profiles', {}).items():
                all_profiles.setdefault(vn, []).append(data_2d)

        merged_time = np.concatenate(all_times)
        sort_idx = np.arange(len(merged_time))  # already sorted

        n_runs = len(self._pfile_cache)
        shots = sorted(set(str(p['shot']) for p in self._pfile_cache.values()))
        shot = ', '.join(shots)
        suffix = f"{n_runs} files"

        pseudo_cdf = {'label': shot, 'profiles': {},
                      'shot_labels': [f"#{s}" for s in all_shots]}
        for vn, data_list in all_profiles.items():
            merged = np.concatenate(data_list, axis=0)[sort_idx]
            unit = all_units.get(vn, '')
            pseudo_cdf['profiles'][vn] = {
                'long_name': f"{vn} [{unit}]" if unit else vn,
                'units': unit,
                'time': merged_time,
                'x': psin,
                'x_dim': 'PSIN',
                'data': merged,
            }

        if not pseudo_cdf['profiles']:
            return

        from ui.widgets.preview_dialog import TranspProfilePreviewDialog
        dlg = TranspProfilePreviewDialog(
            self.frame, pseudo_cdf, shot, suffix,
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
            key = item.data(Qt.UserRole)
            self.selected_listbox.takeItem(self.selected_listbox.row(item))
            self._pfile_cache.pop(key, None)
        self._rebuild_var_index()
        self._update_var_combo()
        if not self._pfile_cache:
            self._select_group.setEnabled(False)
            self._var_group.setEnabled(False)

    def _parse_entry(self, entry):
        """Parse '23879_011700 (kin_0)' -> (cache_key, time_s, shot, suffix, ms)."""
        m = re.match(r'^(\d+)_(\d+)\s*\((.+)\)$', entry)
        shot = m.group(1)
        ms = int(m.group(2))
        suffix = m.group(3)
        cache_key = self._cache_key(shot, suffix)
        time_s = ms / 1000.0
        return cache_key, time_s, shot, suffix, ms

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
        keys = [self.selected_listbox.item(i).data(Qt.UserRole)
                for i in range(self.selected_listbox.count())]
        entries = [self.selected_listbox.item(i).text()
                   for i in range(self.selected_listbox.count())]
        if not keys:
            return
        var_name = self.var_combo.currentData()
        if not var_name:
            return

        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)

        display_text = self.var_combo.currentText()
        ax.set_ylabel(display_text, fontsize=self.label_fontsize)
        ax.set_xlabel(r'$\psi_N$', fontsize=self.label_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(ls='--', lw=0.3, color='#444444')

        colors = self._get_plot_colors(len(keys))

        for idx, key in enumerate(keys):
            try:
                pdata = self._pfile_cache.get(key)
                if pdata is None or var_name not in pdata['profiles']:
                    continue
                shot = pdata['shot']
                suffix = pdata['tree']
                t_ms = int(round(pdata['time'][0] * 1000))
                prof = pdata['profiles'][var_name]
                x_plot = pdata['psin']
                ax.plot(x_plot, prof[0, :],
                        '-', color=colors[idx], lw=1.5, marker='.', markersize=4,
                        label=f'#{shot} {t_ms:06d}ms ({suffix})')
            except Exception as e:
                print(f"[EFIT] Error plotting {key}: {e}")

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
        entries = [self.selected_listbox.item(i).text()
                   for i in range(self.selected_listbox.count())]
        var_name = self.var_combo.currentData()
        if not entries or not var_name:
            return None

        units = self._all_units.get(var_name, '')

        lines = []
        lines.append(f"# p-File Profile: {var_name}\n")
        lines.append(f"#%10s,%10s,%10s,%10s,%14s,%14s\n" % (
            "Shot", "Suffix", "Time[s]", "psi_N",
            f"{var_name}[{units}]", f"d_{var_name}/dpsiN"))

        for entry in entries:
            cache_key, time_s, shot, suffix, ms = self._parse_entry(entry)
            pdata = self._pfile_cache.get(cache_key)
            if pdata is None or var_name not in pdata['profiles']:
                continue
            prof = pdata['profiles'][var_name]
            t_idx = np.argmin(np.abs(pdata['time'] - time_s))
            t_idx = min(t_idx, prof.shape[0] - 1)
            t_actual = pdata['time'][t_idx]
            y = prof[t_idx, :]

            deriv = pdata.get('derivatives', {}).get(var_name)
            if deriv is not None:
                d_idx = min(t_idx, deriv.shape[0] - 1)
                dy = deriv[d_idx, :]
            else:
                dy = np.zeros_like(y)

            psin = pdata['psin']
            for psi, yi, dyi in zip(psin, y, dy):
                lines.append(" %10s,%10s,%10.4f,%10.6f,%14.6e,%14.6e\n" % (
                    shot, suffix, t_actual, psi, yi, dyi))

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
        dlg.setWindowTitle(f"Data Preview \u2014 {title_line}" if title_line else "Data Preview")
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

        # Ctrl+C for selected cells
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
