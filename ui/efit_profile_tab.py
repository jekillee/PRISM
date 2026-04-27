"""
EFIT Profile tab - load EFIT GEQDSK 1D profiles from MDS+ or g-files.
Plots PRES, PPRIME, FPOL, FFPRIM, QPSI, RHOVN at selected time points.
Supports multiple shots/trees for comparison.
"""

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


# Profile variables available in GEQDSK / EFIT MDS+ trees
_EFIT_PROFILE_VARS = {
    'PRES':   'PRES - Pressure [Pa]',
    'PPRIME': 'PPRIME - dP/d\u03c8 [Pa/(Wb/rad)]',
    'FPOL':   'FPOL - Poloidal current F=RB_t [T\u00b7m]',
    'FFPRIM': 'FFPRIM - FF\' [T\u00b2\u00b7m\u00b2/(Wb/rad)]',
    'QPSI':   'QPSI - Safety factor q',
    'RHOVN':  'RHOVN - Normalized toroidal flux \u03c1_tor',
}

_EFIT_PROFILE_UNITS = {
    'PRES': 'Pa', 'PPRIME': 'Pa/(Wb/rad)',
    'FPOL': 'T*m', 'FFPRIM': 'T^2*m^2/(Wb/rad)',
    'QPSI': '', 'RHOVN': '',
}


def _load_efit_mds(shot, tree, mds_ip):
    """Load EFIT 1D profile data from MDS+ tree.

    Returns dict with keys:
        shot, tree, time (1D, seconds), psin (1D grid),
        profiles: {var_name: 2D array [nt, nw]}, ...
    """
    from MDSplus import Connection

    mds = Connection(mds_ip)
    mds.openTree(tree, int(shot))
    try:
        time_arr = np.array(mds.get('\\gtime').data(), dtype=float)
        # efitrt1/efitrt2 store time in seconds; others in ms
        if tree not in ('efitrt1', 'efitrt2'):
            time_arr = time_arr / 1e3

        nw = int(mds.get('\\mw').data()[0])
        psin = np.linspace(0.0, 1.0, nw)

        profiles = {}
        for node, var in [('\\pres', 'PRES'), ('\\pprime', 'PPRIME'),
                          ('\\fpol', 'FPOL'), ('\\ffprim', 'FFPRIM'),
                          ('\\qpsi', 'QPSI'), ('\\rhovn', 'RHOVN')]:
            try:
                d = np.array(mds.get(node).data(), dtype=float)
                if d.ndim == 1:
                    d = d.reshape(1, -1)
                profiles[var] = d
            except Exception:
                pass
    finally:
        mds.closeTree(tree, int(shot))

    return {
        'shot': str(shot), 'tree': tree,
        'time': time_arr, 'psin': psin,
        'profiles': profiles,
    }


def _load_gfile(path):
    """Read a GEQDSK g-file and return the same dict format as _load_efit_mds.

    Supports standard GEQDSK format (fixed-width Fortran output).
    """
    import os
    print(f"[EFIT] Loading g-file: {os.path.basename(path)}")

    with open(path, 'r') as f:
        lines = f.readlines()

    # ---- Header line ----
    header = lines[0]
    # Last 3 integers on header line: idum, nw, nh
    tokens = header.split()
    nw = int(tokens[-2])
    nh = int(tokens[-1])

    def _read_1d(start_line, count):
        """Read *count* floats starting from line index *start_line*,
        5 values per line in 16-char fixed-width fields."""
        vals = []
        li = start_line
        while len(vals) < count:
            line = lines[li]
            # Each value occupies 16 chars
            n_on_line = len(line.rstrip('\n')) // 16
            for k in range(n_on_line):
                vals.append(float(line[16*k:16*(k+1)]))
                if len(vals) == count:
                    break
            li += 1
        return np.array(vals), li

    def _n_lines(count):
        return (count + 4) // 5  # 5 values per line

    # Line 1: header (already read)
    # Line 2: 4 values — rdim, zdim, rcentr, rleft, zmid
    # Actually line index 1 has 5 scalars
    scalar_line = 1
    # Read 20 scalar values (lines 1..4, 5 per line)
    scalars, next_li = _read_1d(1, 20)
    rdim, zdim, rcentr, rleft, zmid = scalars[0:5]
    rmaxis, zmaxis, simag_val, sibry_val, bcentr = scalars[5:10]
    current, simag2, _, rmaxis2, _ = scalars[10:15]
    zmaxis2, _, sibry2, _, _ = scalars[15:20]

    # 1D arrays: fpol, pres, ffprim, pprime — each nw long
    fpol, next_li = _read_1d(next_li, nw)
    pres, next_li = _read_1d(next_li, nw)  # this is actually ffprim in some formats
    # Standard order: fpol, pres, ffprim, pprime
    # Re-read properly:
    # After 20 scalars come: FPOL(nw), PRES(nw), FFPRIM(nw), PPRIME(nw), PSIRZ(nw*nh), QPSI(nw)
    # Let's re-read from scratch after scalars
    li = 1 + _n_lines(20)
    fpol_arr, li = _read_1d(li, nw)
    pres_arr, li = _read_1d(li, nw)
    ffprim_arr, li = _read_1d(li, nw)
    pprime_arr, li = _read_1d(li, nw)
    psirz_flat, li = _read_1d(li, nw * nh)
    qpsi_arr, li = _read_1d(li, nw)

    psin = np.linspace(0.0, 1.0, nw)

    # Compute RHOVN from q profile (toroidal flux integration)
    phi = np.zeros(nw)
    for i in range(nw - 1):
        phi[i + 1] = np.trapz(qpsi_arr[:i + 2], x=psin[:i + 2])
    if phi[-1] != 0:
        rhovn_arr = np.sqrt(phi / phi[-1])
    else:
        rhovn_arr = psin.copy()

    profiles = {
        'PRES': pres_arr.reshape(1, -1),
        'PPRIME': pprime_arr.reshape(1, -1),
        'FPOL': fpol_arr.reshape(1, -1),
        'FFPRIM': ffprim_arr.reshape(1, -1),
        'QPSI': qpsi_arr.reshape(1, -1),
        'RHOVN': rhovn_arr.reshape(1, -1),
    }

    # Try to extract shot/time/suffix from filename  e.g. g040848.003000_kin_0
    basename = os.path.basename(path)
    m = re.match(r'[ga](\d+)\.(\d+)(?:_(.+))?', basename)
    if m:
        shot_str = str(int(m.group(1)))
        time_ms = int(m.group(2))
        time_s = time_ms / 1000.0
        tree = m.group(3) if m.group(3) else 'gfile'
    else:
        shot_str = basename
        time_s = 0.0
        tree = 'gfile'
    return {
        'shot': shot_str, 'tree': tree,
        'time': np.array([time_s]),
        'psin': psin,
        'profiles': profiles,
    }


class EfitProfileTab:
    """EFIT Profile tab: load EFIT MDS+ / g-file -> select times -> plot 1D profiles."""

    def __init__(self, main_window):
        self.main = main_window
        self._efit_cache = {}          # {cache_key: efit_data}
        self._current_key = None       # cache key of last loaded dataset
        self._all_profile_vars = {}    # {var_name: display_text}

        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10
        self._settings_key = "efit_profile"

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
        from config.app_config import AppConfig
        group = QGroupBox("1. Load Data")
        grid = QGridLayout(group)
        grid.setColumnStretch(1, 1)

        # Shot entry
        grid.addWidget(QLabel("Shot"), 0, 0)
        self.shot_entry = QLineEdit()
        self.shot_entry.returnPressed.connect(self._fetch)
        grid.addWidget(self.shot_entry, 0, 1)

        # Shot up/down buttons
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

        # EFIT tree dropdown (between up/down and Fetch)
        self.tree_combo = QComboBox()
        cfg = AppConfig()
        for display, tree in cfg.EFIT_TREES.items():
            self.tree_combo.addItem(tree, tree)
        grid.addWidget(self.tree_combo, 0, 3)

        # Fetch button
        self.fetch_button = QPushButton("Fetch")
        self.fetch_button.setFixedWidth(70)
        self.fetch_button.clicked.connect(self._fetch)
        grid.addWidget(self.fetch_button, 0, 4)

        # Open g-file button row with "Or" prefix
        or_label = QLabel("Or")
        grid.addWidget(or_label, 1, 0)
        gfile_btn = QPushButton("Open g-File...")
        gfile_btn.clicked.connect(self._open_gfile)
        grid.addWidget(gfile_btn, 1, 1, 1, 4)

        # Status label
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        grid.addWidget(self.status_label, 2, 0, 1, 5)

        parent.layout().addWidget(group)

    # ---- 2. Select Data ----

    def _create_select_section(self, parent):
        group = QGroupBox("2. Select Data")
        group.setEnabled(False)
        self._select_group = group
        layout = QVBoxLayout(group)

        # Preview / Browse button
        self.preview_button = QPushButton("Load data to browse")
        self.preview_button.setEnabled(False)
        self.preview_button.clicked.connect(self._open_preview)
        layout.addWidget(self.preview_button)

        # Available / Selected time lists
        lists_row = QHBoxLayout()

        avail_col = QVBoxLayout()
        avail_col.addWidget(QLabel("Available"))
        self.available_listbox = QListWidget()
        self.available_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.available_listbox.setFixedHeight(180)
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
        self.selected_listbox.setFixedHeight(180)
        sel_col.addWidget(self.selected_listbox)
        lists_row.addLayout(sel_col, stretch=1)

        layout.addLayout(lists_row)

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

        # X-axis radios (psiN default, R disabled for EFIT profiles)
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

    def _adjust_shot(self, delta):
        txt = self.shot_entry.text().strip()
        try:
            self.shot_entry.setText(str(int(txt) + delta))
        except ValueError:
            pass

    @staticmethod
    def _cache_key(shot, tree, time_s=None):
        if time_s is not None:
            ms = int(round(time_s * 1000))
            return f"{shot}_{ms:06d}_{tree}"
        return f"{shot}_{tree}"

    def _fetch(self):
        shot = self.shot_entry.text().strip()
        if not shot:
            return
        tree = self.tree_combo.currentData()
        key = self._cache_key(shot, tree)

        if key in self._efit_cache:
            self._current_key = key
        else:
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                from config.app_config import AppConfig
                cfg = AppConfig()
                efit = _load_efit_mds(shot, tree, cfg.MDS_IP)
                self._efit_cache[key] = efit
                self._current_key = key
            except Exception as e:
                QApplication.restoreOverrideCursor()
                self.status_label.setStyleSheet(
                    "color: red; font-weight: bold; font-size: 9pt;")
                self.status_label.setText(str(e))
                return
            finally:
                QApplication.restoreOverrideCursor()

        self._rebuild_var_index()
        self._update_var_combo()
        self._update_status()
        self._populate_available()
        self._select_group.setEnabled(bool(self._efit_cache))
        self._var_group.setEnabled(bool(self._efit_cache))

        if self._efit_cache:
            efit = self._efit_cache[self._current_key]
            self.preview_button.setEnabled(True)
            self.preview_button.setText(
                f"Browse #{efit['shot']} ({efit['tree']})")

    def _open_gfile(self):
        start_dir = getattr(self, '_last_gfile_dir', '')
        if not start_dir:
            import os
            start_dir = os.path.expanduser("~")
        dlg = QFileDialog(self.frame, "Open GEQDSK g-File", start_dir,
                          "GEQDSK files (*.geq g*);;All files (*)")
        dlg.setFileMode(QFileDialog.ExistingFiles)
        dlg.setOption(QFileDialog.DontUseNativeDialog, True)
        dlg.setWindowModality(Qt.NonModal)
        dlg.filesSelected.connect(self._on_gfiles_selected)
        dlg.show()

    def _on_gfiles_selected(self, paths):
        import os
        if not paths:
            return
        self._last_gfile_dir = os.path.dirname(paths[0])

        new_keys = []
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            for path in paths:
                efit = _load_gfile(path)
                key = self._cache_key(efit['shot'], efit['tree'], efit['time'][0])
                self._efit_cache[key] = efit
                self._current_key = key
                if key not in new_keys:
                    new_keys.append(key)
        except Exception as e:
            self.status_label.setStyleSheet(
                "color: red; font-weight: bold; font-size: 9pt;")
            self.status_label.setText(f"g-file error: {e}")
        finally:
            QApplication.restoreOverrideCursor()

        self._rebuild_var_index()
        self._update_var_combo()
        self._update_status('green')
        # Show only the just-opened g-files in Available list
        self.available_listbox.clear()
        for key in new_keys:
            efit = self._efit_cache[key]
            shot = efit['shot']
            tree = efit['tree']
            for t in efit['time']:
                ms = int(round(t * 1000))
                self.available_listbox.addItem(f"{shot}_{ms:06d} ({tree})")
        self._select_group.setEnabled(True)
        self._var_group.setEnabled(True)

        self.preview_button.setEnabled(True)
        self.preview_button.setText(
            f"Browse ({len(self._efit_cache)} loaded)")

    def _update_status(self, color='green'):
        if not self._efit_cache:
            self.status_label.setStyleSheet(
                "color: #888; font-weight: bold; font-size: 9pt;")
            self.status_label.setText("No EFIT data loaded")
            return
        n_runs = len(self._efit_cache)
        n_v = max(len(e['profiles']) for e in self._efit_cache.values())
        total_t = sum(len(e['time']) for e in self._efit_cache.values())
        self.status_label.setStyleSheet(
            f"color: {color}; font-weight: bold; font-size: 9pt;")
        self.status_label.setText(
            f"{n_runs} run(s) loaded  \u2014  {n_v} profiles, {total_t} time points")

    def _rebuild_var_index(self):
        """Rebuild intersection of profile variable names across all cached datasets."""
        self._all_profile_vars.clear()
        var_sets = [set(efit['profiles'].keys()) for efit in self._efit_cache.values()]
        common = set.intersection(*var_sets) if var_sets else set()
        for vn in common:
            self._all_profile_vars[vn] = _EFIT_PROFILE_VARS.get(vn, vn)

    def _update_var_combo(self):
        """Refresh combo items based on filter text."""
        filt = self.var_filter.text().strip().upper()
        current_var = self.var_combo.currentData()

        self.var_combo.blockSignals(True)
        self.var_combo.clear()
        for vn in sorted(self._all_profile_vars):
            display = self._all_profile_vars[vn]
            if filt and filt not in vn.upper() and filt not in display.upper():
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

    def _populate_available(self):
        """Fill available list with time points from all cached datasets."""
        self.available_listbox.clear()
        for key in sorted(self._efit_cache):
            efit = self._efit_cache[key]
            shot = efit['shot']
            tree = efit['tree']
            for t in efit['time']:
                ms = int(round(t * 1000))
                self.available_listbox.addItem(f"{shot}_{ms:06d} ({tree})")

    def _open_preview(self):
        """Browse EFIT profiles — merge all cached data into one preview."""
        if not self._efit_cache:
            return

        # Merge all cached datasets: concatenate times and profile data
        all_times = []
        all_shots = []
        all_profiles = {}
        psin = None
        for efit in sorted(self._efit_cache.values(), key=lambda e: (str(e['shot']), e['time'][0])):
            if psin is None:
                psin = efit.get('psin', np.linspace(0, 1, 129))
            all_times.append(efit['time'])
            all_shots.extend([str(efit['shot'])] * len(efit['time']))
            for vn, data_2d in efit.get('profiles', {}).items():
                all_profiles.setdefault(vn, []).append(data_2d)

        merged_time = np.concatenate(all_times)
        # Already sorted by shot→time via sorted() above
        sort_idx = np.arange(len(merged_time))

        n_runs = len(self._efit_cache)
        shots = sorted(set(str(e['shot']) for e in self._efit_cache.values()))
        shot = ', '.join(shots)
        tree = f"{n_runs} runs"

        pseudo_cdf = {'label': f"{shot}", 'profiles': {},
                      'shot_labels': [f"#{s}" for s in all_shots]}
        for vn, data_list in all_profiles.items():
            merged = np.concatenate(data_list, axis=0)[sort_idx]
            pseudo_cdf['profiles'][vn] = {
                'long_name': _EFIT_PROFILE_VARS.get(vn, vn),
                'units': '',
                'time': merged_time,
                'x': psin,
                'x_dim': 'PSIN',
                'data': merged,
            }

        if not pseudo_cdf['profiles']:
            return

        from ui.widgets.preview_dialog import TranspProfilePreviewDialog
        dlg = TranspProfilePreviewDialog(
            self.frame, pseudo_cdf, shot, tree,
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

    def _parse_entry(self, entry):
        """Parse '40848_003000 (efit01)' -> (cache_key, time_s, shot, tree, ms)."""
        m = re.match(r'^(\d+)_(\d+)\s*\((.+)\)$', entry)
        shot = m.group(1)
        ms = int(m.group(2))
        tree = m.group(3)
        time_s = ms / 1000.0
        # Try time-based key first (g-file), then without time (MDS+)
        cache_key = self._cache_key(shot, tree, time_s)
        if cache_key not in self._efit_cache:
            cache_key = self._cache_key(shot, tree)
        return cache_key, time_s, shot, tree, ms

    # ================================================================
    # X-axis conversion
    # ================================================================

    @staticmethod
    def _get_x_axis(efit, var_name, t_idx, xaxis_key):
        """Return x-axis array for the requested coordinate.

        EFIT profiles are defined on a uniform psi_N grid [0, 1].
        - psiN: native grid
        - rho_pol: sqrt(psi_N)
        - rho_tor: RHOVN profile data at t_idx (if available)
        """
        psin = efit['psin']

        if xaxis_key == 'psiN':
            return psin

        if xaxis_key == 'rho_pol':
            return np.sqrt(np.clip(psin, 0, None))

        if xaxis_key == 'rho_tor':
            rhovn = efit['profiles'].get('RHOVN')
            if rhovn is not None:
                idx = min(t_idx, rhovn.shape[0] - 1)
                return rhovn[idx, :]
            # Fallback: sqrt(psin)
            return np.sqrt(np.clip(psin, 0, None))

        return psin

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
        entries = [self.selected_listbox.item(i).text()
                   for i in range(self.selected_listbox.count())]
        if not entries:
            return
        var_name = self.var_combo.currentData()
        if not var_name:
            return

        # Determine x-axis mode
        xaxis_key = 'psiN'
        checked = self._xaxis_group.checkedButton()
        if checked:
            xaxis_key = checked.property('xaxis_key')
        _XAXIS_LABELS = {
            'psiN': r'$\psi_N$',
            'rho_pol': r'$\rho_{pol}$',
            'rho_tor': r'$\rho_{tor}$',
        }

        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)

        display_text = self.var_combo.currentText()
        ax.set_ylabel(display_text, fontsize=self.label_fontsize)
        ax.set_xlabel(_XAXIS_LABELS.get(xaxis_key, xaxis_key),
                      fontsize=self.label_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(ls='--', lw=0.3, color='#444444')

        colors = self._get_plot_colors(len(entries))

        for idx, entry in enumerate(entries):
            try:
                cache_key, time_s, shot, tree, ms = self._parse_entry(entry)
                efit = self._efit_cache.get(cache_key)
                if efit is None or var_name not in efit['profiles']:
                    continue
                prof = efit['profiles'][var_name]
                t_idx = np.argmin(np.abs(efit['time'] - time_s))
                t_idx = min(t_idx, prof.shape[0] - 1)
                x_plot = self._get_x_axis(efit, var_name, t_idx, xaxis_key)
                ax.plot(x_plot, prof[t_idx, :],
                        '-', color=colors[idx], lw=1.5, marker='.', markersize=4,
                        label=f'#{shot} {ms:06d}ms ({tree})')
            except Exception as e:
                print(f"[EFIT] Error plotting {entry}: {e}")

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

        units = _EFIT_PROFILE_UNITS.get(var_name, '')

        lines = []
        lines.append(f"# EFIT Profile: {var_name}\n")
        lines.append(f"#%10s,%10s,%10s,%10s,%10s,%10s,%14s\n" % (
            "Shot", "Tree", "Time[s]", "psi_N", "rho_pol", "rho_tor",
            f"{var_name}[{units}]"))

        for entry in entries:
            cache_key, time_s, shot, tree, ms = self._parse_entry(entry)
            efit = self._efit_cache.get(cache_key)
            if efit is None or var_name not in efit['profiles']:
                continue
            prof = efit['profiles'][var_name]
            t_idx = np.argmin(np.abs(efit['time'] - time_s))
            t_idx = min(t_idx, prof.shape[0] - 1)
            t_actual = efit['time'][t_idx]
            y = prof[t_idx, :]

            x_psi = self._get_x_axis(efit, var_name, t_idx, 'psiN')
            x_rpol = self._get_x_axis(efit, var_name, t_idx, 'rho_pol')
            x_rtor = self._get_x_axis(efit, var_name, t_idx, 'rho_tor')

            for psi, rpol, rtor, yi in zip(x_psi, x_rpol, x_rtor, y):
                lines.append(" %10s,%10s,%10.4f,%10.6f,%10.6f,%10.6f,%14.6e\n" % (
                    shot, tree, t_actual, psi, rpol, rtor, yi))

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
        saved_shot = s.get("shot", "")
        if saved_shot:
            self.shot_entry.setText(saved_shot)

    def save_settings(self):
        from config.user_settings import get_tab_settings, set_tab_settings
        s = get_tab_settings(self._settings_key)
        s["color_mode"] = getattr(self, '_color_mode', 'Gradient(viridis)')
        s["label_fontsize"] = self.label_fontsize
        s["legend_fontsize"] = self.legend_fontsize
        s["tick_fontsize"] = self.tick_fontsize
        s["variable"] = self.var_combo.currentData() or ""
        s["shot"] = self.shot_entry.text().strip()
        set_tab_settings(self._settings_key, s)
