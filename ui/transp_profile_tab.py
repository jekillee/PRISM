"""
TRANSP CDF Profile tab - load CDF files and plot profile data at selected time points.
Supports multiple CDFs for comparison.
"""

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QScrollArea, QLayout,
    QGroupBox, QLabel, QLineEdit, QPushButton, QFrame,
    QListWidget, QAbstractItemView, QComboBox, QFileDialog,
    QApplication, QStyle, QSpinBox,
    QDialog, QDialogButtonBox,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QShortcut, QKeySequence

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style, get_icon
from ui.theme import ThemeManager


class TranspProfileTab:
    """TRANSP CDF Profile tab: load CDF → select variable + times → plot."""

    def __init__(self, main_window):
        self.main = main_window
        self._cdf_cache = {}  # {label: cdf_data}
        self._all_profile_vars = {}  # {var_name: display_text}

        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10
        self._settings_key = "transp_profile"

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
        self._create_plot_controls(control_frame)
        control_layout.addStretch()

    # ---- 1. Load CDF ----

    def _create_load_section(self, parent):
        group = QGroupBox("1. Load CDF")
        layout = QVBoxLayout(group)

        open_btn = QPushButton("Open CDF...")
        open_btn.clicked.connect(self._open_cdf)
        layout.addWidget(open_btn)

        parent.layout().addWidget(group)

    # ---- 2. Select Data ----

    def _create_select_section(self, parent):
        group = QGroupBox("2. Select Data")
        group.setEnabled(False)
        self._select_group = group
        layout = QVBoxLayout(group)

        # Loaded CDF status
        self.cdf_status = QLabel("No CDF loaded")
        self.cdf_status.setWordWrap(True)
        layout.addWidget(self.cdf_status)

        # Variable filter + combo
        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Filter"))
        self.var_filter = QLineEdit()
        self.var_filter.setPlaceholderText("Type to filter variables...")
        self.var_filter.textChanged.connect(self._update_var_combo)
        filter_row.addWidget(self.var_filter, 1)
        layout.addLayout(filter_row)

        var_row = QHBoxLayout()
        var_row.addWidget(QLabel("Variable"))
        self.var_combo = QComboBox()
        self.var_combo.setMaxVisibleItems(20)
        self.var_combo.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.var_combo.setMinimumContentsLength(10)
        self.var_combo.setStyleSheet("QComboBox { combobox-popup: 0; }")
        self.var_combo.view().setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.var_combo.currentIndexChanged.connect(self._on_var_changed)
        var_row.addWidget(self.var_combo, 1)
        layout.addLayout(var_row)

        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setFrameShadow(QFrame.Sunken)
        layout.addWidget(sep)

        # Available / Selected time lists
        lists_row = QHBoxLayout()

        avail_col = QVBoxLayout()
        avail_col.addWidget(QLabel("Time [s]"))
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
                          "CDF files (*.CDF *.cdf *.nc);;All files (*)")
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

        self.available_listbox.clear()

        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            from data_loaders.transp_cdf_loader import load_transp_cdf
            cdf = load_transp_cdf(path)
            self._cdf_cache[cdf['label']] = cdf
            self._current_label = cdf['label']
        finally:
            QApplication.restoreOverrideCursor()

        self._rebuild_var_index()
        self._update_var_combo()
        self._update_cdf_status()
        self._select_group.setEnabled(bool(self._cdf_cache))

    def _update_cdf_status(self, color='green'):
        if not self._cdf_cache:
            self.cdf_status.setStyleSheet("color: #888; font-weight: bold; font-size: 9pt;")
            self.cdf_status.setText("No CDF loaded")
            return
        label = getattr(self, '_current_label', list(self._cdf_cache.keys())[-1])
        cdf = self._cdf_cache[label]
        shot, run = self._parse_run_label(label)
        n_p = len(cdf['profiles'])
        n_t = len(cdf['timetraces'])
        self.cdf_status.setStyleSheet(f"color: {color}; font-weight: bold; font-size: 9pt;")
        self.cdf_status.setText(
            f"#{shot} ({run})  —  {n_p} profiles, {n_t} time traces")

    def _rebuild_var_index(self):
        """Rebuild union of profile variable names across all CDFs."""
        self._all_profile_vars.clear()
        for cdf in self._cdf_cache.values():
            for vn, vd in cdf['profiles'].items():
                if vn not in self._all_profile_vars:
                    parts = [vn]
                    if vd['long_name'] and vd['long_name'] != vn:
                        parts.append(f"- {vd['long_name']}")
                    if vd['units']:
                        parts.append(f"[{vd['units']}]")
                    self._all_profile_vars[vn] = ' '.join(parts)

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

        # Restore previous selection
        if current_var:
            for i in range(self.var_combo.count()):
                if self.var_combo.itemData(i) == current_var:
                    self.var_combo.setCurrentIndex(i)
                    break
        self.var_combo.blockSignals(False)
        self._on_var_changed()

    def _on_var_changed(self):
        """Populate available time list for selected variable (current CDF only)."""
        self.available_listbox.clear()
        var_name = self.var_combo.currentData()
        if not var_name:
            return

        label = getattr(self, '_current_label', None)
        if label is None or label not in self._cdf_cache:
            return
        cdf = self._cdf_cache[label]
        if var_name not in cdf['profiles']:
            return
        shot, run = self._parse_run_label(label)
        prof = cdf['profiles'][var_name]
        for t in prof['time']:
            ms = int(round(t * 1000))
            self.available_listbox.addItem(f"{shot}_{ms:06d} ({run})")

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
        """Parse '39551_001500 (W09)' → (cdf_label, time_s, shot, run, ms)."""
        import re
        m = re.match(r'^(\d+)_(\d+)\s*\((.+)\)$', entry)
        shot = m.group(1)
        ms = int(m.group(2))
        run = m.group(3)
        cdf_label = f"{shot}{run}"
        time_s = ms / 1000.0
        return cdf_label, time_s, shot, run, ms

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

        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)

        # Axis labels
        display_text = self.var_combo.currentText()
        ax.set_ylabel(display_text, fontsize=self.label_fontsize)

        x_label = 'X'
        for cdf in self._cdf_cache.values():
            if var_name in cdf['profiles']:
                x_label = cdf['profiles'][var_name].get('x_label', 'X')
                break
        ax.set_xlabel(x_label, fontsize=self.label_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.grid(ls='--', lw=0.3, color='#444444')

        colors = self._get_plot_colors(len(entries))

        for idx, entry in enumerate(entries):
            try:
                cdf_label, time_s, shot, run, ms = self._parse_entry(entry)
                cdf = self._cdf_cache.get(cdf_label)
                if cdf is None or var_name not in cdf['profiles']:
                    continue
                prof = cdf['profiles'][var_name]
                t_idx = np.argmin(np.abs(prof['time'] - time_s))
                ax.plot(prof['x'], prof['data'][t_idx, :],
                        '-', color=colors[idx], lw=1.5,
                        label=f'#{shot} {ms:06d}ms ({run})')
            except Exception as e:
                print(f"[TRANSP] Error plotting {entry}: {e}")

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

    def save_settings(self):
        from config.user_settings import get_tab_settings, set_tab_settings
        s = get_tab_settings(self._settings_key)
        s["color_mode"] = getattr(self, '_color_mode', 'Gradient(viridis)')
        s["label_fontsize"] = self.label_fontsize
        s["legend_fontsize"] = self.legend_fontsize
        s["tick_fontsize"] = self.tick_fontsize
        set_tab_settings(self._settings_key, s)
