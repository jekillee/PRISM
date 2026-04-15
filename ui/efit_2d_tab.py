"""
EFIT 2D tab - load EFIT GEQDSK 2D equilibrium from MDS+ or g-files.
Plots 2D psi_N contours, separatrix, boundary, and magnetic axis.
Supports multiple shots/trees for boundary overlay comparison.
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


class Efit2DTab:
    """EFIT 2D tab: load EFIT MDS+ / g-file -> select times -> plot 2D equilibrium."""

    def __init__(self, main_window):
        self.main = main_window
        self._efit_cache = {}          # {cache_key: efit_data}
        self._current_key = None       # cache key of last loaded dataset

        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10
        self._cmap = 'RdBu_r'
        self._n_contours = 9
        self._contour_type = 'psi_N'
        self._sol_contours = 4
        self._sol_distance = 0.01  # meters
        self._settings_key = "efit_2d"

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

        # Fetch button
        self.fetch_button = QPushButton("Fetch")
        self.fetch_button.setFixedWidth(70)
        self.fetch_button.clicked.connect(self._fetch)
        grid.addWidget(self.fetch_button, 0, 3)

        # EFIT tree dropdown
        grid.addWidget(QLabel("EFIT Tree"), 1, 0)
        self.tree_combo = QComboBox()
        cfg = AppConfig()
        for display, tree in cfg.EFIT_TREES.items():
            self.tree_combo.addItem(display, tree)
        grid.addWidget(self.tree_combo, 1, 1, 1, 3)

        # Open g-file button
        gfile_btn = QPushButton("Open g-File...")
        gfile_btn.clicked.connect(self._open_gfile)
        grid.addWidget(gfile_btn, 2, 0, 1, 4)

        # Status label
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        grid.addWidget(self.status_label, 3, 0, 1, 4)

        parent.layout().addWidget(group)

    # ---- 2. Select Data ----

    def _create_select_section(self, parent):
        group = QGroupBox("2. Select Data")
        group.setEnabled(False)
        self._select_group = group
        layout = QVBoxLayout(group)

        # Browse button
        self.browse_button = QPushButton("Load data to browse")
        self.browse_button.setEnabled(False)
        self.browse_button.clicked.connect(self._open_browse)
        layout.addWidget(self.browse_button)

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

    # ---- 3. Plot ----

    def _create_plot_controls(self, parent):
        from PySide6.QtWidgets import QRadioButton, QButtonGroup

        group = QGroupBox("3. Plot")
        layout = QVBoxLayout(group)

        # Contour type radios (same layout as PRISM-Profiles X-axis)
        ct_row = QHBoxLayout()
        ct_row.addWidget(QLabel("Contour"))
        ct_row.addStretch(2)
        self._ct_group = QButtonGroup()
        for i, (key, label) in enumerate([
            ('psi_N', '\u03c8\u2099'),                                 # ψₙ
            ('rho_pol', '\u03c1\u209a\u2092\u2097'),                  # ρₚₒₗ
            ('rho_tor', '\u03c1\u209c\u2092\u1d63'),                  # ρₜₒᵣ
        ]):
            rb = QRadioButton(label)
            rb.setProperty('ct_key', key)
            self._ct_group.addButton(rb, i)
            ct_row.addWidget(rb, 1)
            if key == self._contour_type:
                rb.setChecked(True)
        layout.addLayout(ct_row)

        btn_row = QHBoxLayout()
        plot_btn = QPushButton("Plot")
        plot_btn.clicked.connect(self._plot)
        btn_row.addWidget(plot_btn, 3)
        opt_btn = QPushButton("Option")
        opt_btn.clicked.connect(self._show_style_dialog)
        btn_row.addWidget(opt_btn, 1)
        layout.addLayout(btn_row)

        parent.layout().addWidget(group)

    # ---- 4. Save Data ----

    def _create_save_section(self, parent):
        group = QGroupBox("4. Save Data")
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
            self._efit_cache.clear()
            self.available_listbox.clear()
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                from data_loaders.efit_viewer_loader import load_efit_mds
                efit_data = load_efit_mds(int(shot), tree, load_scalars=False, load_profiles=False)
                self._efit_cache[key] = efit_data
                self._current_key = key
            except Exception as e:
                QApplication.restoreOverrideCursor()
                self.status_label.setStyleSheet(
                    "color: red; font-weight: bold; font-size: 9pt;")
                self.status_label.setText(str(e))
                return
            finally:
                QApplication.restoreOverrideCursor()

        self._update_status()
        self._populate_available()
        self._select_group.setEnabled(bool(self._efit_cache))
        if self._current_key:
            shot, tree = self._current_key.split('_', 1)
            self.browse_button.setEnabled(True)
            self.browse_button.setText(f"Browse #{shot} ({tree})")

    def _open_gfile(self):
        start_dir = getattr(self, '_last_gfile_dir', '')
        if not start_dir:
            import os
            start_dir = os.path.expanduser("~")
        dlg = QFileDialog(self.frame, "Open GEQDSK g-File", start_dir,
                          "GEQDSK (g*);;All files (*)")
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

        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            from data_loaders.efit_viewer_loader import load_gfile
            for path in paths:
                efit_data = load_gfile(path)
                shot = efit_data['shot']
                tree = efit_data.get('tree') or 'gfile'
                if shot is None:
                    shot = os.path.splitext(os.path.basename(path))[0]
                eq2d = efit_data.get('eq2d', {})
                t0 = eq2d['time'][0] if eq2d.get('time') is not None and len(eq2d.get('time', [])) > 0 else None
                key = self._cache_key(shot, tree, t0)
                if key in self._efit_cache:
                    continue
                self._efit_cache[key] = efit_data
                self._current_key = key
        except Exception as e:
            self.status_label.setStyleSheet(
                "color: red; font-weight: bold; font-size: 9pt;")
            self.status_label.setText(f"g-file error: {e}")
        finally:
            QApplication.restoreOverrideCursor()

        self._update_status('green')
        self._populate_available()
        self._select_group.setEnabled(bool(self._efit_cache))
        if self._current_key:
            shot, tree = self._current_key.split('_', 1)
            self.browse_button.setEnabled(True)
            self.browse_button.setText(f"Browse #{shot} ({tree})")

    def _update_status(self, color='green'):
        if not self._efit_cache:
            self.status_label.setStyleSheet(
                "color: #888; font-weight: bold; font-size: 9pt;")
            self.status_label.setText("No EFIT data loaded")
            return
        n_runs = len(self._efit_cache)
        total_t = sum(
            len(e.get('eq2d', {}).get('time', []))
            for e in self._efit_cache.values())
        self.status_label.setStyleSheet(
            f"color: {color}; font-weight: bold; font-size: 9pt;")
        self.status_label.setText(
            f"{n_runs} run(s) loaded  \u2014  {total_t} time points")

    def _populate_available(self):
        """Fill available list with time points from all cached datasets."""
        self.available_listbox.clear()
        for key in sorted(self._efit_cache):
            efit = self._efit_cache[key]
            eq2d = efit.get('eq2d', {})
            time_arr = eq2d.get('time')
            if time_arr is None:
                continue
            shot = efit.get('shot', '?')
            tree = efit.get('tree', '?')
            for t in time_arr:
                ms = int(round(t * 1000))
                self.available_listbox.addItem(f"{shot}_{ms:06d} ({tree})")

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
        m = re.match(r'^(.+?)_(\d+)\s*\((.+)\)$', entry)
        shot = m.group(1)
        ms = int(m.group(2))
        tree = m.group(3)
        time_s = ms / 1000.0
        cache_key = self._cache_key(shot, tree, time_s)
        if cache_key not in self._efit_cache:
            cache_key = self._cache_key(shot, tree)
        return cache_key, time_s, shot, tree, ms

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

    def _open_browse(self):
        """Open a slider-based 2D equilibrium browser — merges all cached data."""
        if not self._efit_cache:
            return

        # Merge all cached eq2d into a list of (efit_data, t_idx, time_s) tuples
        merged_slices = []  # [(efit_data, t_idx, time_s, shot, tree)]
        for key, efit in self._efit_cache.items():
            eq2d = efit.get('eq2d', {})
            time_arr = eq2d.get('time')
            if time_arr is None or 'psirz' not in eq2d:
                continue
            shot_str = str(efit.get('shot', '?'))
            tree_str = str(efit.get('tree', '?'))
            for ti, t in enumerate(time_arr):
                merged_slices.append((efit, ti, t, shot_str, tree_str))

        if not merged_slices:
            return
        # Sort by shot, then time
        merged_slices.sort(key=lambda x: (x[3], x[2]))

        from ui.widgets.preview_dialog import _PreviewBase

        n_times = len(merged_slices)
        time_arr = np.array([s[2] for s in merged_slices])
        shot = merged_slices[0][3]
        tree = f"{len(self._efit_cache)} runs"

        class _Efit2DBrowse(_PreviewBase):
            def __init__(s, parent, slices, shot, tree, selected_listbox):
                s.slices = slices  # [(efit_data, t_idx, time_s, shot, tree)]
                s.shot = shot
                s.tree = tree
                s._n_contours = self._n_contours
                s._contour_type = self._contour_type
                s._sol_contours = self._sol_contours
                s._cmap = self._cmap
                super().__init__(parent, n_times,
                                 f"Browse #{shot} ({tree})",
                                 selected_listbox)
                s._build_ui()
                s._init_plot()
                s._update_plot(0)

            def _build_ui(s):
                from PySide6.QtWidgets import QHBoxLayout, QLabel, QLineEdit, QPushButton, QGroupBox, QVBoxLayout, QStyle
                main_layout = QHBoxLayout(s)
                main_layout.setContentsMargins(4, 4, 4, 4)
                main_layout.setSpacing(4)
                s._build_left(main_layout, "Use mouse wheel or arrow keys to navigate time.")
                def _nav(layout):
                    row1 = QHBoxLayout()
                    row1.addWidget(QLabel("Frame"))
                    s.frame_entry = QLineEdit("1")
                    s.frame_entry.returnPressed.connect(s._goto_frame)
                    row1.addWidget(s.frame_entry)
                    s.frame_max_label = QLabel(f"/ {n_times}")
                    row1.addWidget(s.frame_max_label)
                    go = QPushButton()
                    go.setIcon(s.style().standardIcon(QStyle.SP_DialogOkButton))
                    go.setFixedSize(24, 24)
                    go.clicked.connect(s._goto_frame)
                    row1.addWidget(go)
                    layout.addLayout(row1)
                    row2 = QHBoxLayout()
                    row2.addWidget(QLabel("Time [s]"))
                    s.time_entry = QLineEdit("0.0")
                    s.time_entry.returnPressed.connect(s._goto_time)
                    row2.addWidget(s.time_entry)
                    go2 = QPushButton()
                    go2.setIcon(s.style().standardIcon(QStyle.SP_DialogOkButton))
                    go2.setFixedSize(24, 24)
                    go2.clicked.connect(s._goto_time)
                    row2.addWidget(go2)
                    layout.addLayout(row2)
                s._build_right(main_layout, _nav)

            def _goto_frame(s):
                try:
                    s.slider.setValue(max(0, min(n_times - 1, int(s.frame_entry.text()) - 1)))
                except ValueError:
                    pass

            def _goto_time(s):
                try:
                    t = float(s.time_entry.text())
                    s.slider.setValue(int(np.argmin(np.abs(time_arr - t))))
                except ValueError:
                    pass

            def _init_plot(s):
                s.figure.clear()
                s._axes = []
                ax = s.figure.add_subplot(1, 1, 1)
                s._axes.append(ax)
                s.figure.subplots_adjust(left=0.12, right=0.97, top=0.90, bottom=0.10)
                apply_dark_figure_style(s.figure)

            def _update_plot(s, t_idx):
                efit_data, real_t_idx, t_actual, sl_shot, sl_tree = s.slices[t_idx]
                eq = efit_data.get('eq2d', {})
                ms = int(round(t_actual * 1000))
                s.figure.suptitle(f"#{sl_shot} ({sl_tree})  {ms:06d}ms", fontsize=12)
                s.frame_entry.setText(str(t_idx + 1))
                s.time_entry.setText(f"{t_actual:.4f}")

                for a in s._artists:
                    a.remove()
                s._artists = []

                ax = s._axes[0]
                ax.clear()
                psirz = eq['psirz']
                rgrid = eq['rgrid']
                zgrid = eq['zgrid']
                ssimag = eq.get('ssimag')
                ssibry = eq.get('ssibry')
                ti = real_t_idx
                psi_2d = psirz[ti] if psirz.ndim == 3 else psirz
                if ssimag is not None and ssibry is not None:
                    pm = float(ssimag[ti]) if ssimag.ndim > 0 else float(ssimag)
                    pb = float(ssibry[ti]) if ssibry.ndim > 0 else float(ssibry)
                    d = pb - pm
                    psi_N = (psi_2d - pm) / d if abs(d) > 1e-12 else psi_2d
                else:
                    psi_N = psi_2d

                n = s._n_contours
                levels = np.linspace(0, 1, n + 1)[1:]
                ax.contour(rgrid, zgrid, psi_N, levels=levels,
                           colors=['#1f77b4'], linewidths=0.8, alpha=0.8)
                ax.contour(rgrid, zgrid, psi_N, levels=[1.0],
                           colors=['#1f77b4'], linewidths=2.0)
                # SOL
                sol_n = s._sol_contours
                if sol_n > 0:
                    sol = [1.0 + (i+1)*0.05 for i in range(sol_n)]
                    ax.contour(rgrid, zgrid, psi_N, levels=sol,
                               colors=['#1f77b4'], linewidths=0.5,
                               linestyles='dashed', alpha=0.5)
                # Boundary
                rb = eq.get('rbbbs')
                zb = eq.get('zbbbs')
                if rb is not None and zb is not None:
                    r = rb[ti] if rb.ndim == 2 else rb
                    z = zb[ti] if zb.ndim == 2 else zb
                    mask = (r != 0) | (z != 0)
                    if np.any(mask):
                        ax.plot(r[mask], z[mask], '-', color='#1f77b4', lw=2)
                # Axis + X-point
                rm = eq.get('rmaxis')
                zm = eq.get('zmaxis')
                if rm is not None and zm is not None:
                    ax.plot(float(rm[ti]) if rm.ndim > 0 else float(rm),
                            float(zm[ti]) if zm.ndim > 0 else float(zm),
                            '+', color='red', markersize=10, mew=2)
                for xr, xz in [('rxpt1','zxpt1'), ('rxpt2','zxpt2')]:
                    xrd = eq.get(xr); xzd = eq.get(xz)
                    if xrd is not None and xzd is not None:
                        xi = min(ti, len(xrd)-1)
                        rx, zx = float(xrd[xi]), float(xzd[xi])
                        if rx > 0 and abs(rx) < 5:
                            ax.plot(rx, zx, 'x', color='red', markersize=8, mew=2)
                # Limiter
                rl = eq.get('rlim'); zl = eq.get('zlim')
                if rl is not None and zl is not None:
                    ax.plot(np.append(rl, rl[0]), np.append(zl, zl[0]),
                            '-', color='gray', lw=1.5, alpha=0.5)
                ax.set_aspect('equal')
                ax.set_xlabel('R [m]', fontsize=12)
                ax.set_ylabel('Z [m]', fontsize=12)
                ax.tick_params(labelsize=10)
                ax.grid(ls='--', lw=0.3, color='#444444')
                apply_dark_figure_style(s.figure)
                s.canvas.draw_idle()

            def _add_current(s):
                _, _, t_actual, sl_shot, sl_tree = s.slices[s.slider.value()]
                ms = int(round(t_actual * 1000))
                entry = f"{sl_shot}_{ms:06d} ({sl_tree})"
                if s.selected_listbox is not None:
                    existing = {s.selected_listbox.item(i).text()
                                for i in range(s.selected_listbox.count())}
                    if entry in existing:
                        s.add_btn.setText(f"{ms:06d}ms already selected")
                        return
                    s.selected_listbox.addItem(entry)
                s.added_count += 1
                s.add_btn.setText(f"Added {ms:06d}ms  ({s.added_count} total)")

        dlg = _Efit2DBrowse(self.frame, merged_slices, shot, tree,
                            selected_listbox=self.selected_listbox)
        dlg.show()

    def _plot(self):
        entries = [self.selected_listbox.item(i).text()
                   for i in range(self.selected_listbox.count())]
        if not entries:
            return

        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)

        n_contours = getattr(self, '_n_contours', 9)
        checked = self._ct_group.checkedButton()
        contour_type = checked.property('ct_key') if checked else 'psi_N'
        sol_n = getattr(self, '_sol_contours', 4)
        sol_dist = getattr(self, '_sol_distance', 0.01)

        # Interior contour levels (0 to 1)
        interior_levels = np.linspace(0, 1, n_contours + 1)[1:]
        # SOL levels based on distance (psi_N increment ~ sol_distance per step)
        sol_levels = [1.0 + (i + 1) * sol_dist for i in range(sol_n)] if sol_n > 0 else []

        colors = self._get_plot_colors(len(entries))
        _limiters = {}

        for idx, entry in enumerate(entries):
            try:
                cache_key, time_s, shot, tree, ms = self._parse_entry(entry)
            except Exception as e:
                print(f"[EFIT] Error parsing {entry}: {e}")
                continue

            efit = self._efit_cache.get(cache_key)
            if efit is None:
                continue
            eq2d = efit.get('eq2d', {})
            if 'psirz' not in eq2d or 'rgrid' not in eq2d or 'zgrid' not in eq2d:
                continue

            time_arr = eq2d['time']
            t_idx = int(np.argmin(np.abs(time_arr - time_s)))
            color = colors[idx]

            psirz = eq2d['psirz']
            rgrid = eq2d['rgrid']
            zgrid = eq2d['zgrid']
            ssimag = eq2d.get('ssimag')
            ssibry = eq2d.get('ssibry')

            psi_2d = psirz[t_idx] if psirz.ndim == 3 else psirz
            if ssimag is not None and ssibry is not None:
                psi_mag = float(ssimag[t_idx]) if ssimag.ndim > 0 else float(ssimag)
                psi_bry = float(ssibry[t_idx]) if ssibry.ndim > 0 else float(ssibry)
                denom = psi_bry - psi_mag
                psi_N = (psi_2d - psi_mag) / denom if abs(denom) > 1e-12 else psi_2d
            else:
                psi_N = psi_2d

            # Convert contour grid based on type
            if contour_type == 'rho_pol':
                contour_2d = np.sqrt(np.clip(psi_N, 0, None))
                levels = np.sqrt(np.clip(interior_levels, 0, None))
                sep_level = 1.0
                sol = [np.sqrt(v) for v in sol_levels] if sol_levels else []
            elif contour_type == 'rho_tor':
                # Use RHOVN lookup: rho_tor = interp(psi_N -> rhovn)
                rhovn = eq2d.get('rhovn')
                if rhovn is not None:
                    from scipy.interpolate import interp1d
                    rhovn_1d = rhovn[t_idx] if rhovn.ndim == 2 else rhovn
                    psin_grid = np.linspace(0, 1, len(rhovn_1d))
                    rhovn_interp = interp1d(psin_grid, rhovn_1d,
                                            fill_value='extrapolate')
                    contour_2d = rhovn_interp(np.clip(psi_N, 0, None))
                    levels = rhovn_interp(interior_levels)
                    sep_level = float(rhovn_interp(1.0))
                    sol = [float(rhovn_interp(v)) for v in sol_levels] if sol_levels else []
                else:
                    contour_2d = psi_N
                    levels = interior_levels
                    sep_level = 1.0
                    sol = sol_levels
            else:  # psi_N (default)
                contour_2d = psi_N
                levels = interior_levels
                sep_level = 1.0
                sol = sol_levels

            # Interior contour lines
            label = f'#{shot} {ms:06d}ms ({tree})'
            ax.contour(rgrid, zgrid, contour_2d, levels=levels,
                       colors=[color], linewidths=0.8, alpha=0.8)

            # Separatrix (thicker)
            ax.contour(rgrid, zgrid, contour_2d, levels=[sep_level],
                       colors=[color], linewidths=2.0)

            # SOL contours (dashed)
            if sol:
                ax.contour(rgrid, zgrid, contour_2d, levels=sol,
                           colors=[color], linewidths=0.5, linestyles='dashed',
                           alpha=0.5)

            # Boundary
            rbbbs = eq2d.get('rbbbs')
            zbbbs = eq2d.get('zbbbs')
            if rbbbs is not None and zbbbs is not None:
                rb = rbbbs[t_idx] if rbbbs.ndim == 2 else rbbbs
                zb = zbbbs[t_idx] if zbbbs.ndim == 2 else zbbbs
                mask = (rb != 0.0) | (zb != 0.0)
                if np.any(mask):
                    ax.plot(rb[mask], zb[mask], '-', color=color, lw=2, label=label)

            # Magnetic axis (+)
            rmaxis = eq2d.get('rmaxis')
            zmaxis = eq2d.get('zmaxis')
            if rmaxis is not None and zmaxis is not None:
                rm = float(rmaxis[t_idx]) if rmaxis.ndim > 0 else float(rmaxis)
                zm = float(zmaxis[t_idx]) if zmaxis.ndim > 0 else float(zmaxis)
                ax.plot(rm, zm, '+', color=color, markersize=10, mew=2)

            # X-points (x marker)
            for xr_key, xz_key in [('rxpt1', 'zxpt1'), ('rxpt2', 'zxpt2')]:
                xr = eq2d.get(xr_key)
                xz = eq2d.get(xz_key)
                if xr is not None and xz is not None:
                    # X-point time arrays may differ in length from gtime
                    xi = min(t_idx, len(xr) - 1)
                    rx = float(xr[xi])
                    zx = float(xz[xi])
                    if rx > 0 and abs(rx) < 5:  # valid (not -9.99)
                        ax.plot(rx, zx, 'x', color=color, markersize=8, mew=2)

            # Collect limiter data per cache_key
            if cache_key not in _limiters:
                rlim = eq2d.get('rlim')
                zlim = eq2d.get('zlim')
                if rlim is not None and zlim is not None:
                    _limiters[cache_key] = (rlim, zlim, [color])
                else:
                    _limiters[cache_key] = (None, None, [color])
            else:
                _limiters[cache_key][2].append(color)

        # Draw limiters: if all same shape → gray, if different → colored
        lim_shapes = {}
        for key, (rl, zl, cols) in _limiters.items():
            if rl is not None:
                shape_key = (tuple(np.round(rl, 5)), tuple(np.round(zl, 5)))
                lim_shapes.setdefault(shape_key, []).extend(cols)

        if len(lim_shapes) == 1:
            rl, zl = list(lim_shapes.keys())[0]
            ax.plot(list(rl) + [rl[0]], list(zl) + [zl[0]],
                    '-', color='gray', lw=1.5, alpha=0.5)
        elif len(lim_shapes) > 1:
            for (rl, zl), cols in lim_shapes.items():
                for c in set(cols):
                    ax.plot(list(rl) + [rl[0]], list(zl) + [zl[0]],
                            '-', color=c, lw=1.2, alpha=0.4)

        ax.set_aspect('equal')
        ax.set_xlabel('R [m]', fontsize=self.label_fontsize)
        ax.set_ylabel('Z [m]', fontsize=self.label_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        if ax.get_legend_handles_labels()[1]:
            ax.legend(fontsize=self.legend_fontsize, loc='upper right', frameon=False)

        self.figure.subplots_adjust(
            left=0.12, right=0.97, top=0.95, bottom=0.10)
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()

    def _plot_multi(self, entries):
        """Plot 2D contour of last entry + boundary overlays from all entries."""
        last_entry = entries[-1]
        try:
            cache_key, time_s, shot, tree, ms = self._parse_entry(last_entry)
        except Exception as e:
            print(f"[EFIT] Error parsing entry {last_entry}: {e}")
            return

        efit = self._efit_cache.get(cache_key)
        if efit is None:
            return
        eq2d = efit.get('eq2d', {})

        ax = self.figure.add_subplot(1, 1, 1)

        # Plot 2D contour of the last entry as background if available
        if 'psirz' in eq2d and 'rgrid' in eq2d and 'zgrid' in eq2d:
            time_arr = eq2d['time']
            t_idx = int(np.argmin(np.abs(time_arr - time_s)))

            psirz = eq2d['psirz']
            rgrid = eq2d['rgrid']
            zgrid = eq2d['zgrid']
            ssimag = eq2d.get('ssimag')
            ssibry = eq2d.get('ssibry')

            if psirz.ndim == 3:
                psi_2d = psirz[t_idx]
            else:
                psi_2d = psirz

            if ssimag is not None and ssibry is not None:
                psi_mag = float(ssimag[t_idx]) if len(ssimag) > t_idx else float(ssimag)
                psi_bry = float(ssibry[t_idx]) if len(ssibry) > t_idx else float(ssibry)
                denom = psi_bry - psi_mag
                if abs(denom) > 1e-12:
                    psi_N = (psi_2d - psi_mag) / denom
                else:
                    psi_N = psi_2d
            else:
                psi_N = psi_2d

            cf = ax.contourf(rgrid, zgrid, psi_N, levels=50, cmap=self._cmap,
                             alpha=0.4)
            ax.contour(rgrid, zgrid, psi_N, levels=[1.0],
                       colors='gray', linewidths=0.5, alpha=0.5)

        # Overlay boundaries from all entries in different colors
        colors = self._get_plot_colors(len(entries))

        for idx, entry in enumerate(entries):
            try:
                ck, ts, sh, tr, m_s = self._parse_entry(entry)
                ef = self._efit_cache.get(ck)
                if ef is None:
                    continue
                eq = ef.get('eq2d', {})
                ta = eq.get('time')
                if ta is None:
                    continue
                ti = int(np.argmin(np.abs(ta - ts)))

                # Boundary
                rbbbs = eq.get('rbbbs')
                zbbbs = eq.get('zbbbs')
                if rbbbs is not None and zbbbs is not None:
                    if rbbbs.ndim == 2:
                        rb = rbbbs[ti]
                        zb = zbbbs[ti]
                    else:
                        rb = rbbbs
                        zb = zbbbs
                    mask = (rb != 0.0) | (zb != 0.0)
                    if np.any(mask):
                        ax.plot(rb[mask], zb[mask], '-', color=colors[idx], lw=2,
                                label=f'#{sh} {m_s}ms ({tr})')

                # Magnetic axis
                rmaxis = eq.get('rmaxis')
                zmaxis = eq.get('zmaxis')
                if rmaxis is not None and zmaxis is not None:
                    rm = float(rmaxis[ti]) if rmaxis.ndim > 0 and len(rmaxis) > ti else float(rmaxis)
                    zm = float(zmaxis[ti]) if zmaxis.ndim > 0 and len(zmaxis) > ti else float(zmaxis)
                    ax.plot(rm, zm, 'x', color=colors[idx], markersize=8, mew=2)

            except Exception as e:
                print(f"[EFIT] Error plotting {entry}: {e}")

        ax.set_aspect('equal')
        ax.set_xlabel('R [m]', fontsize=self.label_fontsize)
        ax.set_ylabel('Z [m]', fontsize=self.label_fontsize)
        ax.tick_params(labelsize=self.tick_fontsize)
        ax.set_title('EFIT 2D Boundary Comparison', fontsize=self.label_fontsize)

        if ax.get_legend_handles_labels()[1]:
            ax.legend(fontsize=self.legend_fontsize, loc='upper right', frameon=False)

        self.figure.subplots_adjust(
            left=0.12, right=0.97, top=0.95, bottom=0.10)

    # ================================================================
    # Style dialog
    # ================================================================

    def _show_style_dialog(self):
        W = 150
        dlg = QDialog(self.frame)
        dlg.setWindowTitle("Plot Options")
        dlg.setMinimumWidth(300)
        dl = QVBoxLayout(dlg)

        # Colormap selection
        cmap_row = QHBoxLayout()
        cmap_row.addWidget(QLabel("Colormap"))
        cmap_combo = QComboBox()
        cmap_combo.setFixedWidth(W)
        cmap_combo.addItems([
            "RdBu_r", "plasma", "inferno", "viridis", "jet",
            "coolwarm", "seismic", "hot",
        ])
        cmap_combo.setCurrentText(self._cmap)
        cmap_row.addWidget(cmap_combo)
        dl.addLayout(cmap_row)

        # Boundary color mode
        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Boundary Color"))
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

        # Number of contours
        nc_row = QHBoxLayout()
        nc_row.addWidget(QLabel("Contour levels"))
        nc_spin = QSpinBox()
        nc_spin.setFixedWidth(W)
        nc_spin.setRange(3, 50)
        nc_spin.setValue(getattr(self, '_n_contours', 9))
        nc_row.addWidget(nc_spin)
        dl.addLayout(nc_row)

        # SOL contours
        sol_row = QHBoxLayout()
        sol_row.addWidget(QLabel("SOL contours"))
        sol_spin = QSpinBox()
        sol_spin.setFixedWidth(W)
        sol_spin.setRange(0, 10)
        sol_spin.setValue(getattr(self, '_sol_contours', 4))
        sol_row.addWidget(sol_spin)
        dl.addLayout(sol_row)

        # SOL distance
        from PySide6.QtWidgets import QDoubleSpinBox
        sd_row = QHBoxLayout()
        sd_row.addWidget(QLabel("SOL distance [m]"))
        sd_spin = QDoubleSpinBox()
        sd_spin.setFixedWidth(W)
        sd_spin.setRange(0.001, 0.1)
        sd_spin.setDecimals(3)
        sd_spin.setSingleStep(0.005)
        sd_spin.setValue(getattr(self, '_sol_distance', 0.01))
        sd_row.addWidget(sd_spin)
        dl.addLayout(sd_row)

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
            self._cmap = cmap_combo.currentText()
            self._color_mode = color_combo.currentText()
            self._n_contours = nc_spin.value()
            self._sol_contours = sol_spin.value()
            self._sol_distance = sd_spin.value()
            for spin in dlg.findChildren(QSpinBox):
                if spin.property('attr'):
                    setattr(self, spin.property('attr'), spin.value())
            self._plot()

        def _reset():
            cmap_combo.setCurrentText("RdBu_r")
            color_combo.setCurrentText("Gradient(viridis)")
            nc_spin.setValue(9)
            sol_spin.setValue(4)
            sd_spin.setValue(0.01)
            for spin in dlg.findChildren(QSpinBox):
                if spin.property('default') is not None:
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
        """Build CSV lines for boundary R,Z export."""
        entries = [self.selected_listbox.item(i).text()
                   for i in range(self.selected_listbox.count())]
        if not entries:
            return None

        lines = []
        lines.append("# EFIT 2D Boundary Data\n")
        lines.append("#%10s,%10s,%10s,%14s,%14s\n" % (
            "Shot", "Tree", "Time[s]", "R_boundary[m]", "Z_boundary[m]"))

        has_data = False
        for entry in entries:
            try:
                cache_key, time_s, shot, tree, ms = self._parse_entry(entry)
                efit = self._efit_cache.get(cache_key)
                if efit is None:
                    continue
                eq2d = efit.get('eq2d', {})
                time_arr = eq2d.get('time')
                if time_arr is None:
                    continue
                t_idx = int(np.argmin(np.abs(time_arr - time_s)))
                t_actual = time_arr[t_idx]

                rbbbs = eq2d.get('rbbbs')
                zbbbs = eq2d.get('zbbbs')
                if rbbbs is None or zbbbs is None:
                    continue

                if rbbbs.ndim == 2:
                    rb = rbbbs[t_idx]
                    zb = zbbbs[t_idx]
                else:
                    rb = rbbbs
                    zb = zbbbs

                mask = (rb != 0.0) | (zb != 0.0)
                tree_str = tree if tree else 'gfile'
                for r, z in zip(rb[mask], zb[mask]):
                    lines.append(" %10s,%10s,%10.4f,%14.6e,%14.6e\n" % (
                        shot, tree_str, t_actual, r, z))
                    has_data = True
            except Exception as e:
                print(f"[EFIT] Error building CSV for {entry}: {e}")

        return lines if has_data else None

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
        self._cmap = s.get("cmap", "RdBu_r")
        self._color_mode = s.get("color_mode", "Gradient(viridis)")
        self.label_fontsize = s.get("label_fontsize", 12)
        self.legend_fontsize = s.get("legend_fontsize", 8)
        self.tick_fontsize = s.get("tick_fontsize", 10)
        saved_shot = s.get("shot", "")
        if saved_shot:
            self.shot_entry.setText(saved_shot)

    def save_settings(self):
        from config.user_settings import get_tab_settings, set_tab_settings
        s = get_tab_settings(self._settings_key)
        s["cmap"] = self._cmap
        s["color_mode"] = getattr(self, '_color_mode', 'Gradient(viridis)')
        s["label_fontsize"] = self.label_fontsize
        s["legend_fontsize"] = self.legend_fontsize
        s["tick_fontsize"] = self.tick_fontsize
        s["shot"] = self.shot_entry.text().strip()
        set_tab_settings(self._settings_key, s)
