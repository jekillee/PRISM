"""
BiProfile Derived Profile tab.

For each selected time slice, plot one or more derived plasma quantities
(p_e, dT/dR, R/L_T, beta, nu*, omega_pe, c_s, rho_i, E_r, ...) vs ψ_N.

Reuses the Available/Selected time-picker pattern from biprofile_profile_tab,
adds a quantity-selector list and options for Zeff / Er mode / ion species.
"""

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QScrollArea, QLayout,
    QGroupBox, QGridLayout, QLabel, QLineEdit, QPushButton, QFrame,
    QListWidget, QListWidgetItem, QAbstractItemView, QFileDialog,
    QApplication, QMessageBox, QStyle, QSpinBox, QDoubleSpinBox,
    QComboBox, QRadioButton, QButtonGroup, QCheckBox,
    QDialog, QDialogButtonBox,
    QTableWidget, QTableWidgetItem, QHeaderView,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QShortcut, QKeySequence

from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, apply_dark_figure_style,
    apply_listbox_arrow_icons, apply_shot_arrow_icons,
)
from ui.theme import ThemeManager
from core.derived_quantities import (
    DerivedQuantities, PROFILE_QUANTITIES, get_quantity_meta,
)


class BiProfileDerivedTab:
    """BiProfile Derived Profile tab — multi-quantity, multi-time."""

    def __init__(self, main_window, params=None):
        # params unused (derived uses Ti, vT, Te, ne all together)
        self.main = main_window
        self._current_shot = None
        self._current_data = None
        self._derived = None         # cached DerivedQuantities instance

        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10
        self._color_mode = 'Gradient(viridis)'

        # Options (defaults)
        self.zeff = 2.0
        self.ion_species = 'D'
        self.imp_species = 'C6'
        self.er_mode = 'neo'   # 'neo' or 'vp_zero'

        self._settings_key = "bi_derived_profile"

        self.frame = QWidget()
        self.canvas = None
        self._create_widgets()
        self._restore_settings()

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------

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
        self._create_time_selector(control_frame)
        self._create_quantity_selector(control_frame)
        self._create_plot_controls(control_frame)
        self._create_save_controls(control_frame)
        control_layout.addStretch()

    # ---- 1. Load Data ----

    def _create_shot_input(self, parent):
        group = QGroupBox("1. Load Shot (Ti, vT, Te, ne)")
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

    # ---- 2. Select Quantity ----

    def _create_quantity_selector(self, parent):
        from PySide6.QtWidgets import QSizePolicy
        group = QGroupBox("3. Select Variable")
        layout = QVBoxLayout(group)

        self.quantity_combo = QComboBox()
        # Constrain to the parent control panel's column width: combo fills
        # the available horizontal space but does NOT grow to fit long items.
        # Long text is elided in the closed combo; full text appears in the
        # popup list (popup width follows items, not the closed combo).
        self.quantity_combo.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.quantity_combo.setSizeAdjustPolicy(
            QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.quantity_combo.setMinimumContentsLength(8)
        self.quantity_combo.setStyleSheet("QComboBox { combobox-popup: 0; }")
        self.quantity_combo.view().setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        for key, name, label, unit, _req in PROFILE_QUANTITIES:
            self.quantity_combo.addItem(f"{key}  —  {name}  [{unit}]", userData=key)
        self.quantity_combo.currentIndexChanged.connect(self._on_quantity_changed)
        layout.addWidget(self.quantity_combo)
        parent.layout().addWidget(group)

    def _on_quantity_changed(self):
        # No conditional options anymore (Zeff is hard-coded to 2.0, vθ-mode
        # is selected via separate quantity entries E_r_vp0 / E_r_neo).
        return

    # ---- 3. Times ----

    def _create_time_selector(self, parent):
        group = QGroupBox("2. Select Times")
        layout = QVBoxLayout(group)

        self.browse_button = QPushButton("Fetch a shot to browse")
        self.browse_button.setEnabled(False)
        self.browse_button.clicked.connect(self._open_browse)
        layout.addWidget(self.browse_button)

        lists_row = QHBoxLayout()

        avail_col = QVBoxLayout()
        avail_col.addWidget(QLabel("Time [ms]"))
        self.available_listbox = QListWidget()
        self.available_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.available_listbox.setFixedHeight(140)
        avail_col.addWidget(self.available_listbox)
        lists_row.addLayout(avail_col, stretch=1)

        btn_col = QVBoxLayout()
        btn_col.addStretch()
        add_btn = QPushButton()
        add_btn.setFixedWidth(30); add_btn.clicked.connect(self._add_times)
        btn_col.addWidget(add_btn)
        rm_btn = QPushButton()
        rm_btn.setFixedWidth(30); rm_btn.clicked.connect(self._remove_times)
        btn_col.addWidget(rm_btn)
        btn_col.addStretch()
        apply_listbox_arrow_icons(add_btn, rm_btn)
        lists_row.addLayout(btn_col)

        sel_col = QVBoxLayout()
        sel_col.addWidget(QLabel("Selected"))
        self.selected_listbox = QListWidget()
        self.selected_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.selected_listbox.setFixedHeight(140)
        sel_col.addWidget(self.selected_listbox)
        lists_row.addLayout(sel_col, stretch=1)

        layout.addLayout(lists_row)
        QShortcut(QKeySequence(Qt.Key_Delete), self.selected_listbox
                  ).activated.connect(self._remove_times)
        parent.layout().addWidget(group)

    # ---- 4. Options (Zeff / Er) ----

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

        sep = QFrame(); sep.setFrameShape(QFrame.HLine); sep.setFrameShadow(QFrame.Sunken)
        layout.addWidget(sep)

        self.info_label = QLabel("")
        self.info_label.setWordWrap(True)
        self.info_label.setStyleSheet("color: #888; font-size: 11px;")
        layout.addWidget(self.info_label)
        parent.layout().addWidget(group)

    def _create_save_controls(self, parent):
        group = QGroupBox("5. Save Data")
        layout = QVBoxLayout(group)
        btn = QPushButton("Preview && Save")
        btn.clicked.connect(self._preview_and_save)
        layout.addWidget(btn)
        parent.layout().addWidget(group)

    # ------------------------------------------------------------------
    # Settings persistence
    # ------------------------------------------------------------------

    def _restore_settings(self):
        from config.user_settings import get_tab_settings
        s = get_tab_settings(self._settings_key)
        if s.get("shot"):
            self.shot_entry.setText(str(s["shot"]))
        self._color_mode = s.get("color_mode", "Gradient(viridis)")
        self.label_fontsize = s.get("label_fontsize", 12)
        self.legend_fontsize = s.get("legend_fontsize", 8)
        self.tick_fontsize = s.get("tick_fontsize", 10)
        # Restore previously selected quantity if any
        if 'quantity' in s and hasattr(self, 'quantity_combo'):
            idx = self.quantity_combo.findData(s['quantity'])
            if idx >= 0:
                self.quantity_combo.setCurrentIndex(idx)
        # Apply conditional enable based on the (possibly restored) quantity
        self._on_quantity_changed()

    def save_settings(self):
        from config.user_settings import get_tab_settings, set_tab_settings
        s = get_tab_settings(self._settings_key)
        s["shot"] = self.shot_entry.text()
        s["color_mode"] = self._color_mode
        s["label_fontsize"] = self.label_fontsize
        s["legend_fontsize"] = self.legend_fontsize
        s["tick_fontsize"] = self.tick_fontsize
        if hasattr(self, 'quantity_combo'):
            s["quantity"] = self.quantity_combo.currentData()
        set_tab_settings(self._settings_key, s)

    # ------------------------------------------------------------------
    # Data fetch
    # ------------------------------------------------------------------

    def _adjust_shot(self, delta):
        try:
            self.shot_entry.setText(str(max(1, int(self.shot_entry.text()) + delta)))
        except ValueError:
            pass

    def _update_status(self, message, color='blue'):
        from ui.ui_constants import show_status
        show_status(self.frame, "BiProfile Derived", message, color)

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
            sdata = self.main.fetch_biprofile_shot(shot, params=('Ti', 'vT', 'Te', 'ne'))
        finally:
            self.fetch_button.setText("Fetch"); self.fetch_button.setEnabled(True)
            QApplication.restoreOverrideCursor()

        if not sdata or not sdata.get('bi'):
            self._update_status(f"#{shot}: No data", color='red')
            return

        bi = sdata['bi']
        # Need all four params for full derived suite; warn but allow partial
        missing = [p for p in ('Ti', 'vT', 'Te', 'ne') if p not in bi]
        loaded = [p for p in ('Ti', 'vT', 'Te', 'ne') if p in bi]
        if not loaded:
            self._update_status(f"#{shot}: No usable data", color='red')
            return

        self._current_shot = shot
        self._current_data = sdata
        self._derived = None    # invalidate cache; recompute on plot

        # Populate available times from first available param's time grid
        ref = loaded[0]
        self.available_listbox.clear()
        for t in bi[ref]['time']:
            self.available_listbox.addItem(f"{shot:06d}_{t*1e3:06.0f}")

        # EFIT is loaded by `fetch_biprofile_shot` (which calls load_efit_psin).
        # v2.5.6's load_efit_psin always returns the extended geometry keys
        # (q_psin_grid, q_t_psi, fpol_t_psi, bcentr, rmaxis), so no manual
        # re-load is required here. If EFIT load failed entirely (sdata['efit']
        # is None), derived quantities that need EFIT will be NaN; user is
        # notified via status below.
        efit_used = bi[ref].get('efit_used', '?')
        if sdata.get('efit') is None:
            self._update_status(
                f"#{shot}: EFIT load failed — derived quantities that need q/B "
                f"will be NaN.", color='orange')

        msg = f"#{shot}: loaded {', '.join(loaded)}"
        if missing:
            msg += f" (missing {', '.join(missing)})"
        self._update_status(msg, color='green' if not missing else 'orange')
        self.info_label.setText(
            f"EFIT: {efit_used}\n"
            f"{len(bi[ref]['time'])} time pts, {len(bi[ref]['psin'])} ψₙ pts")

        # Enable browse button
        self.browse_button.setEnabled(True)
        self.browse_button.setText(f"Browse #{shot}")

    def _add_times(self):
        existing = {self.selected_listbox.item(i).text()
                    for i in range(self.selected_listbox.count())}
        for item in self.available_listbox.selectedItems():
            if item.text() not in existing:
                self.selected_listbox.addItem(item.text())

    def _remove_times(self):
        for item in reversed(self.selected_listbox.selectedItems()):
            self.selected_listbox.takeItem(self.selected_listbox.row(item))

    def _open_browse(self):
        """Open the Derived browse dialog — slider through times, dropdown to
        switch quantity. Adds selected times to this tab's selected_listbox."""
        if not self._current_data:
            return
        dq = self._ensure_derived()
        if dq is None:
            return
        from ui.widgets.preview_dialog import BiProfileDerivedPreviewDialog
        dlg = BiProfileDerivedPreviewDialog(
            self.frame, dq, self._current_shot,
            selected_listbox=self.selected_listbox)
        dlg.show()

    # ------------------------------------------------------------------
    # Compute & cache
    # ------------------------------------------------------------------

    def _ensure_derived(self):
        """Create or update the DerivedQuantities computer with current options.

        EFIT is expected to have been loaded by `_fetch` already (eager load).
        We don't re-check / re-load here — `DerivedQuantities` has defensive
        fallbacks that return NaN for quantities depending on missing EFIT
        fields, so partial data still works.
        """
        if not self._current_data:
            QMessageBox.warning(self.frame, "Warning", "Fetch a shot first.")
            return None
        efit = self._current_data.get('efit')
        if efit is None:
            QMessageBox.warning(
                self.frame, "Warning",
                "EFIT data not available — derived quantities require EFIT.\n"
                "Re-fetch the shot.")
            return None
        # E_r mode is now selected via separate quantities (E_r_vp0 / E_r_neo).
        # Z_eff is hard-coded to 2.0; D / C6 are fixed for KSTAR carbon-wall.
        opts = (self.zeff, 'D', 'C6')
        if (self._derived is not None and
            getattr(self._derived, '_opts_signature', None) == opts):
            return self._derived
        try:
            self._derived = DerivedQuantities(
                self._current_data['bi'], efit,
                zeff=opts[0], ion_species=opts[1], imp_species=opts[2])
        except Exception as e:
            import traceback; traceback.print_exc()
            QMessageBox.critical(self.frame, "Derived setup error",
                                 f"Failed to build derived computer:\n{e}")
            return None
        self._derived._opts_signature = opts
        self.zeff = opts[0]
        self.er_mode = opts[3]
        return self._derived

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------

    def _parse_entry(self, entry):
        main = entry.split('(')[0].strip()
        parts = main.split('_')
        return int(parts[0]), float(parts[1]) / 1e3

    def _selected_quantities(self):
        key = self.quantity_combo.currentData() if hasattr(self, 'quantity_combo') else None
        return [key] if key else []

    def _get_plot_colors(self, n):
        import matplotlib.pyplot as plt
        mode = self._color_mode
        start = mode.find('('); end = mode.find(')')
        cmap_name = mode[start+1:end] if start != -1 and end != -1 else 'viridis'
        cmap = plt.get_cmap(cmap_name)
        if mode.startswith('Fixed'):
            return [cmap.colors[i % len(cmap.colors)] for i in range(n)]
        return [cmap(0.5)] if n == 1 else [cmap(i / (n - 1)) for i in range(n)]

    def _plot(self):
        qkeys = self._selected_quantities()
        entries = [self.selected_listbox.item(i).text()
                   for i in range(self.selected_listbox.count())]
        if not qkeys:
            QMessageBox.information(self.frame, "Info", "Select at least one quantity.")
            return
        if not entries:
            QMessageBox.information(self.frame, "Info", "Add at least one time to Selected.")
            return
        dq = self._ensure_derived()
        if dq is None:
            return

        try:
            self._do_plot(dq, qkeys, entries)
        except Exception as e:
            import traceback; traceback.print_exc()
            QMessageBox.critical(self.frame, "Plot error",
                                 f"Plot failed:\n{e}\n(See console for traceback.)")

    def _do_plot(self, dq, qkeys, entries):
        # Lay out subplots: ceil(sqrt) grid
        n = len(qkeys)
        ncols = 2 if n >= 2 else 1
        nrows = (n + ncols - 1) // ncols

        self.figure.clear()
        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'
        axes = []
        for idx, key in enumerate(qkeys):
            ax = self.figure.add_subplot(nrows, ncols, idx + 1)
            meta = get_quantity_meta(key)
            if meta:
                name, label, unit, _ = meta
            else:
                name, label, unit = key, key, ''
            ax.set_xlabel(r'$\psi_\mathrm{N}$', fontsize=self.label_fontsize)
            ax.set_ylabel(f'{name}, {label}  [{unit}]', fontsize=self.label_fontsize)
            ax.tick_params(labelsize=self.tick_fontsize)
            ax.set_xlim(0, 1.1)
            ax.grid(ls='--', lw=0.3, color='#444444')
            ax.axvline(x=1.0, color=zc, ls='--', lw=0.8, alpha=0.5)
            if key.startswith('nu_star'):
                ax.set_yscale('log')
            axes.append((ax, key))

        n_lines = 0
        colors = self._get_plot_colors(len(entries))
        for t_idx, entry in enumerate(entries):
            s, t_target = self._parse_entry(entry)
            color = colors[t_idx]
            ms = int(round(t_target * 1000))
            for ax, key in axes:
                d = dq.compute(key)
                if d is None:
                    continue
                ti = int(np.argmin(np.abs(d['time'] - t_target)))
                psin = d['psin']; vals = d['value'][:, ti]
                valid = np.isfinite(vals)
                if not np.any(valid):
                    continue
                ax.plot(psin[valid], vals[valid], '-', color=color, lw=2,
                        label=f'#{s} {ms:06d}ms')
                n_lines += 1

        for ax, _key in axes:
            if ax.get_legend_handles_labels()[1]:
                ax.legend(fontsize=self.legend_fontsize, loc='best', frameon=False)

        self.figure.subplots_adjust(
            left=0.10, right=0.97, top=0.95, bottom=0.10,
            wspace=0.27, hspace=0.32)
        apply_dark_figure_style(self.figure)
        self.canvas.draw_idle()

        if n_lines == 0:
            self._update_status(
                "No valid data — all selected quantities returned NaN. "
                "Check console for warnings.", color='red')
        else:
            self._update_status(f"Plotted {n_lines} curve(s).", color='green')

    # ------------------------------------------------------------------
    # Save data
    # ------------------------------------------------------------------

    def _preview_and_save(self):
        lines = self._build_csv_lines()
        if not lines:
            QMessageBox.warning(self.frame, "Warning", "No data to preview")
            return
        self._show_data_preview(lines)

    def _build_csv_lines(self):
        qkeys = self._selected_quantities()
        entries = [self.selected_listbox.item(i).text()
                   for i in range(self.selected_listbox.count())]
        if not qkeys or not entries:
            return None
        dq = self._ensure_derived()
        if dq is None:
            return None

        lines = ["# BiProfile Derived Quantities\n"]
        lines.append(f"# Zeff={self.zeff}, ion={self.ion_species}, "
                     f"imp={self.imp_species}, Er v_theta={self.er_mode}\n")
        header_parts = ["Shot", "Time[s]", "psi_N"]
        for k in qkeys:
            meta = get_quantity_meta(k)
            unit = meta[2] if meta else ''
            header_parts.append(f"{k}[{unit}]")
        lines.append("#" + ','.join(f"{h:>14s}" for h in header_parts) + "\n")

        for entry in entries:
            s, t_target = self._parse_entry(entry)
            # Compute all quantities once per time slice
            psin_ref = None
            col_values = {}
            t_actual = None
            for k in qkeys:
                d = dq.compute(k)
                if d is None:
                    col_values[k] = None
                    continue
                ti = int(np.argmin(np.abs(d['time'] - t_target)))
                t_actual = d['time'][ti]
                psin_ref = d['psin']
                col_values[k] = d['value'][:, ti]
            if psin_ref is None:
                continue
            for i, psi in enumerate(psin_ref):
                row = f" {s:>14d},{t_actual:>14.4f},{psi:>14.6f}"
                for k in qkeys:
                    v = col_values[k]
                    if v is None or not np.isfinite(v[i]):
                        row += f",{'':>14s}"
                    else:
                        row += f",{v[i]:>14.6e}"
                lines.append(row + "\n")
        return lines

    def _show_data_preview(self, lines):
        from PySide6.QtGui import QFont
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
                    title_line = content if not title_line else f"{title_line} | {content}"
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
        dlg.setWindowTitle(f"Data Preview — Derived")
        dlg.resize(min(120 * n_cols, 1200), min(30 * n_rows + 100, 700))
        dl = QVBoxLayout(dlg)

        info = QLabel(f"{n_rows} rows × {n_cols} columns — {title_line}")
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

        btn_row = QHBoxLayout(); btn_row.addStretch()
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

        close_btn = QPushButton("Close"); close_btn.clicked.connect(dlg.close)
        btn_row.addWidget(close_btn)
        dl.addLayout(btn_row)
        dlg.show()

    # ------------------------------------------------------------------
    # Style dialog
    # ------------------------------------------------------------------

    def _show_style_dialog(self):
        W = 150
        dlg = QDialog(self.frame)
        dlg.setWindowTitle("Plot Options")
        dlg.setMinimumWidth(300)
        dl = QVBoxLayout(dlg)

        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Color"))
        color_combo = QComboBox(); color_combo.setFixedWidth(W)
        color_combo.addItems([
            "Gradient(viridis)", "Gradient(hot)", "Gradient(jet)", "Gradient(coolwarm)",
            "Fixed(tab10)", "Fixed(tab20)", "Fixed(Set1)", "Fixed(Set2)", "Fixed(Set3)",
        ])
        color_combo.setCurrentText(self._color_mode)
        color_row.addWidget(color_combo); dl.addLayout(color_row)

        for name, attr, lo, hi, default in [
            ("Label font size", 'label_fontsize', 6, 24, 12),
            ("Legend font size", 'legend_fontsize', 4, 20, 8),
            ("Tick font size", 'tick_fontsize', 6, 20, 10),
        ]:
            row = QHBoxLayout(); row.addWidget(QLabel(name))
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
        dlg.show()
