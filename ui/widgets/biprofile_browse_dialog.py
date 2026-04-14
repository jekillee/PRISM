"""
BiProfile Browse Dialog - slider-based browsing for BiProfile data.
Profile mode: browse time points (fitted profile vs psi_N at each time).
TimeTrace mode: browse psi_N positions (time evolution at each psi_N).
"""

import numpy as np
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QSlider, QPushButton,
    QLineEdit, QComboBox, QStyle, QGroupBox, QWidget,
)
from PySide6.QtCore import Qt, QTimer

from ui.ui_constants import apply_dark_figure_style, get_icon
from ui.theme import ThemeManager

UNITS = {'Ti': 'keV', 'vT': 'km/s', 'Te': 'keV', 'ne': r'$10^{19}$/m$^3$'}
LABELS = {'Ti': r'T$_i$', 'vT': r'v$_T$', 'Te': r'T$_e$', 'ne': r'n$_e$'}
_ZERO_BOTTOM = {'Ti', 'Te', 'ne'}


class BiProfileBrowseDialog(QDialog):
    """Non-modal dialog for browsing BiProfile data.

    mode='profile': slider browses time, plots profile vs psi_N
    mode='timetrace': slider browses psi_N, plots time trace
    """

    def __init__(self, parent, bi_data, shot_number, params, mode='profile',
                 selected_listbox=None, sdata=None, plot_color='#1f77b4'):
        """
        Args:
            bi_data: dict {param: {psin, time, mean, unc, ...}}
            shot_number: int
            params: tuple ('Ti','vT') or ('ne','Te')
            mode: 'profile' or 'timetrace'
            selected_listbox: QListWidget to add entries to
            sdata: full shot data dict (for raw overlay in profile mode)
            plot_color: color string
        """
        super().__init__(parent)
        self.bi_data = bi_data
        self.shot_number = shot_number
        self.params = params
        self.mode = mode
        self.selected_listbox = selected_listbox
        self.sdata = sdata
        self.plot_color = plot_color
        self.added_count = 0

        self._artists = []
        self._axes = []
        self._fix_axes = False

        ref = next((p for p in params if p in bi_data), None)
        d = bi_data[ref]
        if mode == 'profile':
            self.browse_arr = d['time']
            self.browse_label = "Time"
            self.browse_unit = "s"
        else:
            self.browse_arr = d['psin']
            self.browse_label = "\u03c8\u2099"
            self.browse_unit = ""

        self._play_timer = QTimer()
        self._play_timer.timeout.connect(self._play_step)
        self._playing = False

        mode_str = "Profile" if mode == 'profile' else "TimeTrace"
        self.setWindowTitle(
            f"Browse {mode_str}  #{shot_number} ({', '.join(params)})")
        self.resize(1000, 550)
        self._build_ui()
        self._init_plot()
        self._update_plot(0)

    def _build_ui(self):
        from ui.widgets.toggle_switch import ToggleSwitch

        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(4)

        # Left: Canvas + slider
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)

        self.figure = Figure((8, 4.5))
        apply_dark_figure_style(self.figure)
        self.canvas = FigureCanvasQTAgg(self.figure)
        toolbar = NavigationToolbar2QT(self.canvas, self)
        left_layout.addWidget(toolbar)
        left_layout.addWidget(self.canvas, stretch=1)

        slider_row = QHBoxLayout()
        prev_btn = QPushButton()
        prev_btn.setIcon(get_icon(QStyle.SP_ArrowBack))
        prev_btn.setFixedWidth(28)
        prev_btn.clicked.connect(lambda: self._step(-1))
        slider_row.addWidget(prev_btn)

        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(len(self.browse_arr) - 1)
        self.slider.valueChanged.connect(self._on_slider)
        slider_row.addWidget(self.slider, stretch=1)

        next_btn = QPushButton()
        next_btn.setIcon(get_icon(QStyle.SP_ArrowForward))
        next_btn.setFixedWidth(28)
        next_btn.clicked.connect(lambda: self._step(1))
        slider_row.addWidget(next_btn)
        left_layout.addLayout(slider_row)

        hint = QLabel("Use mouse wheel or arrow keys to navigate.")
        hint.setStyleSheet("color: gray; font-size: 8pt;")
        hint.setAlignment(Qt.AlignCenter)
        left_layout.addWidget(hint)
        main_layout.addWidget(left, stretch=1)

        # Right: Controls
        ctrl = QWidget()
        ctrl.setFixedWidth(240)
        cl = QVBoxLayout(ctrl)
        cl.setContentsMargins(4, 0, 4, 4)

        # Navigation
        nav_group = QGroupBox("Navigation")
        nav_layout = QVBoxLayout(nav_group)

        idx_row = QHBoxLayout()
        idx_row.addWidget(QLabel("Index"))
        self.idx_entry = QLineEdit("0")
        self.idx_entry.returnPressed.connect(self._goto_idx)
        idx_row.addWidget(self.idx_entry)
        idx_row.addWidget(QLabel(f"/ {len(self.browse_arr)-1}"))
        go_idx = QPushButton()
        go_idx.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
        go_idx.setFixedSize(24, 24)
        go_idx.clicked.connect(self._goto_idx)
        idx_row.addWidget(go_idx)
        nav_layout.addLayout(idx_row)

        val_row = QHBoxLayout()
        val_label = "Time [s]" if self.mode == 'profile' else "\u03c8\u2099"
        val_row.addWidget(QLabel(val_label))
        self.val_entry = QLineEdit("0.0")
        self.val_entry.returnPressed.connect(self._goto_val)
        val_row.addWidget(self.val_entry)
        go_val = QPushButton()
        go_val.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
        go_val.setFixedSize(24, 24)
        go_val.clicked.connect(self._goto_val)
        val_row.addWidget(go_val)
        nav_layout.addLayout(val_row)
        cl.addWidget(nav_group)

        # Playback
        play_group = QGroupBox("Playback")
        play_layout = QHBoxLayout(play_group)
        self.speed_combo = QComboBox()
        self.speed_combo.addItems(['1 fps', '2 fps', '5 fps', '10 fps', '20 fps', '50 fps'])
        self.speed_combo.setCurrentText('5 fps')
        play_layout.addWidget(self.speed_combo)
        self.play_btn = QPushButton()
        self.play_btn.setIcon(get_icon(QStyle.SP_MediaPlay))
        self.play_btn.clicked.connect(self._toggle_play)
        play_layout.addWidget(self.play_btn)
        stop_btn = QPushButton()
        stop_btn.setIcon(get_icon(QStyle.SP_MediaStop))
        stop_btn.clicked.connect(self._stop_play)
        play_layout.addWidget(stop_btn)
        cl.addWidget(play_group)

        # Options
        opt_group = QGroupBox("Options")
        opt_layout = QVBoxLayout(opt_group)

        fix_row = QHBoxLayout()
        self.fix_axes_toggle = ToggleSwitch()
        self.fix_axes_toggle.toggled.connect(self._on_fix_axes)
        fix_row.addWidget(self.fix_axes_toggle)
        fix_row.addWidget(QLabel("Fix Axes"))
        fix_row.addStretch()
        opt_layout.addLayout(fix_row)

        # Deapply TS Scale toggle (profile mode, ne/Te only) — scale applied by default
        self._apply_scale = True
        if self.mode == 'profile' and ('ne' in self.params or 'Te' in self.params):
            scale_row = QHBoxLayout()
            self.scale_toggle = ToggleSwitch()
            self.scale_toggle.toggled.connect(self._on_scale_toggled)
            scale_row.addWidget(self.scale_toggle)
            scale_row.addWidget(QLabel("Unapply TS ne Scale"))
            scale_row.addStretch()
            opt_layout.addLayout(scale_row)

        cl.addWidget(opt_group)

        cl.addStretch()

        self.add_btn = QPushButton("Add to Selected")
        self.add_btn.clicked.connect(self._add_current)
        cl.addWidget(self.add_btn)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self._on_close)
        cl.addWidget(close_btn)

        main_layout.addWidget(ctrl)

    # ---- Plot ----

    def _init_plot(self):
        self.figure.clear()
        self._axes = []
        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'

        for col, param in enumerate(self.params):
            ax = self.figure.add_subplot(1, 2, col + 1)
            if self.mode == 'profile':
                ax.set_xlabel(r'$\psi_N$', fontsize=12)
                ax.set_xlim(0, 1.1)
                ax.axvline(x=1.0, color=zc, ls='--', lw=0.8, alpha=0.5)
            else:
                ax.set_xlabel('Time [s]', fontsize=12)
            ax.set_ylabel(f'{LABELS.get(param, param)} [{UNITS.get(param, "")}]',
                          fontsize=12)
            ax.tick_params(labelsize=10)
            ax.grid(ls='--', lw=0.3, color='#444444')
            if param == 'vT':
                ax.axhline(y=0, color=zc, ls='--', gid='zero_ref')
            self._axes.append(ax)

        self.figure.subplots_adjust(
            left=0.12, right=0.97, top=0.90, bottom=0.14, wspace=0.30)
        apply_dark_figure_style(self.figure)

    def _update_plot(self, idx):
        val = self.browse_arr[idx]
        bi = self.bi_data
        color = self.plot_color

        if self.mode == 'profile':
            ms = int(round(val * 1000))
            title = f"#{self.shot_number}  {ms:06d}ms"
        else:
            title = f"#{self.shot_number}  \u03c8\u2099={val:.4f}"

        self.figure.suptitle(title, fontsize=12)
        self.idx_entry.setText(str(idx))
        self.val_entry.setText(f"{val:.4f}")

        for a in self._artists:
            a.remove()
        self._artists = []

        for col, param in enumerate(self.params):
            ax = self._axes[col]
            if param not in bi:
                continue
            d = bi[param]

            if self.mode == 'profile':
                # Profile at time idx
                psin = d['psin']
                t_idx = np.argmin(np.abs(d['time'] - val))
                mean_p = d['mean'][:, t_idx]
                unc_p = d['unc'][:, t_idx]
                valid = ~np.isnan(mean_p)
                if np.any(valid):
                    line, = ax.plot(psin[valid], mean_p[valid], '-',
                                    color=color, lw=2)
                    self._artists.append(line)
                    fill = ax.fill_between(psin[valid],
                                           mean_p[valid] - unc_p[valid],
                                           mean_p[valid] + unc_p[valid],
                                           alpha=0.2, color=color)
                    self._artists.append(fill)

                    if not self._fix_axes:
                        y_min = np.nanmin(mean_p[valid] - unc_p[valid])
                        y_max = np.nanmax(mean_p[valid] + unc_p[valid])
                        margin = (y_max - y_min) * 0.1 if y_max > y_min else 0.1
                        if param in _ZERO_BOTTOM:
                            ax.set_ylim(0, y_max + margin)
                        else:
                            ax.set_ylim(y_min - margin, y_max + margin)

                # Raw data overlay
                if self.sdata:
                    self._overlay_raw(ax, param, val, color)
            else:
                # TimeTrace at psi_N idx
                time = d['time']
                mean_t = d['mean'][idx, :]
                unc_t = d['unc'][idx, :]
                valid = ~np.isnan(mean_t)
                if np.any(valid):
                    line, = ax.plot(time[valid], mean_t[valid], '-',
                                    color=color, lw=2)
                    self._artists.append(line)
                    fill = ax.fill_between(time[valid],
                                           mean_t[valid] - unc_t[valid],
                                           mean_t[valid] + unc_t[valid],
                                           alpha=0.2, color=color)
                    self._artists.append(fill)

                    if not self._fix_axes:
                        y_min = np.nanmin(mean_t[valid] - unc_t[valid])
                        y_max = np.nanmax(mean_t[valid] + unc_t[valid])
                        margin = (y_max - y_min) * 0.1 if y_max > y_min else 0.1
                        if param in _ZERO_BOTTOM:
                            ax.set_ylim(0, y_max + margin)
                        else:
                            ax.set_ylim(y_min - margin, y_max + margin)
                    ax.set_xlim(time[0], time[-1])

        self.canvas.draw_idle()

    # ---- Navigation ----

    def _on_slider(self, value):
        self._update_plot(value)

    def _step(self, delta):
        self.slider.setValue(
            max(0, min(self.slider.maximum(), self.slider.value() + delta)))

    def _goto_idx(self):
        try:
            self.slider.setValue(
                max(0, min(len(self.browse_arr) - 1, int(self.idx_entry.text()))))
        except ValueError:
            pass

    def _goto_val(self):
        try:
            self.slider.setValue(
                int(np.argmin(np.abs(self.browse_arr - float(self.val_entry.text())))))
        except ValueError:
            pass

    def _on_fix_axes(self, checked):
        self._fix_axes = checked
        if not checked:
            self._update_plot(self.slider.value())

    # ---- Playback ----

    def _get_fps(self):
        return int(self.speed_combo.currentText().split()[0])

    def _toggle_play(self):
        if self._playing:
            self._stop_play()
        else:
            self._playing = True
            self.play_btn.setIcon(get_icon(QStyle.SP_MediaPause))
            self._play_timer.start(max(1, 1000 // self._get_fps()))

    def _stop_play(self):
        self._playing = False
        self.play_btn.setIcon(get_icon(QStyle.SP_MediaPlay))
        self._play_timer.stop()

    def _play_step(self):
        if self.slider.value() >= self.slider.maximum():
            self._stop_play()
            return
        self.slider.setValue(self.slider.value() + 1)
        self._play_timer.setInterval(max(1, 1000 // self._get_fps()))

    # ---- Raw overlay (profile mode only) ----

    def _on_scale_toggled(self, checked):
        self._apply_scale = not checked
        self._update_plot(self.slider.value())

    def _overlay_raw(self, ax, param, t_actual, color):
        from data_loaders.biprofile_loader import map_R_to_psin
        sdata = self.sdata
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
        bi = self.bi_data
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

    def _plot_raw(self, ax, psin, y, yerr, use_flags, color, marker):
        yerr = np.where(yerr < 0, 0, yerr)
        valid = (psin >= 0) & (psin <= 1.2) & np.isfinite(y)
        used = valid & (use_flags == 1)
        excluded = valid & (use_flags == 0)
        unknown = valid & np.isnan(use_flags)

        if np.any(used):
            c = ax.errorbar(psin[used], y[used], yerr=yerr[used],
                            fmt=marker, markersize=5, color=color,
                            capsize=5, zorder=10, markeredgewidth=1)
            self._artists.append(c)
        if np.any(excluded):
            c = ax.errorbar(psin[excluded], y[excluded], yerr=yerr[excluded],
                            fmt=marker, markersize=5, color=(0.6, 0.6, 0.6, 0.35),
                            capsize=5, zorder=5, markeredgewidth=1)
            self._artists.append(c)
        if np.any(unknown):
            c = ax.errorbar(psin[unknown], y[unknown], yerr=yerr[unknown],
                            fmt=marker, markersize=5, color=color,
                            capsize=5, zorder=10, markeredgewidth=1)
            self._artists.append(c)

    # ---- Add to Selected ----

    def _add_current(self):
        idx = self.slider.value()
        val = self.browse_arr[idx]

        if self.mode == 'profile':
            ms = int(round(val * 1000))
            entry = f"{self.shot_number:06d}_{val*1e3:06.0f} (Bi)"
            display = f"{ms:06d}ms"
        else:
            entry = f"{self.shot_number:06d}_{idx:03d} (Bi)"
            display = f"idx={idx} (\u03c8\u2099={val:.4f})"

        if self.selected_listbox is not None:
            existing = {self.selected_listbox.item(i).text()
                        for i in range(self.selected_listbox.count())}
            if entry in existing:
                self.add_btn.setText(f"{display} already selected")
                print(f"[Browse] {entry} already in Selected, skipping")
                return
            self.selected_listbox.addItem(entry)

        self.added_count += 1
        self.add_btn.setText(f"Added {display}  ({self.added_count} total)")
        print(f"[Browse] Added to Selected: {entry}")

    def _on_close(self):
        self._stop_play()
        self.accept()
