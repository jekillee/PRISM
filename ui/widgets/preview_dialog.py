"""
Profile Browse Dialog - slider-based data browsing before Select Data.
Canvas+slider left, controls right. Non-modal.
"""

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QSlider, QPushButton,
    QLineEdit, QComboBox, QStyle, QGroupBox, QWidget, QFrame,
)
from PySide6.QtCore import Qt, QTimer

from ui.ui_constants import apply_dark_figure_style, get_icon
from ui.theme import ThemeManager

_ZERO_BOTTOM = {'Ti', 'Te', 'ne'}


class ProfileBrowseDialog(QDialog):
    """Non-modal dialog for browsing profile data with slider and playback."""

    def __init__(self, parent, data, shot_number, source,
                 param1_name, param2_name, param1_label, param2_label,
                 selected_listbox=None, extra_data=None):
        super().__init__(parent)
        self.data = data
        self.extra_data = extra_data
        self.shot_number = shot_number
        self.source = source
        self.param1_name = param1_name
        self.param2_name = param2_name
        self.param1_label = param1_label
        self.param2_label = param2_label
        self.selected_listbox = selected_listbox
        self.added_count = 0

        self._artists = []
        self._axes = []
        self._fix_axes = False
        self._is_ece_only = (source == 'ECE')
        self.time_arr = getattr(data, 'time_prof', data.time)

        self._play_timer = QTimer()
        self._play_timer.timeout.connect(self._play_step)
        self._playing = False

        self.setWindowTitle(f"Browse  #{shot_number} ({source})")
        self.resize(1000, 550)
        self._build_ui()
        self._init_plot()
        self._update_plot(0)

    def _build_ui(self):
        from ui.widgets.toggle_switch import ToggleSwitch

        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(4)

        # === Left: Canvas + slider below ===
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)

        self.figure = Figure((8, 4.5))
        apply_dark_figure_style(self.figure)
        self.canvas = FigureCanvasQTAgg(self.figure)
        toolbar = NavigationToolbar2QT(self.canvas, self)
        left_layout.addWidget(toolbar)
        left_layout.addWidget(self.canvas, stretch=1)

        # Slider row below canvas
        slider_row = QHBoxLayout()
        prev_btn = QPushButton()
        prev_btn.setIcon(get_icon(QStyle.SP_ArrowBack))
        prev_btn.setFixedWidth(28)
        prev_btn.clicked.connect(lambda: self._step(-1))
        slider_row.addWidget(prev_btn)

        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(len(self.time_arr) - 1)
        self.slider.valueChanged.connect(self._on_slider)
        slider_row.addWidget(self.slider, stretch=1)

        next_btn = QPushButton()
        next_btn.setIcon(get_icon(QStyle.SP_ArrowForward))
        next_btn.setFixedWidth(28)
        next_btn.clicked.connect(lambda: self._step(1))
        slider_row.addWidget(next_btn)

        left_layout.addLayout(slider_row)

        hint = QLabel("Use mouse wheel or arrow keys to navigate frames.")
        hint.setStyleSheet("color: gray; font-size: 8pt;")
        hint.setAlignment(Qt.AlignCenter)
        left_layout.addWidget(hint)

        main_layout.addWidget(left, stretch=1)

        # === Right: Controls ===
        ctrl = QWidget()
        ctrl.setFixedWidth(240)
        cl = QVBoxLayout(ctrl)
        cl.setContentsMargins(4, 0, 4, 4)

        # --- Navigation ---
        nav_group = QGroupBox("Navigation")
        nav_layout = QVBoxLayout(nav_group)

        frame_row = QHBoxLayout()
        frame_row.addWidget(QLabel("Frame"))
        self.frame_entry = QLineEdit("0")
        self.frame_entry.returnPressed.connect(self._goto_frame)
        frame_row.addWidget(self.frame_entry)
        frame_row.addWidget(QLabel(f"/ {len(self.time_arr)-1}"))
        go_frame = QPushButton()
        go_frame.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
        go_frame.setFixedSize(24, 24)
        go_frame.clicked.connect(self._goto_frame)
        frame_row.addWidget(go_frame)
        nav_layout.addLayout(frame_row)

        time_row = QHBoxLayout()
        time_row.addWidget(QLabel("Time [s]"))
        self.time_entry = QLineEdit("0.0")
        self.time_entry.returnPressed.connect(self._goto_time)
        time_row.addWidget(self.time_entry)
        go_time = QPushButton()
        go_time.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
        go_time.setFixedSize(24, 24)
        go_time.clicked.connect(self._goto_time)
        time_row.addWidget(go_time)
        nav_layout.addLayout(time_row)

        cl.addWidget(nav_group)

        # --- Playback ---
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

        # --- Options ---
        opt_group = QGroupBox("Options")
        opt_layout = QHBoxLayout(opt_group)
        self.fix_axes_toggle = ToggleSwitch()
        self.fix_axes_toggle.toggled.connect(self._on_fix_axes)
        opt_layout.addWidget(self.fix_axes_toggle)
        opt_layout.addWidget(QLabel("Fix Axes"))
        opt_layout.addStretch()
        cl.addWidget(opt_group)

        cl.addStretch()

        # --- Actions ---
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
        for col, (pname, plabel) in enumerate([
            (self.param1_name, self.param1_label),
            (self.param2_name, self.param2_label),
        ]):
            ax = self.figure.add_subplot(1, 2, col + 1)
            ax.set_xlabel('R [m]', fontsize=12)
            ax.set_ylabel(plabel, fontsize=12)
            ax.tick_params(labelsize=10)
            ax.grid(ls='--', lw=0.3, color='#444444')
            if 'v' in pname.lower():
                ax.axhline(y=0, color=zc, ls='--', gid='zero_ref')
            self._axes.append(ax)
        self.figure.subplots_adjust(
            left=0.12, right=0.97, top=0.90, bottom=0.14, wspace=0.30)
        apply_dark_figure_style(self.figure)

    def _get_param_data(self, pname, t_actual):
        data = self.data
        if pname not in data.measurements:
            return None, None, None, None
        pdata, eu, el = data.get_parameter_asymmetric(pname)
        if hasattr(data, 'time_prof') and pname in ('q', 'j'):
            t_idx = np.argmin(np.abs(data.time_prof - t_actual))
            meas = data.measurements[pname]
            R = meas['roa'][:, t_idx] if 'roa' in meas else data.radius[:pdata.shape[0]]
        else:
            t_idx = np.argmin(np.abs(data.time - t_actual))
            R = data.radius[:pdata.shape[0]]
        return R, pdata[:, t_idx], eu[:, t_idx], el[:, t_idx]

    def _update_plot(self, t_idx):
        t_actual = self.time_arr[t_idx]
        ms = int(round(t_actual * 1000))

        self.figure.suptitle(
            f"#{self.shot_number}  {ms:06d}ms  ({self.source})", fontsize=12)
        self.frame_entry.setText(str(t_idx))
        self.time_entry.setText(f"{t_actual:.4f}")

        for a in self._artists:
            a.remove()
        self._artists = []

        if self._is_ece_only:
            marker, mfc, mew = 's', 'none', 1.5
        else:
            marker, mfc, mew = 'o', None, 1

        for col, pname in enumerate([self.param1_name, self.param2_name]):
            ax = self._axes[col]
            R, profile, eu, el = self._get_param_data(pname, t_actual)
            if R is None:
                continue
            valid = np.isfinite(profile)
            if not np.any(valid):
                continue

            kwargs = dict(fmt=marker, markersize=5, capsize=5,
                          color='#1f77b4', markeredgewidth=mew, zorder=10)
            if mfc is not None:
                kwargs['markerfacecolor'] = mfc
            c = ax.errorbar(R[valid], profile[valid],
                            yerr=[el[valid], eu[valid]], **kwargs)
            self._artists.append(c)

            if self.extra_data is not None and pname == 'Te':
                self._overlay_extra(ax, pname, t_actual)

            if not self._fix_axes:
                y_vals = profile[valid]
                y_min = np.nanmin(y_vals - el[valid])
                y_max = np.nanmax(y_vals + eu[valid])
                margin = (y_max - y_min) * 0.1 if y_max > y_min else 0.1
                if pname in _ZERO_BOTTOM:
                    ax.set_ylim(0, y_max + margin)
                else:
                    ax.set_ylim(y_min - margin, y_max + margin)

        self.canvas.draw_idle()

    def _overlay_extra(self, ax, pname, t_actual):
        ed = self.extra_data
        if pname not in ed.measurements:
            return
        t_idx = np.argmin(np.abs(ed.time - t_actual))
        te_data = ed.measurements[pname]['data'][:, t_idx]
        R_ece = ed.radius
        valid = np.isfinite(te_data) & (te_data > 0)
        if np.any(valid):
            line, = ax.plot(R_ece[valid], te_data[valid], 's',
                            color='#1f77b4', markersize=5,
                            markerfacecolor='none', markeredgewidth=1.5,
                            zorder=5)
            self._artists.append(line)

    # ---- Navigation ----

    def _on_slider(self, value):
        self._update_plot(value)

    def _step(self, delta):
        self.slider.setValue(
            max(0, min(self.slider.maximum(), self.slider.value() + delta)))

    def _goto_frame(self):
        try:
            self.slider.setValue(
                max(0, min(len(self.time_arr) - 1, int(self.frame_entry.text()))))
        except ValueError:
            pass

    def _goto_time(self):
        try:
            self.slider.setValue(
                int(np.argmin(np.abs(self.time_arr - float(self.time_entry.text())))))
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

    # ---- Add to Selected ----

    def _add_current(self):
        t_val = self.time_arr[self.slider.value()]
        ms = int(round(t_val * 1000))
        entry = f'{self.shot_number:06d}_{t_val*1e3:06.0f} ({self.source})'

        if self.selected_listbox is not None:
            existing = {self.selected_listbox.item(i).text()
                        for i in range(self.selected_listbox.count())}
            if entry in existing:
                self.add_btn.setText(f"{ms:06d}ms already selected")
                print(f"[Browse] {entry} already in Selected, skipping")
                return
            self.selected_listbox.addItem(entry)

        self.added_count += 1
        self.add_btn.setText(f"Added {ms:06d}ms  ({self.added_count} total)")
        print(f"[Browse] Added to Selected: {entry}")

    def _on_close(self):
        self._stop_play()
        self.accept()
