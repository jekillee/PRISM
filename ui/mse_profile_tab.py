"""
MSE Profile tab with j/q selection
"""

import numpy as np
from scipy.interpolate import interp1d
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QMessageBox, QLineEdit, QComboBox, QPushButton, QCheckBox,
    QWidget, QHBoxLayout, QVBoxLayout, QGroupBox, QGridLayout, QLabel,
    QSplitter, QStyle, QScrollArea, QRadioButton, QFrame,
    QSpinBox, QDialog, QDialogButtonBox,
)
from PySide6.QtCore import Qt

from ui.profile_base_tab import ProfileBaseTab
from ui.ui_constants import CONTROL_PANEL_WIDTH, get_icon
from ui.theme import ThemeManager


class MSEProfileTab(ProfileBaseTab):
    """MSE Profile tab with TGAMMA and j/q profiles"""

    TAB_NAME = "MSE"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _get_secondary_loader(self):
        """MSE has no secondary diagnostic"""
        return None

    def create_widgets(self):
        """Create MSE profile tab widgets"""
        self.figure = Figure(self.app_config.FIGURE_SIZE)
        self.figure.subplots_adjust(left=0.10, right=0.97, top=0.93, bottom=0.10, wspace=0.20)

        # Setup 1x2 plot with shared x-axis: TGAMMA (left), j or q (right)
        self.ax1 = self.figure.add_subplot(121)
        self.ax2 = self.figure.add_subplot(122, sharex=self.ax1)

        self.ax1.set_xlabel('R [m]')
        self.ax1.set_ylabel(r'$\gamma$ [rad]')
        self.ax2.set_xlabel('R [m]')
        self.ax2.set_ylabel('q')

        # Create canvas
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.draw()

        # Left side: canvas + toolbar in a vertical layout
        canvas_widget = QWidget()
        canvas_layout = QVBoxLayout(canvas_widget)
        canvas_layout.setContentsMargins(0, 0, 0, 0)
        canvas_layout.addWidget(self.canvas)

        # Right side: scrollable control panel
        scroll_area = QScrollArea()
        scroll_area.setFixedWidth(CONTROL_PANEL_WIDTH)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)

        control_frame = QWidget()
        control_layout = QVBoxLayout(control_frame)
        control_layout.setContentsMargins(9, 9, 9, 9)
        control_layout.setSizeConstraint(QVBoxLayout.SetMinimumSize)

        scroll_area.setWidget(control_frame)
        scroll_area.viewport().setAutoFillBackground(False)
        control_frame.setAutoFillBackground(False)

        # Splitter for left/right layout
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(canvas_widget)
        splitter.addWidget(scroll_area)

        # Set the main layout of self.frame
        main_layout = QVBoxLayout(self.frame)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(splitter)

        self._create_shot_input(control_frame)
        self._create_selection_listboxes(control_frame)
        self._create_plot_controls(control_frame)
        self._create_efit_controls(control_frame)
        self._create_save_controls(control_frame, section_num=5)
        control_layout.addStretch()

    def _create_shot_input(self, parent):
        """Create data loading section"""
        group = QGroupBox("1. Load MSE Data")
        grid = QGridLayout(group)
        grid.setColumnStretch(1, 1)

        # Row 0: Shot label, entry with up/down, Fetch
        grid.addWidget(QLabel('Shot'), 0, 0)

        shot_frame = QWidget()
        shot_layout = QHBoxLayout(shot_frame)
        shot_layout.setContentsMargins(0, 0, 0, 0)

        self.shot_entry = QLineEdit()
        shot_layout.addWidget(self.shot_entry, 1)
        self.shot_entry.returnPressed.connect(self.load_shot_data)

        btn_updown = QWidget()
        btn_updown_layout = QVBoxLayout(btn_updown)
        btn_updown_layout.setContentsMargins(0, 0, 0, 0)
        btn_updown_layout.setSpacing(0)
        mini_btn_style = "padding: 0px; border-radius: 2px;"
        up_btn = QPushButton()
        up_btn.setIcon(get_icon(QStyle.SP_ArrowUp))
        up_btn.setFixedSize(24, 15)
        up_btn.setStyleSheet(mini_btn_style)
        up_btn.clicked.connect(lambda: self._adjust_shot(1))
        btn_updown_layout.addWidget(up_btn)
        down_btn = QPushButton()
        down_btn.setIcon(get_icon(QStyle.SP_ArrowDown))
        down_btn.setFixedSize(24, 15)
        down_btn.setStyleSheet(mini_btn_style)
        down_btn.clicked.connect(lambda: self._adjust_shot(-1))
        btn_updown_layout.addWidget(down_btn)
        shot_layout.addWidget(btn_updown)

        grid.addWidget(shot_frame, 0, 1)

        self.fetch_button = QPushButton('Fetch')
        self.fetch_button.setFixedWidth(70)
        self.fetch_button.clicked.connect(self.load_shot_data)
        grid.addWidget(self.fetch_button, 0, 2)

        parent.layout().addWidget(group)

    def _adjust_shot(self, delta):
        """Adjust shot number by delta"""
        try:
            current = int(self.shot_entry.text())
            new_shot = max(1, current + delta)
            self.shot_entry.setText(str(new_shot))
        except ValueError:
            pass

    def _get_selected_param(self):
        """Get selected parameter from radio buttons"""
        return 'q' if self.param_q_radio.isChecked() else 'j'

    def _create_plot_controls(self, parent):
        """Create plot control buttons"""
        group = QGroupBox("3. Plot")
        group_layout = QVBoxLayout(group)

        # Row 1: q/j radio buttons + Plot button
        row1 = QHBoxLayout()

        self.param_q_radio = QRadioButton(' q')
        self.param_q_radio.setChecked(True)
        row1.addWidget(self.param_q_radio)

        self.param_j_radio = QRadioButton(' j')
        row1.addWidget(self.param_j_radio)

        plot_button = QPushButton("Plot")
        plot_button.clicked.connect(self.plot_data)
        row1.addWidget(plot_button, 2)

        style_btn = QPushButton("Option")
        style_btn.clicked.connect(self._show_style_dialog)
        row1.addWidget(style_btn, 1)

        group_layout.addLayout(row1)

        # Separator line
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        group_layout.addWidget(separator)

        # Hidden combo for color mode (used by _get_plot_colors, saved to settings)
        self.color_mode_combo = QComboBox()
        self.color_mode_combo.addItems([
            "Gradient(viridis)", "Gradient(hot)", "Gradient(jet)", "Gradient(coolwarm)",
            "Fixed(tab10)", "Fixed(tab20)", "Fixed(Set1)", "Fixed(Set2)", "Fixed(Set3)",
        ])
        self.color_mode_combo.setCurrentText("Gradient(viridis)")
        self.color_mode_combo.hide()

        # Default font sizes
        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10

        # Show Nodes + Select Channels
        row2 = QHBoxLayout()

        self.show_channel_checkbox = QCheckBox("Show Nodes")
        self.show_channel_checkbox.setChecked(False)
        row2.addWidget(self.show_channel_checkbox)

        channels_btn = QPushButton("Select Channels")
        channels_btn.setToolTip("Select which channels to enable or dim")
        channels_btn.clicked.connect(self._show_channel_selector)
        row2.addWidget(channels_btn)

        group_layout.addLayout(row2)

        parent.layout().addWidget(group)

    def load_shot_data(self):
        """Load MSE shot data from MDS+"""
        try:
            shot_number = int(self.shot_entry.text())

            data = self.data_loader.load_data(shot_number)

            cache_key = f'{shot_number}_MSE'
            self.data[cache_key] = data

            self.available_listbox.clear()

            # Use profile time for time selection
            for tp in data.time_prof:
                item_str = f'{shot_number:06d}_{tp*1e3:06.0f} (MSE)'
                self.available_listbox.addItem(item_str)

            print(f"[MSE] Data loaded: {len(data.time_prof)} timepoints")
            print(f"[MSE]   NB source: {data.nb_source}")

        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load data: {str(e)}")

    def _parse_entry(self, entry):
        """Parse entry to get shot and time"""
        # entry format: '034142_000730 (MSE)'
        # Remove ' (MSE)' suffix first
        entry_clean = entry.replace(' (MSE)', '')
        parts = entry_clean.split('_')
        shot_number = int(parts[0])
        time_point = float(parts[1]) / 1e3  # ms to s
        return shot_number, time_point, 'MSE'

    def _parse_entry_for_efit(self, entry):
        """Parse entry for EFIT (extract shot and time)"""
        shot_number, time_point, _ = self._parse_entry(entry)
        return shot_number, time_point

    def _get_r_edge_at_time(self, data, time_point):
        """Get R_edge at specific time by interpolation"""
        r_edge_info = data.measurements['r_edge']
        interp_func = interp1d(r_edge_info['time'], r_edge_info['data'],
                               fill_value='extrapolate')
        return float(interp_func(time_point))

    def _get_channel_info(self):
        """Return TGAMMA channel info based on selected listbox entries."""
        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            return []

        # Find loaded MSE data matching any selected entry
        for entry in selected_entries:
            shot_number, _, _ = self._parse_entry(entry)
            cache_key = f'{shot_number}_MSE'
            if cache_key in self.data:
                data = self.data[cache_key]
                if hasattr(data, 'measurements') and 'tgamma' in data.measurements:
                    tgamma = data.measurements['tgamma']
                    good_mask = tgamma['good_mask']
                    info = []
                    for i in range(len(good_mask)):
                        if good_mask[i]:
                            info.append((f"TGAMMA_{i}", f"TGAMMA{i+1:02d}"))
                    return info
        return []

    def plot_data(self):
        """Plot R profiles"""
        self.ax1.clear()
        self.ax2.clear()

        self.ax1.set_xlabel('R [m]')
        self.ax1.set_ylabel(r'$\gamma$ [rad]')
        self.ax2.set_xlabel('R [m]')

        # Set y-label based on selected parameter
        param = self._get_selected_param()
        if param == 'q':
            self.ax2.set_ylabel('q')
        else:
            self.ax2.set_ylabel('j [MA/m$^2$]')

        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            return

        gamma_min, gamma_max = np.inf, -np.inf
        param_min, param_max = np.inf, -np.inf
        colors = self._get_plot_colors(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                shot_number, time_point, _ = self._parse_entry(entry)
                color = colors[i]

                cache_key = f'{shot_number}_MSE'

                if cache_key not in self.data:
                    data = self.data_loader.load_data(shot_number)
                    self.data[cache_key] = data
                else:
                    data = self.data[cache_key]

                # Get time indices
                idx_raw = np.argmin(np.abs(data.time - time_point))
                idx_prof = np.argmin(np.abs(data.time_prof - time_point))

                # Get R_edge and magx at this time
                r_edge = self._get_r_edge_at_time(data, time_point)
                magx = data.measurements['magx']['data'][idx_prof]

                # =============================================================
                # Plot TGAMMA vs R (left plot)
                # =============================================================
                tgamma_meas = data.measurements['tgamma']
                good_mask = tgamma_meas['good_mask']

                R_raw = tgamma_meas['R'][good_mask]
                drr = tgamma_meas['drr'][good_mask]
                tgamma = tgamma_meas['data'][good_mask, idx_raw]
                sgamma = tgamma_meas['error'][good_mask, idx_raw]

                label = f'#{shot_number} {int(time_point*1e3):06d}ms (MSE)'

                # Channel mask for user-disabled channels
                good_indices = [j for j in range(len(good_mask)) if good_mask[j]]
                ch_keys = [f"TGAMMA_{j}" for j in good_indices]
                ch_mask = self._get_channel_mask(ch_keys)

                if ch_mask.any():
                    self.ax1.errorbar(R_raw[ch_mask], tgamma[ch_mask],
                                     xerr=drr[ch_mask]*0.5, yerr=sgamma[ch_mask]*0.5,
                                     fmt='o-', capsize=3, markersize=5, color=color,
                                     linewidth=0.8, label=label)
                if (~ch_mask).any():
                    self.ax1.errorbar(R_raw[~ch_mask], tgamma[~ch_mask],
                                     xerr=drr[~ch_mask]*0.5, yerr=sgamma[~ch_mask]*0.5,
                                     fmt='o', capsize=3, markersize=5,
                                     color=(0.6, 0.6, 0.6, 0.35), linewidth=0.8, label='')

                # Add channel labels for TGAMMA (node: TGAMMA01~25)
                raw_channels = [j+1 for j in range(len(good_mask)) if good_mask[j]]
                self._add_channel_labels(self.ax1, R_raw, tgamma, 'TGAMMA', raw_channels)

                gamma_min = min(gamma_min, np.nanpercentile(tgamma, 2))
                gamma_max = max(gamma_max, np.nanpercentile(tgamma, 98))

                # =============================================================
                # Plot j or q vs R (right plot)
                # - Solid line: pmse profile (20 points)
                # - Markers: interpolated values at TGAMMA R positions
                # =============================================================
                param_meas = data.measurements[param]
                roa = param_meas['roa'][:, idx_prof]
                param_data = param_meas['data'][:, idx_prof]
                param_err = param_meas['error'][:, idx_prof]

                # Convert r/a to R
                R_prof = self.data_loader.get_R_from_roa(roa, magx, r_edge)

                # Sort by R for proper line plotting
                sort_idx = np.argsort(R_prof)
                R_prof_sorted = R_prof[sort_idx]
                param_data_sorted = param_data[sort_idx]
                param_err_sorted = param_err[sort_idx]

                # Plot pmse profile as solid line with error band
                self.ax2.plot(R_prof_sorted, param_data_sorted, '-', color=color,
                             linewidth=1.5, label=label)
                self.ax2.fill_between(R_prof_sorted,
                                      param_data_sorted - param_err_sorted,
                                      param_data_sorted + param_err_sorted,
                                      color=color, alpha=0.2)

                # Interpolate param values at TGAMMA R positions
                # Use R_raw from TGAMMA (good channels only)
                param_interp_func = interp1d(R_prof_sorted, param_data_sorted,
                                             fill_value='extrapolate', bounds_error=False)
                param_at_tgamma = param_interp_func(R_raw)

                # Plot markers at TGAMMA positions (apply channel mask)
                if ch_mask.any():
                    self.ax2.plot(R_raw[ch_mask], param_at_tgamma[ch_mask], 'o', color=color,
                                 markersize=5, markeredgecolor='white', markeredgewidth=0.5)
                if (~ch_mask).any():
                    self.ax2.plot(R_raw[~ch_mask], param_at_tgamma[~ch_mask], 'o',
                                 color=(0.6, 0.6, 0.6, 0.35),
                                 markersize=5, markeredgecolor='white', markeredgewidth=0.5)

                # Add channel labels at TGAMMA positions (node: TGAMMA01~25)
                self._add_channel_labels(self.ax2, R_raw, param_at_tgamma, 'TGAMMA', raw_channels)

                param_max = max(param_max, np.nanmax(param_data))
                param_min = min(param_min, np.nanmin(param_data))

            except Exception as e:
                print(f"[MSE] Error plotting {entry}: {str(e)}")

        # Set limits
        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'
        self.ax1.axhline(0, color=zc, linestyle='--', gid='zero_ref')
        gamma_margin = (gamma_max - gamma_min) * 0.1
        self.ax1.set_ylim(gamma_min - gamma_margin, gamma_max + gamma_margin)

        if param == 'q':
            self.ax2.axhline(1, color=zc, linestyle='--', gid='zero_ref')
        else:
            self.ax2.axhline(0, color=zc, linestyle='--', gid='zero_ref')
        self.ax2.set_ylim(0, 6)

        self.plot_manager.apply_common_styling(
            self.ax1, self.ax2,
            legend_fontsize=self.legend_fontsize,
            label_fontsize=self.label_fontsize,
            tick_fontsize=self.tick_fontsize)
        self.canvas.draw()

        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()

    def plot_efit_profiles(self):
        """Plot profiles with EFIT mapping"""
        if not self.efit_data or self.computed_efit_tree is None:
            QMessageBox.warning(self.frame, "Warning", "Please compute EFIT first.")
            return

        efit_tree = self.computed_efit_tree

        self.ax1.clear()
        self.ax2.clear()

        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            return

        x_axis = self._get_selected_x_axis()

        if x_axis == "psi_N":
            x_label = rf"$\psi_N$ ({efit_tree})"
        elif x_axis == "rho_pol":
            x_label = rf"$\rho_{{pol}}$ ({efit_tree})"
        else:
            x_label = rf"$\rho_{{tor}}$ ({efit_tree})"

        param = self._get_selected_param()

        gamma_max, gamma_min, param_max, param_min = 0, 0, 0, 0
        colors = self._get_plot_colors(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                shot_number, time_point, _ = self._parse_entry(entry)

                if entry not in self.efit_data:
                    continue

                color = colors[i]
                efit_entry = self.efit_data[entry]
                efit_R = efit_entry['R']

                if x_axis == "psi_N":
                    interp_func = interp1d(efit_R, efit_entry['psi_N'], fill_value='extrapolate')
                elif x_axis == "rho_pol":
                    interp_func = interp1d(efit_R, efit_entry['rho_pol'], fill_value='extrapolate')
                else:
                    interp_func = interp1d(efit_R, efit_entry['rho_tor'], fill_value='extrapolate')

                cache_key = f'{shot_number}_MSE'
                data = self.data[cache_key]

                idx_raw = np.argmin(np.abs(data.time - time_point))
                idx_prof = np.argmin(np.abs(data.time_prof - time_point))

                r_edge = self._get_r_edge_at_time(data, time_point)
                magx = data.measurements['magx']['data'][idx_prof]

                # TGAMMA
                tgamma_meas = data.measurements['tgamma']
                good_mask = tgamma_meas['good_mask']

                R_raw = tgamma_meas['R'][good_mask]
                drr = tgamma_meas['drr'][good_mask]
                tgamma = tgamma_meas['data'][good_mask, idx_raw]
                sgamma = tgamma_meas['error'][good_mask, idx_raw]

                x_raw = interp_func(R_raw)
                label = f'#{shot_number} {int(time_point*1e3):06d}ms (MSE)'

                # Channel mask for user-disabled channels
                good_indices = [j for j in range(len(good_mask)) if good_mask[j]]
                ch_keys = [f"TGAMMA_{j}" for j in good_indices]
                ch_mask = self._get_channel_mask(ch_keys)

                if ch_mask.any():
                    self.ax1.errorbar(x_raw[ch_mask], tgamma[ch_mask], yerr=sgamma[ch_mask]*0.5,
                                     fmt='o-', capsize=3, markersize=5, color=color,
                                     linewidth=0.8, label=label)
                if (~ch_mask).any():
                    self.ax1.errorbar(x_raw[~ch_mask], tgamma[~ch_mask], yerr=sgamma[~ch_mask]*0.5,
                                     fmt='o', capsize=3, markersize=5,
                                     color=(0.6, 0.6, 0.6, 0.35), linewidth=0.8, label='')

                # Add channel labels for TGAMMA
                raw_channels = [j+1 for j in range(len(good_mask)) if good_mask[j]]
                self._add_channel_labels(self.ax1, x_raw, tgamma, 'TGAMMA', raw_channels)

                gamma_min = min(gamma_min, np.nanpercentile(tgamma, 2))
                gamma_max = max(gamma_max, np.nanpercentile(tgamma, 98))

                # j or q
                # - Solid line: pmse profile (20 points)
                # - Markers: interpolated values at TGAMMA positions
                param_meas = data.measurements[param]
                roa = param_meas['roa'][:, idx_prof]
                param_data = param_meas['data'][:, idx_prof]
                param_err = param_meas['error'][:, idx_prof]

                R_prof = self.data_loader.get_R_from_roa(roa, magx, r_edge)
                x_prof = interp_func(R_prof)

                # Sort by x coordinate for proper line plotting
                sort_idx = np.argsort(x_prof)
                x_prof_sorted = x_prof[sort_idx]
                param_data_sorted = param_data[sort_idx]
                param_err_sorted = param_err[sort_idx]

                # Plot pmse profile as solid line with error band
                self.ax2.plot(x_prof_sorted, param_data_sorted, '-', color=color,
                             linewidth=1.5, label=label)
                self.ax2.fill_between(x_prof_sorted,
                                      param_data_sorted - param_err_sorted,
                                      param_data_sorted + param_err_sorted,
                                      color=color, alpha=0.2)

                # Interpolate param values at TGAMMA positions (in EFIT coordinates)
                param_interp_func = interp1d(x_prof_sorted, param_data_sorted,
                                             fill_value='extrapolate', bounds_error=False)
                param_at_tgamma = param_interp_func(x_raw)

                # Plot markers at TGAMMA positions (apply channel mask)
                if ch_mask.any():
                    self.ax2.plot(x_raw[ch_mask], param_at_tgamma[ch_mask], 'o', color=color,
                                 markersize=5, markeredgecolor='white', markeredgewidth=0.5)
                if (~ch_mask).any():
                    self.ax2.plot(x_raw[~ch_mask], param_at_tgamma[~ch_mask], 'o',
                                 color=(0.6, 0.6, 0.6, 0.35),
                                 markersize=5, markeredgecolor='white', markeredgewidth=0.5)

                # Add channel labels at TGAMMA positions
                self._add_channel_labels(self.ax2, x_raw, param_at_tgamma, 'TGAMMA', raw_channels)

                param_max = max(param_max, np.nanmax(param_data))
                param_min = min(param_min, np.nanmin(param_data))

            except Exception as e:
                print(f"[MSE] Error plotting {entry}: {str(e)}")

        self.ax1.set_xlabel(x_label)
        self.ax2.set_xlabel(x_label)
        self.ax1.set_ylabel(r'$\gamma$ [rad]')

        if param == 'q':
            self.ax2.set_ylabel('q')
        else:
            self.ax2.set_ylabel('j [MA/m$^2$]')

        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'
        self.ax1.axhline(0, color=zc, linestyle='--', gid='zero_ref')
        gamma_margin = (gamma_max - gamma_min) * 0.1
        self.ax1.set_ylim(gamma_min - gamma_margin, gamma_max + gamma_margin)

        if param == 'q':
            self.ax2.axhline(1, color=zc, linestyle='--', gid='zero_ref')
        else:
            self.ax2.axhline(0, color=zc, linestyle='--', gid='zero_ref')
        self.ax2.set_ylim(0, 6)

        for ax in [self.ax1, self.ax2]:
            ax.axvline(x=0, c='k', ls='--')
            ax.axvline(x=1, c='k', ls='--')

        self.ax1.set_xlim(-0.5, 1.05)
        self.ax2.set_xlim(0, 1.05)

        self.plot_manager.apply_common_styling(
            self.ax1, self.ax2,
            legend_fontsize=self.legend_fontsize,
            label_fontsize=self.label_fontsize,
            tick_fontsize=self.tick_fontsize)
        self.canvas.draw()

        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()

    def _write_data_to_file(self, file_path, selected_entries):
        """Write MSE profile data to text file (gamma, q, j at TGAMMA positions)"""
        with open(file_path, 'w') as f:
            # Write header
            f.write("# MSE Profile Data (q and j interpolated at TGAMMA positions)\n")
            if self.computed_efit_tree:
                f.write(f"# EFIT Source: {self.computed_efit_tree}\n")
            f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "psi_N", "rho_pol", "rho_tor",
                "gamma", "gamma_err", "drr[m]", "q", "q_err", "j[MA/m2]", "j_err"
            ))

            for entry in selected_entries:
                shot_number, time_point, _ = self._parse_entry(entry)
                cache_key = f'{shot_number}_MSE'

                if cache_key not in self.data:
                    data = self.data_loader.load_data(shot_number)
                    self.data[cache_key] = data
                else:
                    data = self.data[cache_key]

                idx_raw = np.argmin(np.abs(data.time - time_point))
                idx_prof = np.argmin(np.abs(data.time_prof - time_point))
                actual_time = data.time[idx_raw]

                # TGAMMA data
                tgamma_meas = data.measurements['tgamma']
                good_mask = tgamma_meas['good_mask']

                R_raw = tgamma_meas['R'][good_mask]
                drr = tgamma_meas['drr'][good_mask]
                tgamma = tgamma_meas['data'][good_mask, idx_raw]
                sgamma = tgamma_meas['error'][good_mask, idx_raw]

                # Get q and j profiles for interpolation
                r_edge = self._get_r_edge_at_time(data, time_point)
                magx = data.measurements['magx']['data'][idx_prof]

                q_meas = data.measurements['q']
                j_meas = data.measurements['j']
                roa = q_meas['roa'][:, idx_prof]
                R_prof = self.data_loader.get_R_from_roa(roa, magx, r_edge)

                # Sort by R for interpolation
                sort_idx = np.argsort(R_prof)
                R_prof_sorted = R_prof[sort_idx]

                q_data_sorted = q_meas['data'][:, idx_prof][sort_idx]
                q_err_sorted = q_meas['error'][:, idx_prof][sort_idx]
                j_data_sorted = j_meas['data'][:, idx_prof][sort_idx]
                j_err_sorted = j_meas['error'][:, idx_prof][sort_idx]

                # Interpolate q and j at TGAMMA R positions
                q_at_tgamma = np.interp(R_raw, R_prof_sorted, q_data_sorted)
                q_err_at_tgamma = np.interp(R_raw, R_prof_sorted, q_err_sorted)
                j_at_tgamma = np.interp(R_raw, R_prof_sorted, j_data_sorted)
                j_err_at_tgamma = np.interp(R_raw, R_prof_sorted, j_err_sorted)

                # Get EFIT values
                psi_n, rho_pol, rho_tor = self._get_efit_values_at_R(entry, R_raw)

                # Write data rows (gamma, q, j at same R positions)
                for i in range(len(R_raw)):
                    f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f\n" % (
                        shot_number, actual_time, R_raw[i],
                        self._format_value(psi_n[i]),
                        self._format_value(rho_pol[i]),
                        self._format_value(rho_tor[i]),
                        tgamma[i], sgamma[i], drr[i],
                        q_at_tgamma[i], q_err_at_tgamma[i],
                        j_at_tgamma[i], j_err_at_tgamma[i]
                    ))
