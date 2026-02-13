"""
MSE Time Trace tab with j/q selection
"""

import numpy as np
from scipy.interpolate import interp1d

from PySide6.QtWidgets import (
    QMessageBox, QLineEdit, QComboBox, QPushButton, QWidget,
    QHBoxLayout, QVBoxLayout, QGroupBox, QGridLayout, QLabel,
    QFileDialog, QStyle, QRadioButton,
)

from ui.timetrace_base_tab import TimeTraceBaseTab
from ui.ui_constants import get_icon


class MSETimeTraceTab(TimeTraceBaseTab):
    """MSE Time Trace tab with TGAMMA and j/q time traces"""

    TAB_NAME = "MSE"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _create_shot_input(self, parent):
        """Create data loading section"""
        group = QGroupBox("1. Load MSE Data")
        grid = QGridLayout(group)

        grid.setColumnStretch(1, 1)

        # Row 0: Shot label, entry with up/down, Fetch
        grid.addWidget(QLabel('Shot'), 0, 0)

        shot_frame = QWidget()
        shot_frame_layout = QHBoxLayout(shot_frame)
        shot_frame_layout.setContentsMargins(0, 0, 0, 0)

        self.shot_entry = QLineEdit()
        shot_frame_layout.addWidget(self.shot_entry, 1)
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
        shot_frame_layout.addWidget(btn_updown)

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
            self.shot_entry.clear()
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

        row1 = QHBoxLayout()

        self.param_q_radio = QRadioButton(' q')
        self.param_q_radio.setChecked(True)
        row1.addWidget(self.param_q_radio)

        self.param_j_radio = QRadioButton(' j')
        row1.addWidget(self.param_j_radio)

        plot_button = QPushButton('Plot')
        plot_button.clicked.connect(self.plot_data)
        row1.addWidget(plot_button, 2)

        group_layout.addLayout(row1)

        parent.layout().addWidget(group)

    def _get_r_edge_at_time(self, data, time_point):
        """Get R_edge at specific time by interpolation"""
        r_edge_info = data.measurements['r_edge']
        interp_func = interp1d(r_edge_info['time'], r_edge_info['data'],
                               fill_value='extrapolate')
        return float(interp_func(time_point))

    def load_shot_data(self):
        """Load MSE shot data from MDS+"""
        try:
            shot_number = int(self.shot_entry.text())

            data = self.data_loader.load_data(shot_number)

            cache_key = f'{shot_number}_MSE'
            self.data[cache_key] = data

            self.available_listbox.clear()

            # Get good channels for TGAMMA
            tgamma_meas = data.measurements['tgamma']
            good_mask = tgamma_meas['good_mask']
            R_raw = tgamma_meas['R']
            raw_channels = tgamma_meas.get('raw_channels', list(range(1, 26)))

            # Add TGAMMA channels (good channels only)
            for ch_idx in np.where(good_mask)[0]:
                R_val = R_raw[ch_idx]
                ch_num = raw_channels[ch_idx] if ch_idx < len(raw_channels) else ch_idx + 1
                item_str = f'{shot_number:06d}_{R_val*1e3:.0f} (MSE{ch_num:02d})'
                self.available_listbox.addItem(item_str)

            print(f"[MSE] Data loaded: {np.sum(good_mask)} good channels")
            print(f"[MSE]   NB source: {data.nb_source}")

        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load data: {str(e)}")

    def _parse_entry(self, entry):
        """Parse time trace entry: XXXXXX_YYYY (MSE##)"""
        # Format: 036053_1800 (MSE01)
        main_part = entry.split('(')[0].strip()
        ch_label = entry.split('(')[1].split(')')[0].strip()
        parts = main_part.split('_')

        shot_number = int(parts[0])
        radius = float(parts[1]) / 1e3  # mm to m

        return shot_number, radius, ch_label

    def plot_data(self):
        """Plot time traces"""
        self.ax1.clear()
        self.ax2.clear()

        self.ax1.set_ylabel(r'$\gamma$ [rad]')
        self.ax2.set_xlabel('Time [s]')

        param = self._get_selected_param()
        if param == 'q':
            self.ax2.set_ylabel('q')
        else:
            self.ax2.set_ylabel('j [MA/m$^2$]')

        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            return

        gamma_min, gamma_max, param_max, param_min = 0, 0, 0, 0
        colors = self.plot_manager.color_manager.get_colors_for_entries(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                shot_number, radius, ch_label = self._parse_entry(entry)
                color = colors[i]

                cache_key = f'{shot_number}_MSE'

                if cache_key not in self.data:
                    data = self.data_loader.load_data(shot_number)
                    self.data[cache_key] = data
                else:
                    data = self.data[cache_key]

                # Find channel index by R position
                tgamma_meas = data.measurements['tgamma']
                R_raw = tgamma_meas['R']
                ch_idx = np.argmin(np.abs(R_raw - radius))
                actual_R = R_raw[ch_idx]

                # Plot TGAMMA time trace with fill_between
                tgamma_trace = tgamma_meas['data'][ch_idx, :]
                sgamma_trace = tgamma_meas['error'][ch_idx, :]

                label = f'#{shot_number} {actual_R*1e3:.0f}mm ({ch_label})'

                self.ax1.plot(data.time, tgamma_trace, '-', color=color,
                             linewidth=0.8, label=label)
                self.ax1.fill_between(data.time,
                                      tgamma_trace - sgamma_trace*0.5,
                                      tgamma_trace + sgamma_trace*0.5,
                                      color=color, alpha=0.3)

                gamma_min = min(gamma_min, np.nanpercentile(tgamma_trace, 2))
                gamma_max = max(gamma_max, np.nanpercentile(tgamma_trace, 98))

                # Plot j or q time trace interpolated at gamma R position
                param_meas = data.measurements[param]
                roa_data = param_meas['roa']
                param_data_all = param_meas['data']
                param_err_all = param_meas['error']

                # Interpolate at each time step
                param_trace = np.zeros(len(data.time_prof))
                param_err_trace = np.zeros(len(data.time_prof))

                for t_idx in range(len(data.time_prof)):
                    magx = data.measurements['magx']['data'][t_idx]
                    r_edge = self._get_r_edge_at_time(data, data.time_prof[t_idx])

                    # Convert r/a to R for this time step
                    roa_profile = roa_data[:, t_idx]
                    R_profile = magx + roa_profile * (r_edge - magx)

                    # Interpolate param at actual_R
                    param_profile = param_data_all[:, t_idx]
                    param_err_profile = param_err_all[:, t_idx]

                    # Sort by R for interpolation
                    sort_idx = np.argsort(R_profile)
                    R_sorted = R_profile[sort_idx]
                    param_sorted = param_profile[sort_idx]
                    param_err_sorted = param_err_profile[sort_idx]

                    param_trace[t_idx] = np.interp(actual_R, R_sorted, param_sorted)
                    param_err_trace[t_idx] = np.interp(actual_R, R_sorted, param_err_sorted)

                label2 = f'#{shot_number} {actual_R*1e3:.0f}mm ({ch_label})'

                self.ax2.plot(data.time_prof, param_trace, '-', color=color,
                             linewidth=0.8, label=label2)
                self.ax2.fill_between(data.time_prof,
                                      param_trace - param_err_trace*0.5,
                                      param_trace + param_err_trace*0.5,
                                      color=color, alpha=0.3)

                param_max = max(param_max, np.nanmax(param_trace))
                param_min = min(param_min, np.nanmin(param_trace))

            except Exception as e:
                print(f"[MSE] Error plotting {entry}: {str(e)}")

        # Set limits and styling
        self.ax1.axhline(0, color='gray', linestyle='--', alpha=0.5)
        if gamma_max > 0:
            gamma_margin = (gamma_max - gamma_min) * 0.1
            self.ax1.set_ylim(gamma_min - gamma_margin, gamma_max + gamma_margin)

        if param == 'q':
            self.ax2.axhline(1, color='gray', linestyle='--', alpha=0.5)
        else:
            self.ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
        self.ax2.set_ylim(0, 6)

        self._finalize_plot()

    def _write_data_to_file(self, file_path, selected_entries):
        """Write MSE time trace data to text file (gamma, q, j at same R position)"""
        with open(file_path, 'w') as f:
            # Write header
            f.write("# MSE Time Trace Data (q and j interpolated at TGAMMA position)\n")
            f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "gamma", "gamma_err", "q", "q_err", "j[MA/m2]", "j_err"
            ))

            for entry in selected_entries:
                shot_number, radius, ch_label = self._parse_entry(entry)
                cache_key = f'{shot_number}_MSE'

                if cache_key not in self.data:
                    data = self.data_loader.load_data(shot_number)
                    self.data[cache_key] = data
                else:
                    data = self.data[cache_key]

                # Find channel index by R position
                tgamma_meas = data.measurements['tgamma']
                R_raw = tgamma_meas['R']
                ch_idx = np.argmin(np.abs(R_raw - radius))
                actual_R = R_raw[ch_idx]

                tgamma_trace = tgamma_meas['data'][ch_idx, :]
                sgamma_trace = tgamma_meas['error'][ch_idx, :]

                # Get q and j measurements
                q_meas = data.measurements['q']
                j_meas = data.measurements['j']
                roa_data = q_meas['roa']

                # Use profile time (same as q/j time)
                for i in range(len(data.time_prof)):
                    # Find closest raw time index for gamma
                    raw_idx = np.argmin(np.abs(data.time - data.time_prof[i]))

                    magx = data.measurements['magx']['data'][i]
                    r_edge = self._get_r_edge_at_time(data, data.time_prof[i])

                    # Convert r/a to R for this time step
                    roa_profile = roa_data[:, i]
                    R_profile = magx + roa_profile * (r_edge - magx)

                    # Sort by R for interpolation
                    sort_idx = np.argsort(R_profile)
                    R_sorted = R_profile[sort_idx]

                    # Interpolate q
                    q_val = np.interp(actual_R, R_sorted, q_meas['data'][:, i][sort_idx])
                    q_err_val = np.interp(actual_R, R_sorted, q_meas['error'][:, i][sort_idx])

                    # Interpolate j
                    j_val = np.interp(actual_R, R_sorted, j_meas['data'][:, i][sort_idx])
                    j_err_val = np.interp(actual_R, R_sorted, j_meas['error'][:, i][sort_idx])

                    f.write(" %10d,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f\n" % (
                        shot_number, data.time_prof[i], actual_R,
                        tgamma_trace[raw_idx], sgamma_trace[raw_idx],
                        q_val, q_err_val, j_val, j_err_val
                    ))
