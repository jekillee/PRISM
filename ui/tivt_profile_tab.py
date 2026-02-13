"""
Ti/vT Profile tab with CES and XICS integration
"""

import numpy as np
from scipy.interpolate import interp1d

from PySide6.QtWidgets import (
    QMessageBox, QLineEdit, QComboBox, QPushButton,
    QWidget, QHBoxLayout, QVBoxLayout, QGroupBox, QGridLayout, QLabel, QFileDialog, QStyle
)

from ui.profile_base_tab import ProfileBaseTab
from ui.ui_constants import get_icon


class TiVTProfileTab(ProfileBaseTab):
    """Ti/vT Profile tab with CES and XICS support"""

    TAB_NAME = "Ti/vT"

    def __init__(self, parent, app_config, diagnostic_name, tab_type,
                 data_loader, efit_loader, plot_manager, file_parser):
        super().__init__(parent, app_config, diagnostic_name, tab_type,
                        data_loader, efit_loader, plot_manager)
        self.file_parser = file_parser
        self.xics_loader = None
        self.xics_data_cache = {}

    def _create_shot_input(self, parent):
        """Create data loading section with analysis type selection"""
        group = QGroupBox("1. Load Ti, vT Data")
        grid = QGridLayout(group)
        grid.setColumnStretch(1, 1)

        # Row 0: Shot label, entry, up/down buttons, dropdown, Fetch, File
        grid.addWidget(QLabel('Shot'), 0, 0)

        self.shot_entry = QLineEdit()
        self.shot_entry.setMinimumWidth(65)
        grid.addWidget(self.shot_entry, 0, 1)
        self.shot_entry.returnPressed.connect(self.load_shot_data)

        btn_updown = QWidget()
        btn_layout = QVBoxLayout(btn_updown)
        btn_layout.setContentsMargins(0, 0, 0, 0)
        btn_layout.setSpacing(0)
        mini_btn_style = "padding: 0px; border-radius: 2px;"
        up_btn = QPushButton()
        up_btn.setIcon(get_icon(QStyle.SP_ArrowUp))
        up_btn.setFixedSize(24, 15)
        up_btn.setStyleSheet(mini_btn_style)
        up_btn.clicked.connect(lambda: self._adjust_shot(1))
        btn_layout.addWidget(up_btn)
        down_btn = QPushButton()
        down_btn.setIcon(get_icon(QStyle.SP_ArrowDown))
        down_btn.setFixedSize(24, 15)
        down_btn.setStyleSheet(mini_btn_style)
        down_btn.clicked.connect(lambda: self._adjust_shot(-1))
        btn_layout.addWidget(down_btn)
        grid.addWidget(btn_updown, 0, 2)

        analysis_types = list(self.diag_config['analysis_types'].keys())
        self.analysis_type_combo = QComboBox()
        self.analysis_type_combo.addItems(analysis_types)
        self.analysis_type_combo.setCurrentText('mod')
        self.analysis_type_combo.setFixedWidth(70)
        grid.addWidget(self.analysis_type_combo, 0, 3)

        btn_frame = QWidget()
        btn_frame_layout = QHBoxLayout(btn_frame)
        btn_frame_layout.setContentsMargins(0, 0, 0, 0)

        self.fetch_button = QPushButton('Fetch')
        self.fetch_button.setFixedWidth(70)
        self.fetch_button.clicked.connect(self.load_shot_data)
        btn_frame_layout.addWidget(self.fetch_button)

        file_button = QPushButton()
        file_button.setIcon(get_icon(QStyle.SP_DirOpenIcon))
        file_button.setFixedWidth(28)
        file_button.clicked.connect(self.load_file_data)
        btn_frame_layout.addWidget(file_button)

        grid.addWidget(btn_frame, 0, 4)

        parent.layout().addWidget(group)

    def _adjust_shot(self, delta):
        """Adjust shot number by delta"""
        try:
            current = int(self.shot_entry.text())
            new_shot = max(1, current + delta)
            self.shot_entry.setText(str(new_shot))
        except ValueError:
            pass

    def _get_secondary_loader(self):
        """Get XICS loader for secondary diagnostic overlay"""
        return self._get_xics_loader()

    def _get_xics_loader(self):
        """Lazy initialization of XICS loader"""
        if self.xics_loader is None:
            from data_loaders.xics_loader import XICSLoader
            from config.diagnostic_config import DIAGNOSTICS
            self.xics_loader = XICSLoader(self.app_config, DIAGNOSTICS.get('XICS', {}))
        return self.xics_loader

    def load_shot_data(self):
        """Load CES and XICS shot data from MDS+"""
        try:
            shot_number = int(self.shot_entry.text())
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
            return

        analysis_type = self.analysis_type_combo.currentText()
        ces_loaded = False
        xics_loaded = False

        # Try to load CES data with selected analysis type only
        try:
            data = self.data_loader.load_data(shot_number, analysis_type)
            cache_key = f'{shot_number}_{analysis_type}'
            self.data[cache_key] = data
            ces_loaded = True
            print(f"[CES] {analysis_type} Data loaded: {len(data.radius)} channels, {len(data.time)} timepoints")
        except Exception as e:
            print(f"[CES] {analysis_type} Not available: {str(e)}")

        # Try to load XICS data
        try:
            loader = self._get_xics_loader()
            xics_data = loader.load_data(shot_number)
            if xics_data is not None:
                xics_cache_key = f'{shot_number}_XICS'
                self.xics_data_cache[xics_cache_key] = xics_data
                xics_loaded = True
                print(f"[XICS] Data loaded: {len(xics_data.time)} timepoints")
        except Exception as e:
            print(f"[XICS] Not available: {str(e)}")

        # Check if any data loaded
        if not ces_loaded and not xics_loaded:
            QMessageBox.critical(self.frame, "Error", f"No CES ({analysis_type}) or XICS data available for shot #{shot_number}")
            return

        # Update listbox
        self.available_listbox.clear()

        if ces_loaded:
            # Show CES timepoints
            for tp in data.time:
                item_str = f'{shot_number:06d}_{tp*1e3:06.0f} ({analysis_type})'
                self.available_listbox.addItem(item_str)
        elif xics_loaded:
            # Show XICS timepoints (fallback when no CES)
            for tp in xics_data.time:
                item_str = f'{shot_number:06d}_{tp*1e3:06.0f} (XICS)'
                self.available_listbox.addItem(item_str)

    def load_file_data(self):
        """Load CES data from result file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self.frame,
            "Select CES Result File",
            self.app_config.CES_RESULT_PATH,
            "Text files (*.txt);;All Files (*.*)"
        )

        if not file_path:
            return

        try:
            data, shot_number = self.file_parser.parse_file(file_path)

            file_key = f'file_{shot_number}_{data.source}'
            self.data[file_key] = data

            self.available_listbox.clear()

            for tp in data.time:
                item_str = f'{shot_number:06d}_{tp*1e3:06.0f} ({data.source})'
                self.available_listbox.addItem(item_str)

            print(f"[CES] Result file loaded: {data.source}")

            # Try to load XICS data for this shot
            loader = self._get_xics_loader()
            xics_data = loader.load_data(shot_number)

            if xics_data is not None:
                xics_cache_key = f'{shot_number}_XICS'
                self.xics_data_cache[xics_cache_key] = xics_data
                print(f"[XICS] Data loaded: {len(xics_data.time)} timepoints")

        except ValueError as e:
            QMessageBox.critical(self.frame, "Invalid File Format", str(e))
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load file: {str(e)}")

    def _parse_entry(self, entry):
        """Parse entry to get shot, time, and source"""
        if '(' in entry and ')' in entry:
            main_part = entry.split('(')[0].strip()
            source = entry.split('(')[1].split(')')[0].strip()
            parts = main_part.split('_')
        else:
            parts = entry.split('_')
            source = 'mod'

        shot_number = int(parts[0])
        time_point = float(parts[1]) / 1e3  # ms to s

        return shot_number, time_point, source

    def _parse_entry_for_efit(self, entry):
        """Parse entry for EFIT (extract shot and time)"""
        shot_number, time_point, _ = self._parse_entry(entry)
        return shot_number, time_point

    def _get_xics_at_time(self, shot_number, time_point):
        """Get XICS data at specific time point"""
        xics_cache_key = f'{shot_number}_XICS'

        if xics_cache_key not in self.xics_data_cache:
            return None

        xics_data = self.xics_data_cache[xics_cache_key]

        # Find nearest time index
        time_idx = np.argmin(np.abs(xics_data.time - time_point))

        # Check if within reasonable range (50ms)
        if np.abs(xics_data.time[time_idx] - time_point) > 0.05:
            return None

        Ti_data, Ti_err = xics_data.get_parameter('Ti')
        vT_data, vT_err = xics_data.get_parameter('vT')

        return {
            'R': xics_data.radius[0],
            'Ti': Ti_data[0, time_idx],
            'Ti_err': Ti_err[0, time_idx],
            'vT': vT_data[0, time_idx],
            'vT_err': vT_err[0, time_idx]
        }

    def _get_marker_style(self, source):
        """Get marker style based on data source"""
        if source == 'mod':
            # Filled circle
            return {'marker': 'o', 'markerfacecolor': None, 'markeredgewidth': 1, 'fillstyle': 'full'}
        elif source == 'nn':
            # Empty circle
            return {'marker': 'o', 'markerfacecolor': 'none', 'markeredgewidth': 1.5, 'fillstyle': 'full'}
        elif source == 'XICS':
            # Filled square
            return {'marker': 's', 'markerfacecolor': None, 'markeredgewidth': 1, 'fillstyle': 'full'}
        else:
            # File data: half-filled circle
            return {'marker': 'o', 'markerfacecolor': None, 'markeredgewidth': 1, 'fillstyle': 'right'}

    def _get_channel_info(self):
        """Return CES channel info based on selected listbox entries."""
        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]

        # Check if any selected entry uses CES data (not XICS-only)
        has_ces = False
        for entry in selected_entries:
            _, _, source = self._parse_entry(entry)
            if source != 'XICS':
                has_ces = True
                break

        if not has_ces:
            return []

        for key, data in self.data.items():
            if hasattr(data, 'radius') and len(data.radius) > 0:
                return [(f"CES_{i}", f"CES Ch{i+1}")
                        for i in range(len(data.radius))]
        return []

    def plot_data(self):
        """Plot R profiles with CES and XICS data"""
        self.ax1.clear()
        self.ax2.clear()

        self.ax1.set_xlabel('R [m]')
        self.ax2.set_xlabel('R [m]')
        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_ylabel(self.param2['label'])

        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            return

        ti_max, vt_max, vt_min = 0, 0, 0
        colors = self._get_plot_colors(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                shot_number, time_point, source = self._parse_entry(entry)
                color = colors[i]

                if source == 'XICS':
                    # XICS-only entry (no CES data)
                    xics_point = self._get_xics_at_time(shot_number, time_point)
                    if xics_point is not None:
                        xics_style = self._get_marker_style('XICS')
                        xics_label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (XICS)'

                        ti_max = max(ti_max, xics_point['Ti'])
                        vt_min = min(vt_min, xics_point['vT'])
                        vt_max = max(vt_max, xics_point['vT'])

                        self.ax1.plot(xics_point['R'], xics_point['Ti'],
                                     marker=xics_style['marker'], color=color,
                                     markersize=5, label=xics_label,
                                     linestyle='none')
                        self.ax2.plot(xics_point['R'], xics_point['vT'],
                                     marker=xics_style['marker'], color=color,
                                     markersize=5, label=xics_label,
                                     linestyle='none')

                        # Add XICS channel label (TXCS node)
                        self._add_channel_labels(self.ax1, [xics_point['R']], [xics_point['Ti']], 'TXCS_TI0', [53])
                        self._add_channel_labels(self.ax2, [xics_point['R']], [xics_point['vT']], 'TXCS_VR0', [53])
                    continue

                # CES data processing
                if source in ['mod', 'nn']:
                    cache_key = f'{shot_number}_{source}'
                else:
                    cache_key = f'file_{shot_number}_{source}'

                # Load data if not cached
                if cache_key not in self.data:
                    if source in ['mod', 'nn']:
                        data = self.data_loader.load_data(shot_number, source)
                    else:
                        continue  # File data should already be loaded
                    self.data[cache_key] = data
                else:
                    data = self.data[cache_key]

                time_idx = np.argmin(np.abs(data.time - time_point))
                R_data = data.radius

                Ti_data, Ti_err = data.get_parameter('Ti')
                vT_data, vT_err = data.get_parameter('vT')

                Ti_profile = Ti_data[:, time_idx]
                Ti_err_profile = Ti_err[:, time_idx]
                vT_profile = vT_data[:, time_idx]
                vT_err_profile = vT_err[:, time_idx]

                valid_idx = np.argmin(np.abs(R_data - self.app_config.R_EDGE))
                ti_max = max(ti_max, np.nanpercentile(Ti_profile[:valid_idx], 98))
                vt_min = min(vt_min, np.nanpercentile(vT_profile[:valid_idx], 2))
                vt_max = max(vt_max, np.nanpercentile(vT_profile[:valid_idx], 98))

                label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms ({source})'

                # Get marker style
                style = self._get_marker_style(source)

                plot_kwargs = {
                    'fmt': style['marker'],
                    'capsize': 5,
                    'label': label,
                    'color': color,
                    'markersize': 5,
                    'fillstyle': style['fillstyle']
                }
                if style['markerfacecolor'] == 'none':
                    plot_kwargs['markerfacecolor'] = 'none'
                    plot_kwargs['markeredgewidth'] = style['markeredgewidth']

                # Split enabled/disabled channels
                ch_keys = [f"CES_{j}" for j in range(len(R_data))]
                mask = self._get_channel_mask(ch_keys)

                if mask.any():
                    self.ax1.errorbar(R_data[mask], Ti_profile[mask], Ti_err_profile[mask], **plot_kwargs)
                    self.ax2.errorbar(R_data[mask], vT_profile[mask], vT_err_profile[mask], **plot_kwargs)
                if (~mask).any():
                    dim_kwargs = {**plot_kwargs, 'color': (0.6, 0.6, 0.6, 0.35), 'label': ''}
                    self.ax1.errorbar(R_data[~mask], Ti_profile[~mask], Ti_err_profile[~mask], **dim_kwargs)
                    self.ax2.errorbar(R_data[~mask], vT_profile[~mask], vT_err_profile[~mask], **dim_kwargs)

                # Add channel labels if enabled (CES_TI or CESNN_TI)
                n_channels = len(R_data)
                channels = list(range(1, n_channels + 1))
                node_prefix = 'CESNN_TI' if source == 'nn' else 'CES_TI'
                self._add_channel_labels(self.ax1, R_data, Ti_profile, node_prefix, channels)
                self._add_channel_labels(self.ax2, R_data, vT_profile, node_prefix.replace('TI', 'VT'), channels)

                # Plot XICS data at the same time point
                xics_point = self._get_xics_at_time(shot_number, time_point)
                if xics_point is not None:
                    xics_style = self._get_marker_style('XICS')
                    xics_label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (XICS)'

                    # Update max values
                    ti_max = max(ti_max, xics_point['Ti'])
                    vt_min = min(vt_min, xics_point['vT'])
                    vt_max = max(vt_max, xics_point['vT'])

                    # Plot XICS Ti (no error bars, closed square)
                    self.ax1.plot(xics_point['R'], xics_point['Ti'],
                                 marker=xics_style['marker'], color=color,
                                 markersize=5, label=xics_label,
                                 linestyle='none')

                    # Plot XICS vT (no error bars, closed square)
                    self.ax2.plot(xics_point['R'], xics_point['vT'],
                                 marker=xics_style['marker'], color=color,
                                 markersize=5, label=xics_label,
                                 linestyle='none')

                    # Add XICS channel label (ch 53)
                    self._add_channel_labels(self.ax1, [xics_point['R']], [xics_point['Ti']], 'TXCS_TI0', [53])
                    self._add_channel_labels(self.ax2, [xics_point['R']], [xics_point['vT']], 'TXCS_TI0', [53])

            except Exception as e:
                print(f"[Ti/vT] Error plotting {entry}: {str(e)}")

        # Set limits
        self.ax1.set_xlim(self.app_config.R_LIMITS)
        self.ax2.set_xlim(self.app_config.R_LIMITS)
        ti_margin = ti_max * 0.1
        self.ax1.set_ylim(0, ti_max + ti_margin)

        vt_margin = (vt_max - vt_min) * 0.1
        self.ax2.set_ylim(vt_min - vt_margin, vt_max + vt_margin)
        self.ax2.axhline(y=0, c='silver', ls='--')

        self.plot_manager.apply_common_styling(self.ax1, self.ax2)
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

        ti_max, vt_max = 0, 0
        colors = self._get_plot_colors(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                shot_number, time_point, source = self._parse_entry(entry)

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

                if source == 'XICS':
                    # XICS-only entry
                    xics_point = self._get_xics_at_time(shot_number, time_point)
                    if xics_point is not None:
                        xics_x = interp_func(xics_point['R'])
                        xics_style = self._get_marker_style('XICS')
                        xics_label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (XICS)'

                        ti_max = max(ti_max, xics_point['Ti'])
                        vt_max = max(vt_max, np.abs(xics_point['vT']))

                        self.ax1.plot(xics_x, xics_point['Ti'],
                                     marker=xics_style['marker'], color=color,
                                     markersize=5, label=xics_label,
                                     linestyle='none')
                        self.ax2.plot(xics_x, xics_point['vT'],
                                     marker=xics_style['marker'], color=color,
                                     markersize=5, label=xics_label,
                                     linestyle='none')

                        # Add XICS channel label (ch 53)
                        self._add_channel_labels(self.ax1, [xics_x], [xics_point['Ti']], 'TXCS_TI0', [53])
                        self._add_channel_labels(self.ax2, [xics_x], [xics_point['vT']], 'TXCS_TI0', [53])
                    continue

                # CES data processing
                if source in ['mod', 'nn']:
                    cache_key = f'{shot_number}_{source}'
                else:
                    cache_key = f'file_{shot_number}_{source}'

                data = self.data[cache_key]

                time_idx = np.argmin(np.abs(data.time - time_point))
                x_data = interp_func(data.radius)

                Ti_data, Ti_err = data.get_parameter('Ti')
                vT_data, vT_err = data.get_parameter('vT')

                Ti_profile = Ti_data[:, time_idx]
                Ti_err_profile = Ti_err[:, time_idx]
                vT_profile = vT_data[:, time_idx]
                vT_err_profile = vT_err[:, time_idx]

                lcfs_idx = np.argmin(np.abs(x_data - 1))
                ti_max = max(ti_max, np.nanmax(Ti_profile[:lcfs_idx]))
                vt_max = max(vt_max, np.nanmax(np.abs(vT_profile[:lcfs_idx])))

                label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms ({source})'

                style = self._get_marker_style(source)

                plot_kwargs = {
                    'fmt': style['marker'],
                    'capsize': 5,
                    'label': label,
                    'color': color,
                    'markersize': 5,
                    'fillstyle': style['fillstyle']
                }
                if style['markerfacecolor'] == 'none':
                    plot_kwargs['markerfacecolor'] = 'none'
                    plot_kwargs['markeredgewidth'] = style['markeredgewidth']

                # Split enabled/disabled channels
                ch_keys = [f"CES_{j}" for j in range(len(x_data))]
                mask = self._get_channel_mask(ch_keys)

                if mask.any():
                    self.ax1.errorbar(x_data[mask], Ti_profile[mask], Ti_err_profile[mask], **plot_kwargs)
                    self.ax2.errorbar(x_data[mask], vT_profile[mask], vT_err_profile[mask], **plot_kwargs)
                if (~mask).any():
                    dim_kwargs = {**plot_kwargs, 'color': (0.6, 0.6, 0.6, 0.35), 'label': ''}
                    self.ax1.errorbar(x_data[~mask], Ti_profile[~mask], Ti_err_profile[~mask], **dim_kwargs)
                    self.ax2.errorbar(x_data[~mask], vT_profile[~mask], vT_err_profile[~mask], **dim_kwargs)

                # Add channel labels if enabled
                n_channels = len(x_data)
                channels = list(range(1, n_channels + 1))
                self._add_channel_labels(self.ax1, x_data, Ti_profile, source, channels)
                self._add_channel_labels(self.ax2, x_data, vT_profile, source, channels)

                # Plot XICS data at the same time point
                xics_point = self._get_xics_at_time(shot_number, time_point)
                if xics_point is not None:
                    xics_x = interp_func(xics_point['R'])
                    xics_style = self._get_marker_style('XICS')
                    xics_label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (XICS)'

                    ti_max = max(ti_max, xics_point['Ti'])
                    vt_max = max(vt_max, np.abs(xics_point['vT']))

                    # Plot XICS Ti (no error bars, closed square)
                    self.ax1.plot(xics_x, xics_point['Ti'],
                                 marker=xics_style['marker'], color=color,
                                 markersize=5, label=xics_label,
                                 linestyle='none')

                    # Plot XICS vT (no error bars, closed square)
                    self.ax2.plot(xics_x, xics_point['vT'],
                                 marker=xics_style['marker'], color=color,
                                 markersize=5, label=xics_label,
                                 linestyle='none')

                    # Add XICS channel label (ch 53)
                    self._add_channel_labels(self.ax1, [xics_x], [xics_point['Ti']], 'TXCS_TI0', [53])
                    self._add_channel_labels(self.ax2, [xics_x], [xics_point['vT']], 'TXCS_TI0', [53])

            except Exception as e:
                print(f"[Ti/vT] Error plotting {entry}: {str(e)}")

        self.ax1.set_xlabel(x_label)
        self.ax2.set_xlabel(x_label)
        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_ylabel(self.param2['label'])

        self.ax1.set_xlim(0, 1.05)
        self.ax2.set_xlim(0, 1.05)
        self.ax1.set_ylim(0, ti_max * 1.1)
        self.ax2.set_ylim(-vt_max * 0.1, vt_max * 1.1)
        self.ax2.axhline(y=0, c='silver', ls='--')

        for ax in [self.ax1, self.ax2]:
            ax.axvline(x=0, c='k', ls='--')
            ax.axvline(x=1, c='k', ls='--')

        self.plot_manager.apply_common_styling(self.ax1, self.ax2)
        self.canvas.draw()

        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()

    def _write_data_to_file(self, file_path, selected_entries):
        """Write Ti/vT profile data to text file"""
        with open(file_path, 'w') as f:
            # Write header
            f.write("# Ti/vT Profile Data\n")
            if self.computed_efit_tree:
                f.write(f"# EFIT Source: {self.computed_efit_tree}\n")
            f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "psi_N", "rho_pol", "rho_tor",
                "Ti[keV]", "Ti_err", "vT[km/s]", "vT_err", "Source"
            ))

            for entry in selected_entries:
                shot_number, time_point, source = self._parse_entry(entry)

                # Handle XICS-only entries
                if source == 'XICS':
                    xics_point = self._get_xics_at_time(shot_number, time_point)
                    if xics_point is not None:
                        xics_psi_n, xics_rho_pol, xics_rho_tor = self._get_efit_values_at_R(
                            entry, np.array([xics_point['R']]))

                        f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10.3f,%10.3f,%10.3f,%10s\n" % (
                            shot_number, time_point, xics_point['R'],
                            self._format_value(xics_psi_n[0]),
                            self._format_value(xics_rho_pol[0]),
                            self._format_value(xics_rho_tor[0]),
                            xics_point['Ti'], xics_point['Ti_err'],
                            xics_point['vT'], xics_point['vT_err'], 'XICS'
                        ))
                    continue

                # Determine cache key for CES data
                if source in ['mod', 'nn']:
                    cache_key = f'{shot_number}_{source}'
                else:
                    cache_key = f'file_{shot_number}_{source}'

                # Load data if not cached
                if cache_key not in self.data:
                    if source in ['mod', 'nn']:
                        data = self.data_loader.load_data(shot_number, source)
                        self.data[cache_key] = data
                    else:
                        continue
                else:
                    data = self.data[cache_key]

                time_idx = np.argmin(np.abs(data.time - time_point))
                actual_time = data.time[time_idx]
                R_data = data.radius

                Ti_data, Ti_err = data.get_parameter('Ti')
                vT_data, vT_err = data.get_parameter('vT')

                Ti_profile = Ti_data[:, time_idx]
                Ti_err_profile = Ti_err[:, time_idx]
                vT_profile = vT_data[:, time_idx]
                vT_err_profile = vT_err[:, time_idx]

                # Get EFIT values
                psi_n, rho_pol, rho_tor = self._get_efit_values_at_R(entry, R_data)

                # Write CES data rows
                for i in range(len(R_data)):
                    f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10.3f,%10.3f,%10.3f,%10s\n" % (
                        shot_number, actual_time, R_data[i],
                        self._format_value(psi_n[i]),
                        self._format_value(rho_pol[i]),
                        self._format_value(rho_tor[i]),
                        Ti_profile[i], Ti_err_profile[i],
                        vT_profile[i], vT_err_profile[i], source
                    ))

                # Write XICS data if available (for CES entries)
                xics_point = self._get_xics_at_time(shot_number, time_point)
                if xics_point is not None:
                    xics_psi_n, xics_rho_pol, xics_rho_tor = self._get_efit_values_at_R(
                        entry, np.array([xics_point['R']]))

                    f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10.3f,%10.3f,%10.3f,%10s\n" % (
                        shot_number, actual_time, xics_point['R'],
                        self._format_value(xics_psi_n[0]),
                        self._format_value(xics_rho_pol[0]),
                        self._format_value(xics_rho_tor[0]),
                        xics_point['Ti'], xics_point['Ti_err'],
                        xics_point['vT'], xics_point['vT_err'], 'XICS'
                    ))
