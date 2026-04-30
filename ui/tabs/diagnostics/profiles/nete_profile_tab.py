"""
ne, Te Profile tab with Thomson/ECE selection
"""

import numpy as np
from scipy.interpolate import interp1d

from PySide6.QtWidgets import (
    QMessageBox, QLineEdit, QComboBox, QPushButton, QApplication,
    QWidget, QHBoxLayout, QVBoxLayout, QGroupBox, QGridLayout, QLabel, QStyle
)

from ui.tabs.profile_base_tab import ProfileBaseTab
from ui.ui_constants import get_icon


class NeTeProfileTab(ProfileBaseTab):
    """ne, Te Profile tab supporting Thomson and/or ECE"""

    TAB_NAME = "ne/Te"

    # TCI tangent radii [m] and round-trip path lengths [m] (from KSTAR TCI spec)
    TCI_CHANNELS = {
        'TCI01': {'Rmid': 1.34, 'path_rt': 7.23},
        'TCI02': {'Rmid': 1.78, 'path_rt': 5.51},
        'TCI03': {'Rmid': 1.91, 'path_rt': 4.76},
        'TCI04': {'Rmid': 2.04, 'path_rt': 3.80},
        'TCI05': {'Rmid': 2.16, 'path_rt': 2.52},
    }
    TCI_REND = 2.25  # Typical R of LCFS [m] (MDS+ path length is based on this value)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ece_loader = None
        self.ece_data_cache = {}
        self.tci_loader = None

    def _get_secondary_loader(self):
        """Get ECE loader for secondary diagnostic overlay"""
        return self._get_ece_loader()

    def _create_shot_input(self, parent):
        """Create data loading section with diagnostic selection"""
        group = QGroupBox("1. Load ne, Te Data")
        grid = QGridLayout(group)
        grid.setColumnStretch(1, 1)

        # Row 0: Shot label, entry, up/down buttons, diagnostic dropdown, Fetch
        grid.addWidget(QLabel('Shot'), 0, 0)

        self.shot_entry = QLineEdit()
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

        diag_options = ['TS+ECE', 'TS', 'ECE (100Hz)', 'ECE (1kHz)']
        self.diag_combo = QComboBox()
        self.diag_combo.addItems(diag_options)
        self.diag_combo.setCurrentText('TS+ECE')
        self.diag_combo.setFixedWidth(110)
        grid.addWidget(self.diag_combo, 0, 3)

        self.fetch_button = QPushButton('Fetch')
        self.fetch_button.setFixedWidth(70)
        self.fetch_button.clicked.connect(self.load_shot_data)
        grid.addWidget(self.fetch_button, 0, 4)

        parent.layout().addWidget(group)

    def _adjust_shot(self, delta):
        """Adjust shot number by delta"""
        try:
            current = int(self.shot_entry.text())
            new_shot = max(1, current + delta)
            self.shot_entry.setText(str(new_shot))
        except ValueError:
            pass


    def _get_preview_info(self):
        if not hasattr(self, '_last_preview_data'):
            return None
        data, shot, source = self._last_preview_data
        # Include ECE data if available (TS+ECE mode)
        ece_data = None
        ece_key = f'{shot}_ECE_100Hz'
        if ece_key in self.ece_data_cache:
            ece_data = self.ece_data_cache[ece_key]
        return (data, shot, source,
                'Te', 'ne', self.param1['label'], self.param2['label'],
                ece_data)

    def _ece_channel_keys(self, ece_data):
        """Return channel-number-based keys for ECE data (e.g. ['ECE_ch62', 'ECE_ch07', ...])"""
        channels = ece_data.measurements['Te'].get('channels',
                   list(range(1, len(ece_data.radius) + 1)))
        return [f"ECE_ch{channels[j]}" for j in range(len(ece_data.radius))]

    def _get_ece_loader(self):
        """On-demand initialization of ECE loader"""
        if self.ece_loader is None:
            from data_loaders.ece_loader import ECELoader
            from config.diagnostic_config import DIAGNOSTICS
            self.ece_loader = ECELoader(self.app_config, DIAGNOSTICS['ECE'])
        return self.ece_loader

    def load_shot_data(self):
        """Load shot data based on diagnostic selection"""
        try:
            shot_number = int(self.shot_entry.text())
            selection = self.diag_combo.currentText()

            self._set_status(f"Loading #{shot_number} ({selection})...", 'blue')
            self.available_listbox.clear()

            # Parse selection to determine mode and ECE sampling
            if selection.startswith('ECE'):
                load_thomson = False
                load_ece = True
                if '100Hz' in selection:
                    sampling_rate = 0.01
                    sampling_key = '100Hz'
                else:  # 1kHz
                    sampling_rate = 0.001
                    sampling_key = '1kHz'
            elif selection == 'TS':
                load_thomson = True
                load_ece = False
                sampling_rate = None
                sampling_key = None
            else:  # TS+ECE
                load_thomson = True
                load_ece = True
                sampling_rate = 0.01
                sampling_key = '100Hz'

            if load_thomson:
                # Load Thomson data
                data = self.data_loader.load_data(shot_number)
                cache_key = f'{shot_number}_TS'
                self.data[cache_key] = data

                # Label depends on mode
                label = selection if selection == 'TS+ECE' else 'TS'

                for tp in data.time:
                    item_str = f'{shot_number:06d}_{tp*1e3:06.0f} ({label})'
                    self.available_listbox.addItem(item_str)

                self._last_preview_data = (data, shot_number, label)
                print(f"[Thomson] Data loaded: {len(data.radius)} channels, {len(data.time)} timepoints")

            if load_ece:
                loader = self._get_ece_loader()
                ece_data = loader.load_data(shot_number, sampling_rate=sampling_rate)

                cache_key = f'{shot_number}_ECE_{sampling_key}'
                self.ece_data_cache[cache_key] = ece_data

                # Add ECE timepoints to listbox only for "ECE only" mode
                if selection.startswith('ECE'):
                    for tp in ece_data.time:
                        item_str = f'{shot_number:06d}_{tp*1e3:06.0f} (ECE)'
                        self.available_listbox.addItem(item_str)

                if selection.startswith('ECE'):
                    self._last_preview_data = (ece_data, shot_number, 'ECE')
                n_valid = np.sum(ece_data.measurements['Te']['valid_mask'])
                n_overlap = np.sum(ece_data.measurements['Te']['overlap_mask'])
                print(f"[ECE] Data loaded: {len(ece_data.radius)} channels, {len(ece_data.time)} timepoints @ {sampling_key}")
                print(f"[ECE]   Valid: {n_valid}, Overlap: {n_overlap}")

            if hasattr(self, 'browse_button'):
                self.browse_button.setEnabled(True)
                self.browse_button.setText(f"Browse #{shot_number}")

            self._set_status(f"#{shot_number} loaded ({selection})", 'green')

        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
        except Exception as e:
            self._set_status(f"Failed to load #{shot_number}: {e}", 'red')
            QMessageBox.critical(self.frame, "Error", f"Failed to load data: {str(e)}")

    def _parse_entry(self, entry):
        """Parse entry to get shot, time, and source"""
        if '(' in entry and ')' in entry:
            main_part = entry.split('(')[0].strip()
            source = entry.split('(')[1].split(')')[0].strip()
            parts = main_part.split('_')
        else:
            parts = entry.split('_')
            source = 'TS'

        shot_number = int(parts[0])
        time_point = float(parts[1]) / 1e3  # ms to s

        return shot_number, time_point, source

    def _parse_entry_for_efit(self, entry):
        """Parse entry for EFIT (extract shot and time)"""
        shot_number, time_point, _ = self._parse_entry(entry)
        return shot_number, time_point

    def _get_ece_cache_key(self, shot_number):
        """Find ECE cache key for given shot"""
        for key in self.ece_data_cache.keys():
            if key.startswith(f'{shot_number}_ECE'):
                return key
        return None

    def _get_rshift_diagnostics(self):
        """Thomson and ECE support R-shift"""
        return ['Thomson', 'ECE']

    def _get_fit_param_types(self):
        """Return parameter types for Te/ne fitting"""
        return ['Te', 'ne']

    def _expand_entries_for_dt(self, selected_entries, dt_s):
        """Expand each TS entry to all time slices within [t-dt, t+dt]."""
        expanded = []
        for entry in selected_entries:
            shot_number, time_point, source = self._parse_entry(entry)
            if source in ['TS', 'TS+ECE']:
                cache_key = f'{shot_number}_TS'
                if cache_key in self.data:
                    data = self.data[cache_key]
                    t_mask = (data.time >= time_point - dt_s) & (data.time <= time_point + dt_s)
                    for t in data.time[t_mask]:
                        expanded.append(f'{shot_number:06d}_{t*1e3:06.0f} ({source})')
                    continue
            elif source == 'ECE':
                cache_key = self._get_ece_cache_key(shot_number)
                if cache_key and cache_key in self.ece_data_cache:
                    ece = self.ece_data_cache[cache_key]
                    t_mask = (ece.time >= time_point - dt_s) & (ece.time <= time_point + dt_s)
                    for t in ece.time[t_mask]:
                        expanded.append(f'{shot_number:06d}_{t*1e3:06.0f} ({source})')
                    continue
            expanded.append(entry)
        seen = set()
        return [e for e in expanded if not (e in seen or seen.add(e))]

    def _get_data_sampling_dt(self, entry):
        shot_number, time_point, source = self._parse_entry(entry)
        if source in ['TS', 'TS+ECE']:
            cache_key = f'{shot_number}_TS'
            if cache_key in self.data and len(self.data[cache_key].time) > 1:
                return float(np.median(np.diff(self.data[cache_key].time)))
        elif source == 'ECE':
            cache_key = self._get_ece_cache_key(shot_number)
            if cache_key and cache_key in self.ece_data_cache:
                ece = self.ece_data_cache[cache_key]
                if len(ece.time) > 1:
                    return float(np.median(np.diff(ece.time)))
        return None

    def _get_efit_profile_data(self, entry, interp_func, dt_s=0.0):
        """Get Te/ne profile data in EFIT coordinates for fitting"""
        shot_number, time_point, source = self._parse_entry(entry)

        if source == 'ECE':
            cache_key = self._get_ece_cache_key(shot_number)
            if cache_key is None:
                return None

            ece_data = self.ece_data_cache[cache_key]
            ece_rshift = self._get_rshift('ECE')
            x_data = interp_func(ece_data.radius + ece_rshift)
            Te_all = ece_data.measurements['Te']['data']
            valid_mask = ece_data.measurements['Te']['valid_mask']

            if dt_s > 0:
                t_mask = (ece_data.time >= time_point - dt_s) & (ece_data.time <= time_point + dt_s)
                n_avg = np.sum(t_mask)
                if n_avg > 0:
                    Te_profile = np.nanmean(Te_all[:, t_mask], axis=1)
                    print(f"[Fitting]   {entry}: ECE averaged {n_avg} slices")
                else:
                    Te_profile = Te_all[:, np.argmin(np.abs(ece_data.time - time_point))]
            else:
                Te_profile = Te_all[:, np.argmin(np.abs(ece_data.time - time_point))]

            ece_keys = self._ece_channel_keys(ece_data)
            ece_ch_mask = self._get_channel_mask(ece_keys)
            combined_mask = valid_mask & ece_ch_mask

            return {
                'Te': {'x': x_data[combined_mask], 'y': Te_profile[combined_mask], 'err': None},
            }

        # Thomson data (TS or TS+ECE)
        cache_key = f'{shot_number}_TS'
        if cache_key not in self.data:
            return None

        data = self.data[cache_key]
        ts_rshift = self._get_rshift('Thomson')
        x_data = interp_func(data.radius + ts_rshift)

        Te_all, Te_eu_all, Te_el_all = data.get_parameter_asymmetric('Te')
        ne_all, ne_eu_all, ne_el_all = data.get_parameter_asymmetric('ne')

        if dt_s > 0:
            t_mask = (data.time >= time_point - dt_s) & (data.time <= time_point + dt_s)
            n_avg = np.sum(t_mask)
            if n_avg > 0:
                Te_profile = np.nanmean(Te_all[:, t_mask], axis=1)
                Te_err = np.sqrt(np.nanmean(((Te_eu_all[:, t_mask] + Te_el_all[:, t_mask]) / 2.0) ** 2, axis=1) / n_avg)
                ne_profile = np.nanmean(ne_all[:, t_mask], axis=1)
                ne_err = np.sqrt(np.nanmean(((ne_eu_all[:, t_mask] + ne_el_all[:, t_mask]) / 2.0) ** 2, axis=1) / n_avg)
                print(f"[Fitting]   {entry}: TS averaged {n_avg} slices in [{time_point-dt_s:.3f}, {time_point+dt_s:.3f}]s")
            else:
                time_idx = np.argmin(np.abs(data.time - time_point))
                Te_profile = Te_all[:, time_idx]
                Te_err = (Te_eu_all[:, time_idx] + Te_el_all[:, time_idx]) / 2.0
                ne_profile = ne_all[:, time_idx]
                ne_err = (ne_eu_all[:, time_idx] + ne_el_all[:, time_idx]) / 2.0
        else:
            time_idx = np.argmin(np.abs(data.time - time_point))
            Te_profile = Te_all[:, time_idx]
            Te_err = (Te_eu_all[:, time_idx] + Te_el_all[:, time_idx]) / 2.0
            ne_profile = ne_all[:, time_idx]
            ne_err = (ne_eu_all[:, time_idx] + ne_el_all[:, time_idx]) / 2.0

        ts_keys = [f"TS_{j}" for j in range(len(data.radius))]
        ts_mask = self._get_channel_mask(ts_keys)

        result = {
            'Te': {'x': x_data[ts_mask], 'y': Te_profile[ts_mask], 'err': Te_err[ts_mask]},
            'ne': {'x': x_data[ts_mask], 'y': ne_profile[ts_mask], 'err': ne_err[ts_mask]},
        }

        # Include ECE data for TS+ECE mode
        if source == 'TS+ECE':
            ece_cache_key = self._get_ece_cache_key(shot_number)
            if ece_cache_key is not None:
                ece_data = self.ece_data_cache[ece_cache_key]
                ece_rshift = self._get_rshift('ECE')
                ece_x_data = interp_func(ece_data.radius + ece_rshift)
                ece_Te_all = ece_data.measurements['Te']['data']
                valid_mask = ece_data.measurements['Te']['valid_mask']

                if dt_s > 0:
                    ece_t_mask = (ece_data.time >= time_point - dt_s) & (ece_data.time <= time_point + dt_s)
                    n_ece = np.sum(ece_t_mask)
                    if n_ece > 0:
                        ece_Te = np.nanmean(ece_Te_all[:, ece_t_mask], axis=1)
                    else:
                        ece_Te = ece_Te_all[:, np.argmin(np.abs(ece_data.time - time_point))]
                else:
                    ece_Te = ece_Te_all[:, np.argmin(np.abs(ece_data.time - time_point))]

                ece_keys = self._ece_channel_keys(ece_data)
                ece_ch_mask = self._get_channel_mask(ece_keys)
                combined_mask = valid_mask & ece_ch_mask

                combined_x = np.concatenate([x_data[ts_mask], ece_x_data[combined_mask]])
                combined_y = np.concatenate([Te_profile[ts_mask], ece_Te[combined_mask]])
                combined_err = np.concatenate([Te_err[ts_mask], np.full(np.sum(combined_mask), np.nan)])

                result['Te'] = {'x': combined_x, 'y': combined_y, 'err': combined_err}

        return result

    def _get_channel_info(self):
        """Return channel info based on selected listbox entries."""
        info = []
        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]

        has_ts = False
        ece_shots = set()
        for entry in selected_entries:
            _, _, source = self._parse_entry(entry)
            if source in ('TS', 'TS+ECE'):
                has_ts = True
                if source == 'TS+ECE':
                    shot_number, _, _ = self._parse_entry(entry)
                    ece_shots.add(shot_number)
            elif source == 'ECE':
                shot_number, _, _ = self._parse_entry(entry)
                ece_shots.add(shot_number)

        # Thomson channels (all shots share the same channel structure)
        if has_ts:
            for key, data in self.data.items():
                if hasattr(data, 'radius') and len(data.radius) > 0:
                    for i in range(len(data.radius)):
                        if i < 14:
                            lbl = f"TS_CORE{i+1}"
                        else:
                            lbl = f"TS_EDGE{i-14+1}"
                        info.append((f"TS_{i}", lbl))
                    break

        # ECE channels from all selected shots (union, no duplicates)
        seen_ece_keys = set()
        for shot in ece_shots:
            cache_key = self._get_ece_cache_key(shot)
            if cache_key is None:
                continue
            ece_data = self.ece_data_cache[cache_key]
            valid_mask = ece_data.measurements['Te']['valid_mask']
            ece_channels = ece_data.measurements['Te'].get('channels',
                           list(range(1, len(ece_data.radius) + 1)))
            ch_keys = self._ece_channel_keys(ece_data)
            for i in range(len(ece_data.radius)):
                if valid_mask[i] and ch_keys[i] not in seen_ece_keys:
                    seen_ece_keys.add(ch_keys[i])
                    lbl = f"ECE{ece_channels[i]:02d}"
                    info.append((ch_keys[i], lbl))

        return info

    def plot_data(self):
        """Plot R profiles"""
        self.ax1.clear()
        self.ax2.clear()
        self._clear_click_points()

        self.ax1.set_xlabel('R [m]')
        self.ax2.set_xlabel('R [m]')
        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_ylabel(self.param2['label'])

        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            return

        te_max, ne_max = 0, 0
        colors = self._get_plot_colors(selected_entries)

        has_thomson = False

        for i, entry in enumerate(selected_entries):
            try:
                shot_number, time_point, source = self._parse_entry(entry)
                color = colors[i]

                if source in ['TS', 'TS+ECE']:
                    # Plot Thomson data
                    has_thomson = True
                    cache_key = f'{shot_number}_TS'
                    data = self.data.get(cache_key) or self.data_loader.load_data(shot_number)
                    if cache_key not in self.data:
                        self.data[cache_key] = data

                    time_idx = np.argmin(np.abs(data.time - time_point))
                    R_data = data.radius

                    Te_data, Te_err_upper, Te_err_lower = data.get_parameter_asymmetric('Te')
                    ne_data, ne_err_upper, ne_err_lower = data.get_parameter_asymmetric('ne')

                    Te_profile = Te_data[:, time_idx]
                    Te_err_upper_profile = Te_err_upper[:, time_idx]
                    Te_err_lower_profile = Te_err_lower[:, time_idx]
                    ne_profile = ne_data[:, time_idx]
                    ne_err_upper_profile = ne_err_upper[:, time_idx]
                    ne_err_lower_profile = ne_err_lower[:, time_idx]

                    # Split enabled/disabled Thomson channels
                    ts_keys = [f"TS_{j}" for j in range(len(R_data))]
                    ts_mask = self._get_channel_mask(ts_keys)

                    valid_idx = np.argmin(np.abs(R_data - self.app_config.R_EDGE))
                    vm = ts_mask[:valid_idx]
                    if vm.any():
                        te_max = max(te_max, np.nanmax(Te_profile[:valid_idx][vm]))
                        ne_max = max(ne_max, np.nanmax(ne_profile[:valid_idx][vm]))

                    label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (TS)'

                    if ts_mask.any():
                        self.ax1.errorbar(R_data[ts_mask], Te_profile[ts_mask],
                                         yerr=[Te_err_lower_profile[ts_mask], Te_err_upper_profile[ts_mask]],
                                         fmt='o', capsize=5, label=label,
                                         color=color, markersize=5, zorder=10)
                        self.ax2.errorbar(R_data[ts_mask], ne_profile[ts_mask],
                                         yerr=[ne_err_lower_profile[ts_mask], ne_err_upper_profile[ts_mask]],
                                         fmt='o', capsize=5, label=label,
                                         color=color, markersize=5, zorder=10)
                    if (~ts_mask).any():
                        self.ax1.errorbar(R_data[~ts_mask], Te_profile[~ts_mask],
                                         yerr=[Te_err_lower_profile[~ts_mask], Te_err_upper_profile[~ts_mask]],
                                         fmt='o', capsize=5, label='',
                                         color=(0.6, 0.6, 0.6, 0.35), markersize=5, zorder=5)
                        self.ax2.errorbar(R_data[~ts_mask], ne_profile[~ts_mask],
                                         yerr=[ne_err_lower_profile[~ts_mask], ne_err_upper_profile[~ts_mask]],
                                         fmt='o', capsize=5, label='',
                                         color=(0.6, 0.6, 0.6, 0.35), markersize=5, zorder=5)

                    # Add Thomson channel labels (TS_CORE1-14 for core, TS_EDGE1-17 for edge)
                    n_channels = len(R_data)
                    for ch_idx in range(n_channels):
                        if self.show_channel_checkbox.isChecked():
                            if ch_idx < 14:
                                lbl = f'TS_CORE{ch_idx + 1}'
                            else:
                                lbl = f'TS_EDGE{ch_idx - 14 + 1}'
                            self.ax1.annotate(lbl, (R_data[ch_idx], Te_profile[ch_idx]),
                                             textcoords='offset points',
                                             xytext=(0, 5), ha='center', fontsize=7, alpha=0.8,
                                             clip_on=True, annotation_clip=True)
                            self.ax2.annotate(lbl, (R_data[ch_idx], ne_profile[ch_idx]),
                                             textcoords='offset points',
                                             xytext=(0, 5), ha='center', fontsize=7, alpha=0.8,
                                             clip_on=True, annotation_clip=True)

                    # Register for double-click toggle
                    self._register_click_points(self.ax1, R_data, Te_profile, ts_keys)
                    self._register_click_points(self.ax2, R_data, ne_profile, ts_keys)

                    # Plot ECE data only for "TS+ECE" mode
                    if source == 'TS+ECE':
                        ece_cache_key = self._get_ece_cache_key(shot_number)
                        if ece_cache_key is not None:
                            ece_data = self.ece_data_cache[ece_cache_key]

                            ece_time_idx = np.argmin(np.abs(ece_data.time - time_point))
                            ece_R_data = ece_data.radius

                            ece_Te_data = ece_data.measurements['Te']['data']
                            ece_Te_profile = ece_Te_data[:, ece_time_idx]

                            valid_mask = ece_data.measurements['Te']['valid_mask']

                            # Split enabled/disabled ECE channels
                            ece_keys = self._ece_channel_keys(ece_data)
                            ece_ch_mask = self._get_channel_mask(ece_keys)

                            valid_in_range = valid_mask & ece_ch_mask & (ece_R_data <= self.app_config.R_EDGE)
                            if np.any(valid_in_range):
                                te_max = max(te_max, np.nanmax(ece_Te_profile[valid_in_range]))

                            ece_time_ms = entry.split("_")[1].split()[0]
                            ece_label = f'#{shot_number} {ece_time_ms}ms (ECE)'
                            # Combine with hardware valid mask
                            enabled = valid_mask & ece_ch_mask
                            disabled = valid_mask & (~ece_ch_mask)

                            if enabled.any():
                                self.ax1.plot(ece_R_data[enabled], ece_Te_profile[enabled],
                                             's', color=color, markersize=5,
                                             markerfacecolor='none', markeredgewidth=1.5,
                                             label=ece_label, zorder=5)
                            if disabled.any():
                                self.ax1.plot(ece_R_data[disabled], ece_Te_profile[disabled],
                                             's', color=(0.6, 0.6, 0.6, 0.35), markersize=5,
                                             markerfacecolor='none', markeredgewidth=1.5,
                                             label='', zorder=3)

                            # Add ECE channel labels (ECE01~76)
                            ece_channels = ece_data.measurements['Te'].get('channels',
                                           list(range(1, len(ece_R_data) + 1)))
                            valid_R = ece_R_data[valid_mask]
                            valid_Te = ece_Te_profile[valid_mask]
                            valid_ch = [ece_channels[j] for j in range(len(valid_mask)) if valid_mask[j]]
                            self._add_channel_labels(self.ax1, valid_R, valid_Te, 'ECE', valid_ch)

                            # Register ECE for double-click toggle
                            valid_ece_keys = [ece_keys[j] for j in range(len(valid_mask)) if valid_mask[j]]
                            self._register_click_points(self.ax1, valid_R, valid_Te, valid_ece_keys)

                elif source == 'ECE':
                    # Plot ECE data (valid channels only)
                    cache_key = self._get_ece_cache_key(shot_number)

                    if cache_key is None:
                        # Load with default sampling rate if not found
                        loader = self._get_ece_loader()
                        ece_data = loader.load_data(shot_number, sampling_rate=0.01)  # 100Hz default
                        cache_key = f'{shot_number}_ECE_100Hz'
                        self.ece_data_cache[cache_key] = ece_data
                    else:
                        ece_data = self.ece_data_cache[cache_key]

                    time_idx = np.argmin(np.abs(ece_data.time - time_point))
                    R_data = ece_data.radius

                    Te_data = ece_data.measurements['Te']['data']
                    Te_profile = Te_data[:, time_idx]

                    valid_mask = ece_data.measurements['Te']['valid_mask']

                    # Split enabled/disabled ECE channels
                    ece_keys = self._ece_channel_keys(ece_data)
                    ece_ch_mask = self._get_channel_mask(ece_keys)

                    valid_in_range = valid_mask & ece_ch_mask & (R_data <= self.app_config.R_EDGE)
                    if np.any(valid_in_range):
                        te_max = max(te_max, np.nanmax(Te_profile[valid_in_range]))

                    label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (ECE)'
                    enabled = valid_mask & ece_ch_mask
                    disabled = valid_mask & (~ece_ch_mask)

                    if enabled.any():
                        self.ax1.plot(R_data[enabled], Te_profile[enabled],
                                     's', color=color, markersize=5,
                                     markerfacecolor='none', markeredgewidth=1.5,
                                     label=label, zorder=5)
                    if disabled.any():
                        self.ax1.plot(R_data[disabled], Te_profile[disabled],
                                     's', color=(0.6, 0.6, 0.6, 0.35), markersize=5,
                                     markerfacecolor='none', markeredgewidth=1.5,
                                     label='', zorder=3)

                    # Add ECE channel labels
                    ece_channels = ece_data.measurements['Te'].get('channels',
                                   list(range(1, len(R_data) + 1)))
                    valid_R = R_data[valid_mask]
                    valid_Te = Te_profile[valid_mask]
                    valid_ch = [ece_channels[j] for j in range(len(valid_mask)) if valid_mask[j]]
                    self._add_channel_labels(self.ax1, valid_R, valid_Te, 'ECE', valid_ch)

                    # Register ECE for double-click toggle
                    valid_ece_keys = [ece_keys[j] for j in range(len(valid_mask)) if valid_mask[j]]
                    self._register_click_points(self.ax1, valid_R, valid_Te, valid_ece_keys)

            except Exception as e:
                print(f"[ne/Te] Error plotting {entry}: {str(e)}")

        # Set limits
        self.ax1.set_xlim(self.app_config.R_LIMITS)
        te_margin = te_max * 0.1
        self.ax1.set_ylim(0, te_max + te_margin)

        if has_thomson:
            self.ax2.set_xlim(self.app_config.R_LIMITS)
            ne_margin = ne_max * 0.1
            self.ax2.set_ylim(0, ne_max + ne_margin)

        # Overlay fit curves if available
        self._overlay_fit_curves(self.ax1, 'Te', selected_entries, colors)
        self._overlay_fit_curves(self.ax2, 'ne', selected_entries, colors)

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
        self._clear_click_points()

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

        te_max, ne_max = 0, 0
        colors = self._get_plot_colors(selected_entries)

        has_thomson = False

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

                if source in ['TS', 'TS+ECE']:
                    has_thomson = True
                    cache_key = f'{shot_number}_TS'
                    data = self.data[cache_key]

                    ts_rshift = self._get_rshift('Thomson')
                    x_data = interp_func(data.radius + ts_rshift)

                    Te_data, Te_err_upper, Te_err_lower = data.get_parameter_asymmetric('Te')
                    ne_data, ne_err_upper, ne_err_lower = data.get_parameter_asymmetric('ne')

                    # Markers match the data that was actually fit:
                    # if a fit exists for this entry on the current axis, use
                    # the same dt-averaged data; otherwise use single slice.
                    has_fit = (entry in self.fit_results
                               and self.fit_x_axis == x_axis
                               and self.fit_results[entry])
                    eff_dt_s = self.fit_dt_s if has_fit else 0.0
                    if eff_dt_s > 0:
                        ts_t_mask = (data.time >= time_point - eff_dt_s) & (data.time <= time_point + eff_dt_s)
                        n_avg = int(np.sum(ts_t_mask))
                    else:
                        n_avg = 0

                    if n_avg > 0:
                        Te_profile = np.nanmean(Te_data[:, ts_t_mask], axis=1)
                        Te_err_upper_profile = np.sqrt(np.nanmean(Te_err_upper[:, ts_t_mask] ** 2, axis=1) / n_avg)
                        Te_err_lower_profile = np.sqrt(np.nanmean(Te_err_lower[:, ts_t_mask] ** 2, axis=1) / n_avg)
                        ne_profile = np.nanmean(ne_data[:, ts_t_mask], axis=1)
                        ne_err_upper_profile = np.sqrt(np.nanmean(ne_err_upper[:, ts_t_mask] ** 2, axis=1) / n_avg)
                        ne_err_lower_profile = np.sqrt(np.nanmean(ne_err_lower[:, ts_t_mask] ** 2, axis=1) / n_avg)
                    else:
                        time_idx = np.argmin(np.abs(data.time - time_point))
                        Te_profile = Te_data[:, time_idx]
                        Te_err_upper_profile = Te_err_upper[:, time_idx]
                        Te_err_lower_profile = Te_err_lower[:, time_idx]
                        ne_profile = ne_data[:, time_idx]
                        ne_err_upper_profile = ne_err_upper[:, time_idx]
                        ne_err_lower_profile = ne_err_lower[:, time_idx]

                    # Split enabled/disabled Thomson channels
                    ts_keys = [f"TS_{j}" for j in range(len(x_data))]
                    ts_mask = self._get_channel_mask(ts_keys)

                    lcfs_idx = np.argmin(np.abs(x_data - 1))
                    vm = ts_mask[:lcfs_idx]
                    if vm.any():
                        te_max = max(te_max, np.nanmax(Te_profile[:lcfs_idx][vm]))
                        ne_max = max(ne_max, np.nanmax(ne_profile[:lcfs_idx][vm]))

                    time_ms = entry.split("_")[1].split()[0]
                    if n_avg > 0:
                        eff_dt_ms = eff_dt_s * 1000.0
                        label = f'#{shot_number} {time_ms}ms±{eff_dt_ms:.0f}ms (TS)'
                    else:
                        label = f'#{shot_number} {time_ms}ms (TS)'

                    mk_kwargs = {}
                    if n_avg > 0:
                        mk_kwargs = {'markeredgecolor': 'black', 'markeredgewidth': 1.0}

                    if ts_mask.any():
                        self.ax1.errorbar(x_data[ts_mask], Te_profile[ts_mask],
                                         yerr=[Te_err_lower_profile[ts_mask], Te_err_upper_profile[ts_mask]],
                                         fmt='o', capsize=5, label=label,
                                         color=color, markersize=5, zorder=10, **mk_kwargs)
                        self.ax2.errorbar(x_data[ts_mask], ne_profile[ts_mask],
                                         yerr=[ne_err_lower_profile[ts_mask], ne_err_upper_profile[ts_mask]],
                                         fmt='o', capsize=5, label=label,
                                         color=color, markersize=5, zorder=10, **mk_kwargs)
                    if (~ts_mask).any():
                        self.ax1.errorbar(x_data[~ts_mask], Te_profile[~ts_mask],
                                         yerr=[Te_err_lower_profile[~ts_mask], Te_err_upper_profile[~ts_mask]],
                                         fmt='o', capsize=5, label='',
                                         color=(0.6, 0.6, 0.6, 0.35), markersize=5, zorder=5)
                        self.ax2.errorbar(x_data[~ts_mask], ne_profile[~ts_mask],
                                         yerr=[ne_err_lower_profile[~ts_mask], ne_err_upper_profile[~ts_mask]],
                                         fmt='o', capsize=5, label='',
                                         color=(0.6, 0.6, 0.6, 0.35), markersize=5, zorder=5)

                    # Add Thomson channel labels
                    n_channels = len(x_data)
                    for ch_idx in range(n_channels):
                        if self.show_channel_checkbox.isChecked():
                            if ch_idx < 14:
                                lbl = f'TS_CORE{ch_idx + 1}'
                            else:
                                lbl = f'TS_EDGE{ch_idx - 14 + 1}'
                            self.ax1.annotate(lbl, (x_data[ch_idx], Te_profile[ch_idx]),
                                             textcoords='offset points',
                                             xytext=(0, 5), ha='center', fontsize=7, alpha=0.8,
                                             clip_on=True, annotation_clip=True)
                            self.ax2.annotate(lbl, (x_data[ch_idx], ne_profile[ch_idx]),
                                             textcoords='offset points',
                                             xytext=(0, 5), ha='center', fontsize=7, alpha=0.8,
                                             clip_on=True, annotation_clip=True)

                    # Register for double-click toggle
                    self._register_click_points(self.ax1, x_data, Te_profile, ts_keys)
                    self._register_click_points(self.ax2, x_data, ne_profile, ts_keys)

                    # Plot ECE data only for "TS+ECE" mode
                    if source == 'TS+ECE':
                        ece_cache_key = self._get_ece_cache_key(shot_number)
                        if ece_cache_key is not None:
                            ece_data = self.ece_data_cache[ece_cache_key]
                            ece_rshift = self._get_rshift('ECE')
                            ece_x_data = interp_func(ece_data.radius + ece_rshift)

                            ece_Te_data = ece_data.measurements['Te']['data']

                            # Match TS branch: average if a fit was done with dt
                            if n_avg > 0:
                                ece_t_mask = (ece_data.time >= time_point - eff_dt_s) & (ece_data.time <= time_point + eff_dt_s)
                                ece_n_avg = int(np.sum(ece_t_mask))
                            else:
                                ece_n_avg = 0
                            if ece_n_avg > 0:
                                ece_Te_profile = np.nanmean(ece_Te_data[:, ece_t_mask], axis=1)
                            else:
                                ece_Te_profile = ece_Te_data[:, np.argmin(np.abs(ece_data.time - time_point))]

                            valid_mask = ece_data.measurements['Te']['valid_mask']

                            # Split enabled/disabled ECE channels
                            ece_keys = self._ece_channel_keys(ece_data)
                            ece_ch_mask = self._get_channel_mask(ece_keys)

                            in_range = (ece_x_data <= 1.0)
                            if np.any(valid_mask & ece_ch_mask & in_range):
                                te_max = max(te_max, np.nanmax(ece_Te_profile[valid_mask & ece_ch_mask & in_range]))

                            ece_time_ms = entry.split("_")[1].split()[0]
                            if ece_n_avg > 0:
                                ece_label = f'#{shot_number} {ece_time_ms}ms±{eff_dt_ms:.0f}ms (ECE)'
                            else:
                                ece_label = f'#{shot_number} {ece_time_ms}ms (ECE)'
                            enabled = valid_mask & ece_ch_mask
                            disabled = valid_mask & (~ece_ch_mask)

                            if enabled.any():
                                self.ax1.plot(ece_x_data[enabled], ece_Te_profile[enabled],
                                             's', color=color, markersize=5,
                                             markerfacecolor='none', markeredgewidth=1.5,
                                             label=ece_label, zorder=5)
                            if disabled.any():
                                self.ax1.plot(ece_x_data[disabled], ece_Te_profile[disabled],
                                             's', color=(0.6, 0.6, 0.6, 0.35), markersize=5,
                                             markerfacecolor='none', markeredgewidth=1.5,
                                             label='', zorder=3)

                            # Add ECE channel labels
                            ece_channels = ece_data.measurements['Te'].get('channels',
                                           list(range(1, len(ece_data.radius) + 1)))
                            valid_x = ece_x_data[valid_mask]
                            valid_Te = ece_Te_profile[valid_mask]
                            valid_ch = [ece_channels[j] for j in range(len(valid_mask)) if valid_mask[j]]
                            self._add_channel_labels(self.ax1, valid_x, valid_Te, 'ECE', valid_ch)

                            # Register ECE for double-click toggle
                            valid_ece_keys = [ece_keys[j] for j in range(len(valid_mask)) if valid_mask[j]]
                            self._register_click_points(self.ax1, valid_x, valid_Te, valid_ece_keys)

                elif source == 'ECE':
                    cache_key = self._get_ece_cache_key(shot_number)

                    if cache_key:
                        ece_data = self.ece_data_cache[cache_key]
                        ece_rshift = self._get_rshift('ECE')
                        x_data = interp_func(ece_data.radius + ece_rshift)

                        Te_data = ece_data.measurements['Te']['data']

                        # Average if a fit was done with dt
                        has_fit = (entry in self.fit_results
                                   and self.fit_x_axis == x_axis
                                   and self.fit_results[entry])
                        eff_dt_s = self.fit_dt_s if has_fit else 0.0
                        if eff_dt_s > 0:
                            t_mask = (ece_data.time >= time_point - eff_dt_s) & (ece_data.time <= time_point + eff_dt_s)
                            n_avg = int(np.sum(t_mask))
                        else:
                            n_avg = 0
                        if n_avg > 0:
                            Te_profile = np.nanmean(Te_data[:, t_mask], axis=1)
                        else:
                            time_idx = np.argmin(np.abs(ece_data.time - time_point))
                            Te_profile = Te_data[:, time_idx]

                        valid_mask = ece_data.measurements['Te']['valid_mask']

                        # Split enabled/disabled ECE channels
                        ece_keys = self._ece_channel_keys(ece_data)
                        ece_ch_mask = self._get_channel_mask(ece_keys)

                        in_range = (x_data <= 1.0)
                        if np.any(valid_mask & ece_ch_mask & in_range):
                            te_max = max(te_max, np.nanmax(Te_profile[valid_mask & ece_ch_mask & in_range]))

                        time_ms = entry.split("_")[1].split()[0]
                        if n_avg > 0:
                            eff_dt_ms = eff_dt_s * 1000.0
                            label = f'#{shot_number} {time_ms}ms±{eff_dt_ms:.0f}ms (ECE)'
                        else:
                            label = f'#{shot_number} {time_ms}ms (ECE)'
                        enabled = valid_mask & ece_ch_mask
                        disabled = valid_mask & (~ece_ch_mask)

                        if enabled.any():
                            self.ax1.plot(x_data[enabled], Te_profile[enabled],
                                         's', color=color, markersize=5,
                                         markerfacecolor='none', markeredgewidth=1.5,
                                         label=label, zorder=5)
                        if disabled.any():
                            self.ax1.plot(x_data[disabled], Te_profile[disabled],
                                         's', color=(0.6, 0.6, 0.6, 0.35), markersize=5,
                                         markerfacecolor='none', markeredgewidth=1.5,
                                         label='', zorder=3)

                        # Add ECE channel labels
                        ece_channels = ece_data.measurements['Te'].get('channels',
                                       list(range(1, len(ece_data.radius) + 1)))
                        valid_x = x_data[valid_mask]
                        valid_Te = Te_profile[valid_mask]
                        valid_ch = [ece_channels[j] for j in range(len(valid_mask)) if valid_mask[j]]
                        self._add_channel_labels(self.ax1, valid_x, valid_Te, 'ECE', valid_ch)

                        # Register ECE for double-click toggle
                        valid_ece_keys = [ece_keys[j] for j in range(len(valid_mask)) if valid_mask[j]]
                        self._register_click_points(self.ax1, valid_x, valid_Te, valid_ece_keys)

            except Exception as e:
                print(f"[ne/Te] Error plotting {entry}: {str(e)}")

        # Overlay fit curves
        self._overlay_fit_curves(self.ax1, 'Te', selected_entries, colors)
        self._overlay_fit_curves(self.ax2, 'ne', selected_entries, colors)

        self.ax1.set_xlabel(x_label)
        self.ax2.set_xlabel(x_label)
        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_ylabel(self.param2['label'])

        x_left, x_right = -0.1, 1.1
        try:
            fit_left = float(self.fit_xmin_entry.text()) - 0.1
            fit_right = float(self.fit_xmax_entry.text()) + 0.1
            if fit_left < x_left:
                x_left = fit_left
            if fit_right > x_right:
                x_right = fit_right
        except (ValueError, AttributeError):
            pass
        self.ax1.set_xlim(x_left, x_right)
        self.ax2.set_xlim(x_left, x_right)
        self.ax1.set_ylim(0, te_max * 1.1)

        if has_thomson:
            self.ax2.set_ylim(0, ne_max * 1.1)

        for ax in [self.ax1, self.ax2]:
            ax.axvspan(-1, 0, color='gray', alpha=0.15, zorder=0)
            ax.axvspan(1, 2, color='gray', alpha=0.15, zorder=0)

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
        """Write ne, Te profile data to text file"""
        dt_s = self._get_dt_seconds()
        dt_ms = dt_s * 1000.0

        with open(file_path, 'w') as f:
            f.write("# ne, Te Profile Data\n")
            if self.computed_efit_tree:
                f.write(f"# EFIT Source: {self.computed_efit_tree}\n")
            if dt_s > 0:
                f.write(f"# Time averaging: dt = {dt_ms:.0f} ms\n")
            f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "psi_N", "rho_pol", "rho_tor",
                "Te[keV]", "Te_err_u", "Te_err_l", "ne[e19]", "ne_err_u", "ne_err_l", "Source"
            ))

            for entry in selected_entries:
                shot_number, time_point, source = self._parse_entry(entry)
                actual_time = time_point

                if source in ['TS', 'TS+ECE']:
                    cache_key = f'{shot_number}_TS'
                    data = self.data.get(cache_key)
                    if data is None:
                        data = self.data_loader.load_data(shot_number)
                        self.data[cache_key] = data

                    R_data = data.radius + self._get_rshift('Thomson')

                    Te_data, Te_err_upper, Te_err_lower = data.get_parameter_asymmetric('Te')
                    ne_data, ne_err_upper, ne_err_lower = data.get_parameter_asymmetric('ne')

                    if dt_s > 0:
                        t_mask = (data.time >= time_point - dt_s) & (data.time <= time_point + dt_s)
                        n_avg = np.sum(t_mask)
                        if n_avg > 0:
                            Te_profile = np.nanmean(Te_data[:, t_mask], axis=1)
                            Te_err_upper_profile = np.sqrt(np.nanmean(Te_err_upper[:, t_mask] ** 2, axis=1) / n_avg)
                            Te_err_lower_profile = np.sqrt(np.nanmean(Te_err_lower[:, t_mask] ** 2, axis=1) / n_avg)
                            ne_profile = np.nanmean(ne_data[:, t_mask], axis=1)
                            ne_err_upper_profile = np.sqrt(np.nanmean(ne_err_upper[:, t_mask] ** 2, axis=1) / n_avg)
                            ne_err_lower_profile = np.sqrt(np.nanmean(ne_err_lower[:, t_mask] ** 2, axis=1) / n_avg)
                        else:
                            time_idx = np.argmin(np.abs(data.time - time_point))
                            Te_profile = Te_data[:, time_idx]
                            Te_err_upper_profile = Te_err_upper[:, time_idx]
                            Te_err_lower_profile = Te_err_lower[:, time_idx]
                            ne_profile = ne_data[:, time_idx]
                            ne_err_upper_profile = ne_err_upper[:, time_idx]
                            ne_err_lower_profile = ne_err_lower[:, time_idx]
                    else:
                        time_idx = np.argmin(np.abs(data.time - time_point))
                        Te_profile = Te_data[:, time_idx]
                        Te_err_upper_profile = Te_err_upper[:, time_idx]
                        Te_err_lower_profile = Te_err_lower[:, time_idx]
                        ne_profile = ne_data[:, time_idx]
                        ne_err_upper_profile = ne_err_upper[:, time_idx]
                        ne_err_lower_profile = ne_err_lower[:, time_idx]

                    # Get EFIT values
                    psi_n, rho_pol, rho_tor = self._get_efit_values_at_R(entry, R_data)

                    # Write Thomson data rows (10 char width each)
                    for i in range(len(R_data)):
                        f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10s\n" % (
                            shot_number, actual_time, R_data[i],
                            self._format_value(psi_n[i]),
                            self._format_value(rho_pol[i]),
                            self._format_value(rho_tor[i]),
                            Te_profile[i], Te_err_upper_profile[i], Te_err_lower_profile[i],
                            ne_profile[i], ne_err_upper_profile[i], ne_err_lower_profile[i], 'TS'
                        ))

                    # Check if ECE data exists for "TS+ECE" mode
                    ece_cache_key = self._get_ece_cache_key(shot_number)
                    if ece_cache_key is not None:
                        ece_data = self.ece_data_cache[ece_cache_key]
                        ece_R_data = ece_data.radius
                        ece_Te_data = ece_data.measurements['Te']['data']

                        if dt_s > 0:
                            ece_t_mask = (ece_data.time >= time_point - dt_s) & (ece_data.time <= time_point + dt_s)
                            if np.sum(ece_t_mask) > 0:
                                ece_Te_profile = np.nanmean(ece_Te_data[:, ece_t_mask], axis=1)
                            else:
                                ece_Te_profile = ece_Te_data[:, np.argmin(np.abs(ece_data.time - time_point))]
                        else:
                            ece_Te_profile = ece_Te_data[:, np.argmin(np.abs(ece_data.time - time_point))]
                        valid_mask = ece_data.measurements['Te']['valid_mask']

                        # Get EFIT values for ECE
                        ece_psi_n, ece_rho_pol, ece_rho_tor = self._get_efit_values_at_R(entry, ece_R_data)

                        # Write ECE data rows (valid channels only, 10 char width each)
                        for i in range(len(ece_R_data)):
                            if valid_mask[i]:
                                f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                                    shot_number, actual_time, ece_R_data[i],
                                    self._format_value(ece_psi_n[i]),
                                    self._format_value(ece_rho_pol[i]),
                                    self._format_value(ece_rho_tor[i]),
                                    ece_Te_profile[i], 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'ECE'
                                ))

                elif source == 'ECE':
                    cache_key = self._get_ece_cache_key(shot_number)
                    if cache_key is None:
                        loader = self._get_ece_loader()
                        ece_data = loader.load_data(shot_number, sampling_rate=0.01)
                        cache_key = f'{shot_number}_ECE_100Hz'
                        self.ece_data_cache[cache_key] = ece_data
                    else:
                        ece_data = self.ece_data_cache[cache_key]

                    time_idx = np.argmin(np.abs(ece_data.time - time_point))
                    actual_time = ece_data.time[time_idx]
                    R_data = ece_data.radius

                    Te_data = ece_data.measurements['Te']['data']
                    Te_profile = Te_data[:, time_idx]
                    valid_mask = ece_data.measurements['Te']['valid_mask']

                    # Get EFIT values
                    psi_n, rho_pol, rho_tor = self._get_efit_values_at_R(entry, R_data)

                    # Write ECE data rows (valid channels only, 10 char width each)
                    for i in range(len(R_data)):
                        if valid_mask[i]:
                            f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                                shot_number, actual_time, R_data[i],
                                self._format_value(psi_n[i]),
                                self._format_value(rho_pol[i]),
                                self._format_value(rho_tor[i]),
                                Te_profile[i], 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'ECE'
                            ))

    # ===== TCI Validation =====

    def _add_extra_fitting_controls(self, grid, row):
        """Add TCI validation toggle switch to fitting section"""
        from ui.widgets.toggle_switch import ToggleSwitch
        from config.user_settings import get_tab_settings

        tci_row = QHBoxLayout()

        self.tci_validate_toggle = ToggleSwitch()
        self.tci_validate_toggle.setToolTip(
            "After fitting, compare synthetic line-averaged ne\n"
            "from fitted profile with measured TCI01~05 values")

        # Restore from settings
        saved = get_tab_settings(self._settings_key).get("tci_validation", False)
        self.tci_validate_toggle.setChecked(saved, animate=False)

        tci_row.addWidget(self.tci_validate_toggle)
        tci_row.addWidget(QLabel("TCI Validation"))
        tci_row.addStretch()

        grid.addLayout(tci_row, row, 0, 1, 4)
        return row + 1

    def _get_tci_loader(self):
        """Lazy initialization of TCI loader"""
        if self.tci_loader is None:
            from data_loaders.tci_loader import TCILoader
            from config.diagnostic_config import DIAGNOSTICS
            self.tci_loader = TCILoader(self.app_config, DIAGNOSTICS.get('TCI', {}))
        return self.tci_loader

    def _post_fitting(self, selected_entries, x_axis):
        """Run TCI validation after fitting if checkbox is checked"""
        if not hasattr(self, 'tci_validate_toggle') or not self.tci_validate_toggle.isChecked():
            return
        if 'ne' not in self.fit_results.get(selected_entries[0], {}):
            print("[TCI] Skipped: no ne fit results available")
            return

        self._validate_tci(selected_entries, x_axis)

    def _validate_tci(self, selected_entries, x_axis):
        """Compare fitted ne profile with TCI line-integrated density"""
        print("\n" + "=" * 60)
        print("[TCI] Line-Integrated Density Validation")
        print("-" * 60)

        for entry in selected_entries:
            if entry not in self.fit_results or 'ne' not in self.fit_results[entry]:
                continue

            ne_result = self.fit_results[entry]['ne']
            if not ne_result.success:
                print(f"[TCI] {entry}: ne fit failed, skipping")
                continue

            shot_number, time_point, source = self._parse_entry(entry)

            if entry not in self.efit_data:
                print(f"[TCI] {entry}: no EFIT data, skipping")
                continue

            efit_entry = self.efit_data[entry]
            efit_R = efit_entry['R']

            # Build interpolation: R -> flux coordinate at midplane
            if x_axis == "psi_N":
                coord_at_R = interp1d(efit_R, efit_entry['psi_N'],
                                      fill_value='extrapolate', bounds_error=False)
            elif x_axis == "rho_pol":
                coord_at_R = interp1d(efit_R, efit_entry['rho_pol'],
                                      fill_value='extrapolate', bounds_error=False)
            else:
                coord_at_R = interp1d(efit_R, efit_entry['rho_tor'],
                                      fill_value='extrapolate', bounds_error=False)

            # Build ne(flux_coord) from fit result
            ne_fit_interp = interp1d(ne_result.x_fit, ne_result.y_fit,
                                     fill_value=0.0, bounds_error=False)

            # Load TCI data for this shot
            try:
                tci_loader = self._get_tci_loader()
                tci_data = tci_loader.load_data(shot_number, quiet=True)
            except Exception as e:
                print(f"[TCI] Failed to load TCI data for #{shot_number}: {e}")
                continue

            # Find nearest TCI time index
            tci_time_idx = np.argmin(np.abs(tci_data.time - time_point))
            if np.abs(tci_data.time[tci_time_idx] - time_point) > 0.05:
                print(f"[TCI] {entry}: no TCI data within ±50ms")
                continue

            ne_tci_measured = tci_data.measurements['ne']['data']  # (5, n_time)

            print(f"[TCI] #{shot_number} {time_point*1e3:.0f}ms ({x_axis}, {ne_result.func_name})")
            print(f"[TCI]   {'Channel':>8s}  {'<ne>_fit':>12s}  {'<ne>_TCI':>12s}  {'fit/TCI ratio':>14s}")
            print(f"[TCI]   {'─'*8}  {'─'*12}  {'─'*12}  {'─'*14}")

            for ch_idx, (ch_name, ch_info) in enumerate(self.TCI_CHANNELS.items()):
                rmid = ch_info['Rmid']

                # Compute synthetic line-averaged ne along tangential chord
                synth_ne_avg = self._compute_line_averaged_ne(
                    rmid, coord_at_R, ne_fit_interp)

                # Measured TCI value: \ne_tci is line-averaged density [1e19 m^-3]
                meas_ne_avg = ne_tci_measured[ch_idx, tci_time_idx]

                if meas_ne_avg != 0 and not np.isnan(meas_ne_avg):
                    ratio = synth_ne_avg / meas_ne_avg
                    print(f"[TCI]   {ch_name:>8s}  {synth_ne_avg:12.4f}  {meas_ne_avg:12.4f}  {ratio:14.3f}")
                else:
                    print(f"[TCI]   {ch_name:>8s}  {synth_ne_avg:12.4f}  {'N/A':>12s}  {'N/A':>14s}")

        print("=" * 60 + "\n")

    def _compute_line_averaged_ne(self, rmid, coord_at_R, ne_fit_interp, n_points=500):
        """Compute line-averaged ne along a tangential TCI chord.

        TCI geometry: horizontal tangential beam with closest approach at R=rmid.
        Along the beam path at distance s from tangent point: R(s) = sqrt(rmid^2 + s^2)
        Integration is symmetric, so we compute one side and multiply by 2.
        Line-averaged: <ne> = ∫ne·ds / L_chord (round-trip factor cancels out).

        Args:
            rmid: Tangent radius [m]
            coord_at_R: Interpolation function R -> flux coordinate
            ne_fit_interp: Interpolation function flux_coord -> ne [1e19 m^-3]
            n_points: Number of integration points per half-chord

        Returns:
            Line-averaged ne [1e19 m^-3]
        """
        rend = self.TCI_REND
        if rmid >= rend:
            return 0.0

        # Half-chord length
        L_half = np.sqrt(rend**2 - rmid**2)

        # Path parameter s: distance from tangent point
        s = np.linspace(0, L_half, n_points)

        # R along the chord
        R_along = np.sqrt(rmid**2 + s**2)

        # Map R to flux coordinate, then evaluate ne
        # Use abs() because EFIT convention makes HFS psi_N negative,
        # but ne profile is symmetric about the magnetic axis
        flux_along = np.abs(coord_at_R(R_along))
        ne_along = ne_fit_interp(flux_along)

        # Clip negative values
        ne_along = np.maximum(ne_along, 0.0)

        # Line-averaged = ∫ne·ds / L_half (symmetric, 2× cancels)
        ne_line_avg = np.trapz(ne_along, x=s) / L_half

        return ne_line_avg
