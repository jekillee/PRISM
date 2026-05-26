"""
Electron (ne, Te) Time Trace tab with unified Thomson/ECE/TCI loading
"""

import numpy as np

from PySide6.QtWidgets import (
    QMessageBox, QLineEdit, QComboBox, QPushButton, QWidget,
    QHBoxLayout, QVBoxLayout, QGroupBox, QGridLayout, QLabel,
    QFileDialog, QStyle
)

from ui.tabs.timetrace_base_tab import TimeTraceBaseTab
from ui.ui_constants import apply_shot_arrow_icons


class ElectronTimeTraceTab(TimeTraceBaseTab):
    """Electron (ne, Te) Time Trace tab with unified Thomson/ECE/TCI loading"""

    TAB_NAME = "Electron"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ece_loader = None
        self.ece_data_cache = {}
        self.tci_loader = None
        self.tci_data_cache = {}
        self.ts_ne_avg_cache = {}  # {shot: {'time': ndarray, 'data': ndarray}}

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
        btn_updown_layout = QVBoxLayout(btn_updown)
        btn_updown_layout.setContentsMargins(0, 0, 0, 0)
        btn_updown_layout.setSpacing(0)
        mini_btn_style = "padding: 0px; border-radius: 2px;"
        up_btn = QPushButton()
        up_btn.setFixedSize(24, 15)
        up_btn.setStyleSheet(mini_btn_style)
        up_btn.clicked.connect(lambda: self._adjust_shot(1))
        btn_updown_layout.addWidget(up_btn)
        down_btn = QPushButton()
        down_btn.setFixedSize(24, 15)
        down_btn.setStyleSheet(mini_btn_style)
        down_btn.clicked.connect(lambda: self._adjust_shot(-1))
        btn_updown_layout.addWidget(down_btn)
        apply_shot_arrow_icons(up_btn, down_btn)
        grid.addWidget(btn_updown, 0, 2)

        diag_options = ['TS', 'ECE (100Hz)', 'ECE (1kHz)', 'TCI (100Hz)', 'TCI (1kHz)']
        self.diag_combo = QComboBox()
        self.diag_combo.addItems(diag_options)
        self.diag_combo.setCurrentText('TS')
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
            self.shot_entry.clear()
            self.shot_entry.setText(str(new_shot))
        except ValueError:
            pass

    def _get_ece_loader(self):
        """On-demand initialization of ECE loader"""
        if self.ece_loader is None:
            from data_loaders.ece_loader import ECELoader
            from config.diagnostic_config import DIAGNOSTICS
            self.ece_loader = ECELoader(self.app_config, DIAGNOSTICS['ECE'])
        return self.ece_loader

    def _get_tci_loader(self):
        """On-demand initialization of TCI loader"""
        if self.tci_loader is None:
            from data_loaders.tci_loader import TCILoader
            from config.diagnostic_config import DIAGNOSTICS
            self.tci_loader = TCILoader(self.app_config, DIAGNOSTICS['TCI'])
        return self.tci_loader

    def load_shot_data(self):
        """Load Thomson/ECE/TCI data based on selection"""
        try:
            shot_number = int(self.shot_entry.text())
            selection = self.diag_combo.currentText()

            self._set_status(f"Loading #{shot_number} ({selection})...", 'blue')
            all_channels = []

            # Parse selection to determine mode
            load_thomson = False
            load_ece = False
            load_tci = False
            sampling_rate = None
            sampling_key = None

            if selection.startswith('ECE'):
                load_ece = True
                if '100Hz' in selection:
                    sampling_rate = 0.01
                    sampling_key = '100Hz'
                else:  # 1kHz
                    sampling_rate = 0.001
                    sampling_key = '1kHz'
            elif selection == 'TS':
                load_thomson = True
            elif selection.startswith('TCI'):
                load_tci = True
                if '100Hz' in selection:
                    sampling_rate = 0.01
                    sampling_key = '100Hz'
                else:  # 1kHz
                    sampling_rate = 0.001
                    sampling_key = '1kHz'

            if load_thomson:
                # Load Thomson data
                try:
                    ts_data = self.data_loader.load_data(shot_number)
                    cache_key = f'{shot_number}_TS'
                    self.data[cache_key] = ts_data

                    # Thomson Core: TSC01-14, Edge: TSE01-17
                    for i, radius in enumerate(ts_data.radius):
                        if i < 14:
                            ch_label = f'TSC{i+1:02d}'
                        else:
                            ch_label = f'TSE{i-14+1:02d}'

                        all_channels.append({
                            'R': radius,
                            'label': f'{shot_number:06d}_{radius*1e3:.0f} ({ch_label})',
                            'source': 'TS',
                            'channel_idx': i,
                            'ch_label': ch_label
                        })

                    print(f"[Thomson] Data loaded: {len(ts_data.radius)} channels")

                    # Add TS_NE_AVG if available from loader
                    if 'ne_avg' in ts_data.measurements:
                        ne_avg = ts_data.measurements['ne_avg']
                        self.ts_ne_avg_cache[shot_number] = {
                            'time': ne_avg['time'],
                            'data': ne_avg['data']
                        }
                        all_channels.append({
                            'R': 0,
                            'label': f'{shot_number:06d} (TS_NE_AVG)',
                            'source': 'TS_NE_AVG',
                            'channel_idx': 0,
                            'ch_label': 'TS_NE_AVG'
                        })

                except Exception as e:
                    print(f"[Thomson] Not available: {str(e)}")

            if load_ece:
                # Load ECE data
                try:
                    loader = self._get_ece_loader()
                    ece_data = loader.load_data(shot_number, sampling_rate=sampling_rate)

                    cache_key = f'{shot_number}_ECE_{sampling_key}'
                    self.ece_data_cache[cache_key] = ece_data

                    channels = ece_data.measurements['Te'].get('channels',
                               list(range(1, len(ece_data.radius) + 1)))

                    for i, radius in enumerate(ece_data.radius):
                        ch_num = channels[i] if i < len(channels) else i + 1
                        ch_label = f'ECE{ch_num:02d}'
                        all_channels.append({
                            'R': radius,
                            'label': f'{shot_number:06d}_{radius*1e3:.0f} ({ch_label}) [{sampling_key}]',
                            'source': 'ECE',
                            'channel_idx': i,
                            'ch_label': ch_label,
                            'sampling_key': sampling_key
                        })

                    n_valid = np.sum(ece_data.measurements['Te']['valid_mask'])
                    n_overlap = np.sum(ece_data.measurements['Te']['overlap_mask'])
                    print(f"[ECE] Data loaded: {len(ece_data.radius)} channels @ {sampling_key}")
                    print(f"[ECE]   Valid: {n_valid}, Overlap: {n_overlap}")
                except Exception as e:
                    print(f"[ECE] Not available: {str(e)}")

            if load_tci:
                # Load TCI data with resampling
                try:
                    loader = self._get_tci_loader()
                    tci_data = loader.load_data(shot_number, sampling_rate=sampling_rate)
                    cache_key = f'{shot_number}_TCI_{sampling_key}'
                    self.tci_data_cache[cache_key] = tci_data

                    channels = tci_data.measurements['ne']['channels']
                    for i, ch in enumerate(channels):
                        ch_label = f'TCI{ch:02d}'
                        all_channels.append({
                            'R': ch,  # Use channel number for sorting
                            'label': f'{shot_number:06d} ({ch_label}) [{sampling_key}]',
                            'source': 'TCI',
                            'channel_idx': i,
                            'ch_label': ch_label,
                            'sampling_key': sampling_key
                        })

                    print(f"[TCI] Data loaded: {len(channels)} channels @ {sampling_key}")
                except Exception as e:
                    print(f"[TCI] Not available: {str(e)}")

            if not all_channels:
                self._set_status(f"No data for #{shot_number}", 'red')
                QMessageBox.critical(self.frame, "Error", "No data available for this shot")
                return

            # Sort by R position (or channel number for TCI)
            all_channels.sort(key=lambda x: x['R'])

            # Update available listbox
            self.available_listbox.clear()
            for channel in all_channels:
                self.available_listbox.addItem(channel['label'])

            print(f"[ne/Te] Total channels loaded: {len(all_channels)}")

            if hasattr(self, 'browse_button'):
                self.browse_button.setEnabled(True)
                self.browse_button.setText(f"Browse #{shot_number}")
                self._last_shot = shot_number

            self._set_status(f"#{shot_number} loaded ({selection}): {len(all_channels)} channels", 'green')

        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
        except Exception as e:
            self._set_status(f"Failed to load: {e}", 'red')
            QMessageBox.critical(self.frame, "Error", f"Failed to load data: {str(e)}")

    def _parse_entry(self, entry):
        """Parse time trace entry with sampling rate info

        Format: XXXXXX_YYYY (DIAG##) [rate] or XXXXXX (TCI##) [rate]
        Returns: shot_number, radius_or_channel, source, channel_label, sampling_key
        """
        # Extract sampling key if present
        sampling_key = None
        if '[' in entry and ']' in entry:
            sampling_key = entry.split('[')[1].split(']')[0].strip()
            entry_base = entry.split('[')[0].strip()
        else:
            entry_base = entry

        # Format: 012345_1800 (ECE01) or 012345 (TCI01)
        main_part = entry_base.split('(')[0].strip()
        diag_label = entry_base.split('(')[1].split(')')[0].strip()

        parts = main_part.split('_')
        shot_number = int(parts[0])

        # Determine source from diag_label
        if diag_label.startswith('TSC') or diag_label.startswith('TSE'):
            source = 'TS'
            radius = float(parts[1]) / 1e3  # mm to m
            return shot_number, radius, source, diag_label, sampling_key
        elif diag_label.startswith('ECE'):
            source = 'ECE'
            radius = float(parts[1]) / 1e3  # mm to m
            return shot_number, radius, source, diag_label, sampling_key
        elif diag_label == 'TS_NE_AVG':
            source = 'TS_NE_AVG'
            return shot_number, 0, source, diag_label, sampling_key
        elif diag_label.startswith('TCI'):
            source = 'TCI'
            channel = int(diag_label.replace('TCI', ''))
            return shot_number, channel, source, diag_label, sampling_key
        else:
            # Fallback
            radius = float(parts[1]) / 1e3 if len(parts) > 1 and parts[1] else 0
            return shot_number, radius, 'TS', diag_label, sampling_key

    def _get_timetrace_preview_channels(self):
        shot = getattr(self, '_last_shot', None)
        if shot is None:
            return []

        channels = []

        # Thomson channels
        ts_key = f'{shot}_TS'
        ts_data = self.data.get(ts_key)
        if ts_data is not None:
            Te_data, _ = ts_data.get_parameter('Te')
            ne_data, _ = ts_data.get_parameter('ne')
            for i, R in enumerate(ts_data.radius):
                ch_label = f'TSC{i+1:02d}' if i < 14 else f'TSE{i-14+1:02d}'
                channels.append({
                    'entry': f'{shot:06d}_{R*1e3:.0f} ({ch_label})',
                    'label': f'{ch_label} R={R*1e3:.0f}mm',
                    'time': ts_data.time,
                    'param1': Te_data[i, :],
                    'param2': ne_data[i, :],
                    'R': float(R),
                    'shot': shot,
                    'source': 'TS',
                })

        # ECE channels (across loaded sampling rates)
        for cache_key, ece_data in self.ece_data_cache.items():
            if not cache_key.startswith(f'{shot}_ECE_'):
                continue
            sampling_key = cache_key.split('_')[-1]
            ece_channels = ece_data.measurements['Te'].get(
                'channels', list(range(1, len(ece_data.radius) + 1)))
            Te_2d = ece_data.measurements['Te']['data']
            for i, R in enumerate(ece_data.radius):
                ch_num = ece_channels[i] if i < len(ece_channels) else i + 1
                ch_label = f'ECE{ch_num:02d}'
                channels.append({
                    'entry': f'{shot:06d}_{R*1e3:.0f} ({ch_label}) [{sampling_key}]',
                    'label': f'{ch_label} R={R*1e3:.0f}mm [{sampling_key}]',
                    'time': ece_data.time,
                    'param1': Te_2d[i, :],
                    'param2': None,
                    'R': float(R),
                    'shot': shot,
                    'source': 'ECE',
                })

        # TCI channels (no R; use channel number for ordering)
        for cache_key, tci_data in self.tci_data_cache.items():
            if not cache_key.startswith(f'{shot}_TCI_'):
                continue
            sampling_key = cache_key.split('_')[-1]
            tci_channels = tci_data.measurements['ne']['channels']
            ne_2d = tci_data.measurements['ne']['data']
            for i, ch in enumerate(tci_channels):
                ch_label = f'TCI{ch:02d}'
                channels.append({
                    'entry': f'{shot:06d} ({ch_label}) [{sampling_key}]',
                    'label': f'{ch_label} [{sampling_key}]',
                    'time': tci_data.time,
                    'param1': None,
                    'param2': ne_2d[i, :],
                    'R': None,
                    'shot': shot,
                    'source': 'TCI',
                })

        # Sort by R when available (TS, ECE), then by entry for the rest (TCI)
        channels.sort(key=lambda c: (c.get('R') is None, c.get('R') or 0, c['entry']))
        return channels

    def plot_data(self):
        """Plot time traces with markers only (no lines)"""
        self.ax1.clear()
        self.ax2.clear()

        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_xlabel('Time [s]')
        self.ax2.set_ylabel(self.param2['label'])

        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            return

        te_max, ne_max = 0, 0
        colors = self._get_plot_colors(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                shot_number, radius_or_ch, source, ch_label, sampling_key = self._parse_entry(entry)
                color = colors[i]

                if source == 'TS':
                    # Plot Thomson time trace
                    cache_key = f'{shot_number}_TS'

                    if cache_key not in self.data:
                        data = self.data_loader.load_data(shot_number)
                        self.data[cache_key] = data
                    else:
                        data = self.data[cache_key]

                    radius_idx = np.argmin(np.abs(data.radius - radius_or_ch))
                    actual_R = data.radius[radius_idx]

                    x_data = data.time

                    Te_data, Te_err_upper, Te_err_lower = data.get_parameter_asymmetric('Te')
                    ne_data, ne_err_upper, ne_err_lower = data.get_parameter_asymmetric('ne')

                    Te_trace = Te_data[radius_idx, :]
                    Te_err_upper_trace = Te_err_upper[radius_idx, :]
                    Te_err_lower_trace = Te_err_lower[radius_idx, :]
                    ne_trace = ne_data[radius_idx, :]
                    ne_err_upper_trace = ne_err_upper[radius_idx, :]
                    ne_err_lower_trace = ne_err_lower[radius_idx, :]

                    te_max = max(te_max, np.nanpercentile(Te_trace, 98))
                    ne_max = max(ne_max, np.nanpercentile(ne_trace, 98))

                    label = f'#{shot_number} {actual_R*1e3:.0f}mm ({ch_label})'

                    # TS: filled circle marker with asymmetric error bars
                    self.ax1.errorbar(x_data, Te_trace,
                                     yerr=[Te_err_lower_trace, Te_err_upper_trace],
                                     fmt='o', capsize=3, color=color,
                                     markersize=3, label=label, zorder=10)
                    self.ax2.errorbar(x_data, ne_trace,
                                     yerr=[ne_err_lower_trace, ne_err_upper_trace],
                                     fmt='o', capsize=3, color=color,
                                     markersize=3, label=label, zorder=10)

                elif source == 'ECE':
                    # Plot ECE time trace with specific sampling rate
                    cache_key = f'{shot_number}_ECE_{sampling_key}'

                    if cache_key not in self.ece_data_cache:
                        # Load with specified sampling rate
                        loader = self._get_ece_loader()
                        sr = 0.01 if sampling_key == '100Hz' else 0.001
                        ece_data = loader.load_data(shot_number, sampling_rate=sr)
                        self.ece_data_cache[cache_key] = ece_data
                    else:
                        ece_data = self.ece_data_cache[cache_key]

                    radius_idx = np.argmin(np.abs(ece_data.radius - radius_or_ch))
                    actual_R = ece_data.radius[radius_idx]

                    x_data = ece_data.time
                    Te_data = ece_data.measurements['Te']['data']
                    Te_trace = Te_data[radius_idx, :]

                    valid_mask = ece_data.measurements['Te']['valid_mask']
                    is_valid = valid_mask[radius_idx]

                    te_max = max(te_max, np.nanmax(Te_trace))

                    label = f'#{shot_number} {actual_R*1e3:.0f}mm ({ch_label}) [{sampling_key}]'
                    if not is_valid:
                        label += ' ovlp'

                    # ECE valid: empty square, ECE overlap: x marker
                    if is_valid:
                        self.ax1.plot(x_data, Te_trace, 's', color=color,
                                     markersize=3, linestyle='',
                                     markerfacecolor='none', markeredgewidth=1.5,
                                     label=label, zorder=5)
                    else:
                        self.ax1.plot(x_data, Te_trace, 'x', color=color,
                                     markersize=3, linestyle='',
                                     alpha=0.5, markeredgewidth=1.5,
                                     label=label, zorder=5)

                elif source == 'TS_NE_AVG':
                    # Plot TS_NE_AVG time trace (line-averaged ne)
                    if shot_number not in self.ts_ne_avg_cache:
                        # Re-load Thomson data to get ne_avg
                        ts_data = self.data_loader.load_data(shot_number)
                        cache_key = f'{shot_number}_TS'
                        self.data[cache_key] = ts_data
                        if 'ne_avg' in ts_data.measurements:
                            ne_avg = ts_data.measurements['ne_avg']
                            self.ts_ne_avg_cache[shot_number] = {
                                'time': ne_avg['time'],
                                'data': ne_avg['data']
                            }

                    cached = self.ts_ne_avg_cache.get(shot_number)
                    if cached is None:
                        continue
                    x_data = cached['time']
                    ne_trace = cached['data']

                    ne_max = max(ne_max, np.nanmax(ne_trace))

                    label = f'#{shot_number} (TS_NE_AVG)'

                    # TS_NE_AVG: solid line on ne axis
                    self.ax2.plot(x_data, ne_trace, '-', color=color,
                                 linewidth=1.5, label=label, zorder=5)

                elif source == 'TCI':
                    # Plot TCI time trace with specific sampling rate
                    cache_key = f'{shot_number}_TCI_{sampling_key}'

                    if cache_key not in self.tci_data_cache:
                        # Load with specified sampling rate
                        loader = self._get_tci_loader()
                        sr = 0.01 if sampling_key == '100Hz' else 0.001
                        tci_data = loader.load_data(shot_number, sampling_rate=sr)
                        self.tci_data_cache[cache_key] = tci_data
                    else:
                        tci_data = self.tci_data_cache[cache_key]

                    channel = int(radius_or_ch)  # For TCI, this is channel number
                    channels = tci_data.measurements['ne']['channels']

                    if channel in channels:
                        ch_idx = channels.index(channel)
                        x_data = tci_data.time
                        ne_data = tci_data.measurements['ne']['data']
                        ne_trace = ne_data[ch_idx, :]

                        ne_max = max(ne_max, np.nanmax(ne_trace))

                        # TCI has no R position
                        label = f'#{shot_number} ({ch_label}) [{sampling_key}]'

                        # TCI: solid line (no markers)
                        self.ax2.plot(x_data, ne_trace, '-', color=color,
                                     linewidth=1.5, label=label, zorder=5)

            except Exception as e:
                print(f"[ne/Te] Error plotting {entry}: {str(e)}")

        # Set limits
        if te_max > 0:
            te_margin = te_max * 0.1
            self.ax1.set_ylim(0, te_max + te_margin)

        if ne_max > 0:
            ne_margin = ne_max * 0.1
            self.ax2.set_ylim(0, ne_max + ne_margin)

        self._finalize_plot()

    def _write_data_to_file(self, file_path, selected_entries):
        """Write ne, Te time trace data to text file"""
        with open(file_path, 'w') as f:
            # Write header with aligned columns (10 char width each)
            f.write("# ne, Te Time Trace Data\n")
            f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "Te[keV]", "Te_err_u", "Te_err_l",
                "ne[e19]", "ne_err_u", "ne_err_l", "Source"
            ))

            for entry in selected_entries:
                shot_number, radius, source, ch_label, sampling_key = self._parse_entry(entry)

                if source == 'TS':
                    cache_key = f'{shot_number}_TS'
                    if cache_key not in self.data:
                        data = self.data_loader.load_data(shot_number)
                        self.data[cache_key] = data
                    else:
                        data = self.data[cache_key]

                    radius_idx = np.argmin(np.abs(data.radius - radius))
                    actual_R = data.radius[radius_idx]

                    Te_data, Te_err_upper, Te_err_lower = data.get_parameter_asymmetric('Te')
                    ne_data, ne_err_upper, ne_err_lower = data.get_parameter_asymmetric('ne')

                    Te_trace = Te_data[radius_idx, :]
                    Te_err_upper_trace = Te_err_upper[radius_idx, :]
                    Te_err_lower_trace = Te_err_lower[radius_idx, :]
                    ne_trace = ne_data[radius_idx, :]
                    ne_err_upper_trace = ne_err_upper[radius_idx, :]
                    ne_err_lower_trace = ne_err_lower[radius_idx, :]

                    # Write Thomson data rows (10 char width each)
                    for i in range(len(data.time)):
                        f.write(" %10d,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10s\n" % (
                            shot_number, data.time[i], actual_R,
                            Te_trace[i], Te_err_upper_trace[i], Te_err_lower_trace[i],
                            ne_trace[i], ne_err_upper_trace[i], ne_err_lower_trace[i], 'TS'
                        ))

                elif source == 'ECE':
                    cache_key = f'{shot_number}_ECE_{sampling_key}'
                    if cache_key not in self.ece_data_cache:
                        loader = self._get_ece_loader()
                        sr = 0.01 if sampling_key == '100Hz' else 0.001
                        ece_data = loader.load_data(shot_number, sampling_rate=sr)
                        self.ece_data_cache[cache_key] = ece_data
                    else:
                        ece_data = self.ece_data_cache[cache_key]

                    radius_idx = np.argmin(np.abs(ece_data.radius - radius))
                    actual_R = ece_data.radius[radius_idx]

                    Te_data = ece_data.measurements['Te']['data']
                    Te_trace = Te_data[radius_idx, :]

                    # Write ECE data rows (no ne for ECE, 10 char width each)
                    for i in range(len(ece_data.time)):
                        f.write(" %10d,%10.3f,%10.3f,%10.3f,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                            shot_number, ece_data.time[i], actual_R,
                            Te_trace[i], 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'ECE'
                        ))

                elif source == 'TS_NE_AVG':
                    if shot_number in self.ts_ne_avg_cache:
                        cached = self.ts_ne_avg_cache[shot_number]
                        for i in range(len(cached['time'])):
                            f.write(" %10d,%10.3f,%10s,%10s,%10s,%10s,%10.3f,%10s,%10s,%10s\n" % (
                                shot_number, cached['time'][i], 'NaN',
                                'NaN', 'NaN', 'NaN', cached['data'][i], 'NaN', 'NaN', 'TS_AVG'
                            ))

                elif source == 'TCI':
                    cache_key = f'{shot_number}_TCI_{sampling_key}'
                    if cache_key not in self.tci_data_cache:
                        loader = self._get_tci_loader()
                        sr = 0.01 if sampling_key == '100Hz' else 0.001
                        tci_data = loader.load_data(shot_number, sampling_rate=sr)
                        self.tci_data_cache[cache_key] = tci_data
                    else:
                        tci_data = self.tci_data_cache[cache_key]

                    channel = int(radius)  # For TCI, radius is channel number
                    channels = tci_data.measurements['ne']['channels']

                    if channel in channels:
                        ch_idx = channels.index(channel)
                        ne_data = tci_data.measurements['ne']['data']
                        ne_trace = ne_data[ch_idx, :]

                        # Write TCI data rows (no Te, no R for TCI)
                        for i in range(len(tci_data.time)):
                            f.write(" %10d,%10.3f,%10s,%10s,%10s,%10s,%10.3f,%10s,%10s,%10s\n" % (
                                shot_number, tci_data.time[i], 'NaN',
                                'NaN', 'NaN', 'NaN', ne_trace[i], 'NaN', 'NaN', 'TCI'
                            ))
