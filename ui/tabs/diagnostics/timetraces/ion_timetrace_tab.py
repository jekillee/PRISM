"""
Ion (Ti, vT) Time Trace tab with CES and XICS integration
"""

import re
import numpy as np

from PySide6.QtWidgets import (
    QMessageBox, QLineEdit, QComboBox, QPushButton, QWidget,
    QHBoxLayout, QVBoxLayout, QGroupBox, QGridLayout, QLabel,
    QFileDialog, QStyle
)

from ui.tabs.timetrace_base_tab import TimeTraceBaseTab
from ui.ui_constants import apply_shot_arrow_icons
from ui.theme import ThemeManager


class IonTimeTraceTab(TimeTraceBaseTab):
    """Ion (Ti, vT) Time Trace tab with CES and XICS support"""

    TAB_NAME = "Ion"
    Y2_VT_LABEL = 'v$_T$ [km/s]'
    Y2_OMEGA_LABEL = r'$\omega_T$ [krad/s]'

    def __init__(self, parent, app_config, diagnostic_name, tab_type,
                 data_loader, efit_loader, plot_manager, file_parser):
        super().__init__(parent, app_config, diagnostic_name, tab_type,
                        data_loader, efit_loader, plot_manager)
        self.file_parser = file_parser
        self.xics_loader = None
        self.xics_data_cache = {}

    # ---------- vT / ω display helpers ----------

    def _is_omega_mode(self):
        return getattr(self, 'param2_combo', None) is not None \
            and self.param2_combo.currentText().startswith('ωT')

    def _y2_label(self):
        return self.Y2_OMEGA_LABEL if self._is_omega_mode() else self.Y2_VT_LABEL

    def _to_y2_trace(self, vT_trace, vT_err_trace, R_scalar):
        """Time trace: divide by scalar R per channel when in omega mode."""
        if not self._is_omega_mode():
            return vT_trace, vT_err_trace
        if R_scalar is None or R_scalar <= 0:
            return np.full_like(vT_trace, np.nan), np.full_like(vT_err_trace, np.nan)
        return vT_trace / R_scalar, vT_err_trace / R_scalar

    def _create_plot_controls(self, parent):
        """Plot controls: Plot/Option + vT/ω toggle (mirrors MSE q/j pattern)."""
        from PySide6.QtWidgets import QFrame as _QFrame, QGroupBox as _QGB
        group = _QGB("3. Plot")
        group_layout = QVBoxLayout(group)

        row1 = QHBoxLayout()
        self.param2_combo = QComboBox()
        self.param2_combo.addItems(['vT [km/s]', 'ωT [krad/s]'])
        self.param2_combo.setCurrentText('vT [km/s]')
        self.param2_combo.setFixedWidth(110)
        self.param2_combo.currentTextChanged.connect(lambda _t: self.plot_data())
        row1.addWidget(self.param2_combo)
        plot_button = QPushButton("Plot")
        plot_button.clicked.connect(self.plot_data)
        row1.addWidget(plot_button, 2)
        style_btn = QPushButton("Option")
        style_btn.clicked.connect(self._show_style_dialog)
        row1.addWidget(style_btn, 1)
        group_layout.addLayout(row1)

        self.color_mode_combo = QComboBox()
        self.color_mode_combo.addItems([
            "Gradient(viridis)", "Gradient(hot)", "Gradient(jet)", "Gradient(coolwarm)",
            "Fixed(tab10)", "Fixed(tab20)", "Fixed(Set1)", "Fixed(Set2)", "Fixed(Set3)",
        ])
        self.color_mode_combo.setCurrentText("Gradient(viridis)")
        self.color_mode_combo.hide()
        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10

        parent.layout().addWidget(group)

    def _restore_shot_from_settings(self):
        super()._restore_shot_from_settings()
        from config.user_settings import get_tab_settings
        s = get_tab_settings(self._settings_key)
        if hasattr(self, 'param2_combo'):
            saved = s.get('y2_mode', 'vT [km/s]')
            idx = self.param2_combo.findText(saved)
            if idx >= 0:
                self.param2_combo.setCurrentIndex(idx)

    def save_settings(self):
        super().save_settings()
        from config.user_settings import get_tab_settings, set_tab_settings
        s = get_tab_settings(self._settings_key)
        if hasattr(self, 'param2_combo'):
            s['y2_mode'] = self.param2_combo.currentText()
        set_tab_settings(self._settings_key, s)

    def _create_shot_input(self, parent):
        """Create data loading section with analysis type selection"""
        group = QGroupBox("1. Load Ti, vT Data")
        grid = QGridLayout(group)

        grid.setColumnStretch(1, 1)

        # Row 0: Shot label, entry, up/down buttons, dropdown, Fetch, File
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

        analysis_types = list(self.diag_config['analysis_types'].keys())
        self.analysis_combo = QComboBox()
        self.analysis_combo.addItems(analysis_types)
        self.analysis_combo.setCurrentText('mod')
        self.analysis_combo.setFixedWidth(110)
        grid.addWidget(self.analysis_combo, 0, 3)

        self.fetch_button = QPushButton('Fetch')
        self.fetch_button.setFixedWidth(70)
        self.fetch_button.clicked.connect(self.load_shot_data)
        grid.addWidget(self.fetch_button, 0, 4)

        from ui.ui_constants import make_or_divider
        grid.addWidget(make_or_divider(), 1, 0, 1, 5)
        file_button = QPushButton("Open CES Result File...")
        file_button.clicked.connect(self.load_file_data)
        grid.addWidget(file_button, 2, 0, 1, 5)

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

    def _get_xics_loader(self):
        """Lazy initialization of XICS loader"""
        if self.xics_loader is None:
            from data_loaders.xics_loader import XICSLoader
            from config.diagnostic_config import DIAGNOSTICS
            self.xics_loader = XICSLoader(self.app_config, DIAGNOSTICS.get('XICS', {}))
        return self.xics_loader

    def load_shot_data(self):
        """Load CES and XICS shot data from MDS+, sorted by R"""
        try:
            shot_number = int(self.shot_entry.text())
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
            return

        analysis_type = self.analysis_combo.currentText()
        self._set_status(f"Loading #{shot_number} ({analysis_type})...", 'blue')
        ces_loaded = False
        xics_loaded = False
        data = None
        xics_data = None

        # Try to load CES data
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
            self._set_status(
                f"No CES ({analysis_type}) or XICS data for #{shot_number}", 'red')
            QMessageBox.critical(self.frame, "Error", f"No CES ({analysis_type}) or XICS data available for shot #{shot_number}")
            return

        sources = []
        if ces_loaded:
            sources.append(f"CES ({analysis_type})")
        if xics_loaded:
            sources.append("XICS")
        self._set_status(f"#{shot_number} loaded: {', '.join(sources)}", 'green')

        # Build list of all channels with R positions
        all_channels = []

        # Add CES channels: mod01-32 or nn01-32
        if ces_loaded:
            for i, radius in enumerate(data.radius):
                ch_label = f'{analysis_type}{i+1:02d}'
                all_channels.append({
                    'R': radius,
                    'label': f'{shot_number:06d}_{radius*1e3:.0f} ({ch_label})'
                })

        # Add XICS channel
        if xics_loaded:
            all_channels.append({
                'R': 1.8,
                'label': f'{shot_number:06d}_1800 (XICS)'
            })

        # Sort by R position
        all_channels.sort(key=lambda x: x['R'])

        # Update listbox
        self.available_listbox.clear()
        for ch in all_channels:
            self.available_listbox.addItem(ch['label'])

        if hasattr(self, 'browse_button'):
            self.browse_button.setEnabled(True)
            self.browse_button.setText(f"Browse #{shot_number}")
            self._last_shot = shot_number
            self._last_source = analysis_type

    def load_file_data(self):
        """Load CES data from result file, sorted by R"""
        file_path = QFileDialog.getOpenFileName(
            self.frame,
            "Select CES Result File",
            self.app_config.CES_RESULT_PATH,
            "Text files (*.txt);;All Files (*)"
        )[0]

        if not file_path:
            return

        try:
            data, shot_number = self.file_parser.parse_file(file_path)

            file_key = f'file_{shot_number}_{data.source}'
            self.data[file_key] = data

            print(f"[CES] Result file loaded: {data.source}")
            self._set_status(f"File loaded: {data.source} for #{shot_number}", 'green')

            # Try to load XICS data for this shot
            loader = self._get_xics_loader()
            xics_data = loader.load_data(shot_number)

            if xics_data is not None:
                xics_cache_key = f'{shot_number}_XICS'
                self.xics_data_cache[xics_cache_key] = xics_data
                print(f"[XICS] Data loaded: {len(xics_data.time)} timepoints")

            # Build list of all channels with R positions
            all_channels = []

            # Add CES file channels
            for i, radius in enumerate(data.radius):
                ch_label = f'{data.source}{i+1:02d}'
                all_channels.append({
                    'R': radius,
                    'label': f'{shot_number:06d}_{radius*1e3:.0f} ({ch_label})'
                })

            # Add XICS channel
            if xics_data is not None:
                all_channels.append({
                    'R': 1.8,
                    'label': f'{shot_number:06d}_1800 (XICS)'
                })

            # Sort by R position
            all_channels.sort(key=lambda x: x['R'])

            # Update listbox
            self.available_listbox.clear()
            for ch in all_channels:
                self.available_listbox.addItem(ch['label'])

        except ValueError as e:
            QMessageBox.critical(self.frame, "Invalid File Format", str(e))
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load file: {str(e)}")

    def _parse_entry(self, entry):
        """Parse time trace entry: XXXXXX_YYYY (DIAG##)

        Returns: shot_number, radius, source, channel_label
        """
        # Format: 012345_1800 (mod01) or 012345_1800 (XICS)
        main_part = entry.split('(')[0].strip()
        diag_label = entry.split('(')[1].split(')')[0].strip()

        parts = main_part.split('_')
        shot_number = int(parts[0])
        radius = float(parts[1]) / 1e3  # mm to m

        # Determine source from diag_label (strip trailing digits)
        if diag_label.startswith('mod'):
            source = 'mod'
        elif diag_label.startswith('nn'):
            source = 'nn'
        elif diag_label == 'XICS':
            source = 'XICS'
        else:
            # File data: strip channel number (e.g., 'tgf01' → 'tgf')
            source = re.sub(r'\d+$', '', diag_label)

        return shot_number, radius, source, diag_label

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

    def _get_timetrace_preview_channels(self):
        shot = getattr(self, '_last_shot', None)
        source = getattr(self, '_last_source', None)
        if shot is None:
            return []

        channels = []
        cache_key = f'{shot}_{source}'
        data = self.data.get(cache_key)
        if data is not None:
            Ti_data, _ = data.get_parameter('Ti')
            vT_data, _ = data.get_parameter('vT')
            for i, R in enumerate(data.radius):
                ch_label = f'{source}{i+1:02d}'
                channels.append({
                    'entry': f'{shot:06d}_{R*1e3:.0f} ({ch_label})',
                    'label': f'{ch_label} R={R*1e3:.0f}mm',
                    'time': data.time,
                    'param1': Ti_data[i, :],
                    'param2': vT_data[i, :],
                    'shot': shot,
                    'source': source,
                })

        xics_key = f'{shot}_XICS'
        xics = self.xics_data_cache.get(xics_key)
        if xics is not None:
            Ti_x, _ = xics.get_parameter('Ti')
            vT_x, _ = xics.get_parameter('vT')
            channels.append({
                'entry': f'{shot:06d}_1800 (XICS)',
                'label': 'XICS R=1800mm',
                'time': xics.time,
                'param1': Ti_x[0, :],
                'param2': vT_x[0, :],
                'shot': shot,
                'source': 'XICS',
            })

        channels.sort(key=lambda c: float(c['entry'].split('_')[1].split()[0]))
        return channels

    def plot_data(self):
        """Plot time traces with markers"""
        self.ax1.clear()
        self.ax2.clear()

        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_xlabel('Time [s]')
        self.ax2.set_ylabel(self._y2_label())

        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            return

        ti_max, vt_max, vt_min = 0, 0, 0
        colors = self._get_plot_colors(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                shot_number, radius, source, ch_label = self._parse_entry(entry)
                color = colors[i]

                if source == 'XICS':
                    # Plot XICS time trace
                    xics_cache_key = f'{shot_number}_XICS'

                    if xics_cache_key not in self.xics_data_cache:
                        continue

                    xics_data = self.xics_data_cache[xics_cache_key]

                    x_data = xics_data.time
                    Ti_data, Ti_err = xics_data.get_parameter('Ti')
                    vT_data, vT_err = xics_data.get_parameter('vT')

                    Ti_trace = Ti_data[0, :]
                    vT_trace = vT_data[0, :]
                    vT_err_x = vT_err[0, :] if vT_err is not None else 0 * vT_trace
                    # XICS R fixed at 1.8 m (TXCS calibration); apply ω conversion if needed
                    vT_trace, _ = self._to_y2_trace(vT_trace, vT_err_x, 1.8)

                    ti_max = max(ti_max, np.nanpercentile(Ti_trace, 98))
                    vt_min = min(vt_min, np.nanpercentile(vT_trace, 2))
                    vt_max = max(vt_max, np.nanpercentile(vT_trace, 98))

                    label = f'#{shot_number} 1800mm ({ch_label})'

                    style = self._get_marker_style('XICS')

                    # XICS: closed square without error bars
                    self.ax1.plot(x_data, Ti_trace,
                                 marker=style['marker'], color=color,
                                 markersize=3, label=label,
                                 linestyle='none')
                    self.ax2.plot(x_data, vT_trace,
                                 marker=style['marker'], color=color,
                                 markersize=3, label=label,
                                 linestyle='none')

                else:
                    # Plot CES time trace
                    if source in ['mod', 'nn']:
                        cache_key = f'{shot_number}_{source}'
                    else:
                        cache_key = f'file_{shot_number}_{source}'

                    # Load data if not cached
                    if cache_key not in self.data:
                        if source in ['mod', 'nn']:
                            data = self.data_loader.load_data(shot_number, source)
                        else:
                            continue
                        self.data[cache_key] = data
                    else:
                        data = self.data[cache_key]

                    radius_idx = np.argmin(np.abs(data.radius - radius))
                    actual_R = data.radius[radius_idx]

                    x_data = data.time

                    Ti_data, Ti_err = data.get_parameter('Ti')
                    vT_data, vT_err = data.get_parameter('vT')

                    Ti_trace = Ti_data[radius_idx, :]
                    Ti_err_trace = Ti_err[radius_idx, :]
                    vT_trace = vT_data[radius_idx, :]
                    vT_err_trace = vT_err[radius_idx, :]

                    # ω conversion: divide by channel's R when omega mode is on
                    vT_trace, vT_err_trace = self._to_y2_trace(vT_trace, vT_err_trace, actual_R)

                    ti_max = max(ti_max, np.nanpercentile(Ti_trace, 98))
                    vt_min = min(vt_min, np.nanpercentile(vT_trace, 2))
                    vt_max = max(vt_max, np.nanpercentile(vT_trace, 98))

                    label = f'#{shot_number} {actual_R*1e3:.0f}mm ({ch_label})'

                    style = self._get_marker_style(source)

                    plot_kwargs = {
                        'fmt': style['marker'],
                        'capsize': 3,
                        'label': label,
                        'color': color,
                        'markersize': 3,
                        'fillstyle': style['fillstyle']
                    }
                    if style['markerfacecolor'] == 'none':
                        plot_kwargs['markerfacecolor'] = 'none'
                        plot_kwargs['markeredgewidth'] = style['markeredgewidth']

                    self.ax1.errorbar(x_data, Ti_trace, Ti_err_trace, **plot_kwargs)
                    self.ax2.errorbar(x_data, vT_trace, vT_err_trace, **plot_kwargs)

            except Exception as e:
                print(f"[Ti/vT] Error plotting {entry}: {str(e)}")

        # Set limits
        ti_margin = ti_max * 0.1
        self.ax1.set_ylim(0, ti_max + ti_margin)

        vt_margin = (vt_max - vt_min) * 0.1
        self.ax2.set_ylim(vt_min - vt_margin, vt_max + vt_margin)
        zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'
        self.ax2.axhline(y=0, c=zc, ls='--', gid='zero_ref')

        self._finalize_plot()

    def _write_data_to_file(self, file_path, selected_entries):
        """Write Ti/vT time trace data to text file"""
        with open(file_path, 'w') as f:
            # Write header
            f.write("# Ti/vT Time Trace Data\n")
            f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "Ti[keV]", "Ti_err", "vT[km/s]", "vT_err", "Source"
            ))

            for entry in selected_entries:
                shot_number, radius, source, ch_label = self._parse_entry(entry)

                if source == 'XICS':
                    # Write XICS data
                    xics_cache_key = f'{shot_number}_XICS'
                    if xics_cache_key not in self.xics_data_cache:
                        continue

                    xics_data = self.xics_data_cache[xics_cache_key]

                    Ti_data, Ti_err = xics_data.get_parameter('Ti')
                    vT_data, vT_err = xics_data.get_parameter('vT')

                    Ti_trace = Ti_data[0, :]
                    Ti_err_trace = Ti_err[0, :]
                    vT_trace = vT_data[0, :]
                    vT_err_trace = vT_err[0, :]

                    for i in range(len(xics_data.time)):
                        f.write(" %10d,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10s\n" % (
                            shot_number, xics_data.time[i], 1.8,
                            Ti_trace[i], Ti_err_trace[i],
                            vT_trace[i], vT_err_trace[i], 'XICS'
                        ))

                else:
                    # Write CES data
                    if source in ['mod', 'nn']:
                        cache_key = f'{shot_number}_{source}'
                    else:
                        cache_key = f'file_{shot_number}_{source}'

                    if cache_key not in self.data:
                        if source in ['mod', 'nn']:
                            data = self.data_loader.load_data(shot_number, source)
                            self.data[cache_key] = data
                        else:
                            continue
                    else:
                        data = self.data[cache_key]

                    radius_idx = np.argmin(np.abs(data.radius - radius))
                    actual_R = data.radius[radius_idx]

                    Ti_data, Ti_err = data.get_parameter('Ti')
                    vT_data, vT_err = data.get_parameter('vT')

                    Ti_trace = Ti_data[radius_idx, :]
                    Ti_err_trace = Ti_err[radius_idx, :]
                    vT_trace = vT_data[radius_idx, :]
                    vT_err_trace = vT_err[radius_idx, :]

                    for i in range(len(data.time)):
                        f.write(" %10d,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10s\n" % (
                            shot_number, data.time[i], actual_R,
                            Ti_trace[i], Ti_err_trace[i],
                            vT_trace[i], vT_err_trace[i], source
                        ))
