"""
Spectrogram tab for high-frequency signal analysis
Supports ECE, Mirnov, and BES diagnostics
"""

import os
import time as tclock
import numpy as np
from scipy.signal import spectrogram
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from MDSplus import Connection

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QGroupBox, QLabel, QLineEdit, QPushButton, QComboBox,
    QSlider, QMessageBox, QFileDialog, QApplication,
    QSplitter, QDialog, QTextEdit, QStyle, QScrollArea,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QColor, QTextCharFormat, QSyntaxHighlighter, QGuiApplication

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style, get_icon
from config.user_settings import get_tab_settings, set_tab_settings
from data_loaders.ecei_loader import ECEILoader
from data_loaders.ece_loader import ECELoader
from data_loaders.bes_loader import BESLoader


# Example script for loading spectrogram NPZ file
SPECTROGRAM_EXAMPLE_SCRIPT = '''"""
Example script for loading and plotting PRISM spectrogram NPZ file
"""

import numpy as np
import matplotlib.pyplot as plt

# Load NPZ file
filepath = 'spectrogram_XXXXXX_ECEYY_0.0-10.0s.npz'
data = np.load(filepath, allow_pickle=True)

# Extract data
time = data['time']           # Time array [s]
frequency = data['frequency'] # Frequency array [kHz]
power = data['power']         # Power array (time, freq) [log10 scale]

# Get metadata
metadata = data['metadata'].item()
print(f"Shot: {metadata['shot']}")
print(f"Channel: {metadata['channel']}")
print(f"Diagnostic: {metadata['diagnostic']}")
print(f"Time range: {metadata['time_range']} s")
print(f"Freq range: {metadata['freq_range']} kHz")

# Plot spectrogram
fig, ax = plt.subplots(figsize=(10, 6))
im = ax.pcolormesh(time, frequency, power.T, shading='auto')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Frequency [kHz]')
ax.set_title(f"Shot #{metadata['shot']} - {metadata['channel']}")
plt.colorbar(im, ax=ax, label='log$_{10}$(Power) [a.u.]')
plt.tight_layout()
plt.show()
'''


class SpectrogramTab:
    """Spectrogram visualization tab"""

    # NFFT options
    NFFT_OPTIONS = ['256', '512', '1024', '2048', '4096']
    DEFAULT_NFFT = '1024'

    # Label column width for consistent alignment
    LABEL_COLUMN_WIDTH = 90

    def __init__(self, parent, app_config, diagnostic_config):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config

        self.frame = QWidget()
        self.toolbar = None

        # ECEI loader
        self.ecei_loader = ECEILoader(app_config.MDS_IP)

        # ECE loader
        self.ece_loader = ECELoader(app_config, None)

        # BES loader
        self.bes_loader = BESLoader(app_config, None)

        # Data cache
        self.ece_info = {}      # {shot: {'channels': [], 'R': [], 'Z': [], 'I_TF': float}}
        self.bes_info = {}      # {shot: {'channels': [], 'R': [], 'Z': []}}
        self.ecei_info = {}     # {(shot, device): {'channels': [], 'R': [], 'Z': []}}
        self.signal_data = None
        self.spectrogram_data = None

        # Current loaded shot
        self.current_shot = None

        # Plot references
        self.im = None
        self.colorbar = None

        # Widget references for enable/disable
        self.signal_widgets = []

    def create_widgets(self):
        """Create spectrogram tab widgets"""
        # Single canvas for spectrogram
        self.figure = Figure(self.app_config.FIGURE_SIZE, tight_layout=True)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xlabel('Time [s]')
        self.ax.set_ylabel('Frequency [kHz]')
        apply_dark_figure_style(self.figure)

        # Create canvas
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.draw()

        # Left side: canvas
        canvas_widget = QWidget()
        canvas_layout = QVBoxLayout(canvas_widget)
        canvas_layout.setContentsMargins(0, 0, 0, 0)
        canvas_layout.addWidget(self.canvas)

        # Right side: scrollable control panel with fixed width
        scroll_area = QScrollArea()
        scroll_area.setFixedWidth(CONTROL_PANEL_WIDTH)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)

        control_frame = QWidget()
        control_layout = QVBoxLayout(control_frame)
        control_layout.setContentsMargins(9, 9, 9, 9)
        control_layout.setSizeConstraint(QVBoxLayout.SetNoConstraint)

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
        self._create_signal_selection(control_frame)
        self._create_spectrogram_params(control_frame)
        self._create_color_controls(control_frame)
        self._create_save_controls(control_frame)

        # Add stretch at the bottom of control panel
        control_layout.addStretch()

        # Initially disable signal selection
        self._set_signal_selection_state(False)

        # Load saved settings
        self.load_settings()

    def _create_shot_input(self, parent):
        """Create shot input section with up/down buttons"""
        frame = QGroupBox("1. Load Shot", parent)
        grid = QGridLayout(frame)

        grid.setColumnMinimumWidth(0, self.LABEL_COLUMN_WIDTH)
        grid.setColumnStretch(1, 1)

        grid.addWidget(QLabel('Shot'), 0, 0)

        # Shot entry with up/down buttons in same row
        shot_widget = QWidget()
        shot_layout = QHBoxLayout(shot_widget)
        shot_layout.setContentsMargins(0, 0, 0, 0)

        self.shot_entry = QLineEdit()
        shot_layout.addWidget(self.shot_entry, 1)
        self.shot_entry.returnPressed.connect(self._load_shot_info)

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

        grid.addWidget(shot_widget, 0, 1)

        self.fetch_button = QPushButton('Fetch')
        self.fetch_button.clicked.connect(self._load_shot_info)
        grid.addWidget(self.fetch_button, 0, 2)

        parent.layout().addWidget(frame)

    def _adjust_shot(self, delta):
        """Adjust shot number by delta"""
        try:
            current = int(self.shot_entry.text())
            new_shot = max(1, current + delta)
            self.shot_entry.setText(str(new_shot))
        except ValueError:
            pass

    def _create_signal_selection(self, parent):
        """Create signal selection section"""
        frame = QGroupBox("2. Select Signal", parent)
        grid = QGridLayout(frame)
        self.signal_frame = frame

        grid.setColumnMinimumWidth(0, self.LABEL_COLUMN_WIDTH)
        grid.setColumnStretch(1, 1)

        # Shot number display (row 0)
        grid.addWidget(QLabel('Loaded Shot'), 0, 0)

        self.shot_display_label = QLabel('Not loaded')
        self.shot_display_label.setStyleSheet('color: gray;')
        grid.addWidget(self.shot_display_label, 0, 1)

        # Diagnostic type dropdown (row 1)
        grid.addWidget(QLabel('Diagnostic'), 1, 0)

        self.diag_options = ['Mirnov (Toroidal)', 'Mirnov (Poloidal)', 'ECE', 'BES', 'TCI',
                             'ECEI-GT', 'ECEI-GR', 'ECEI-HT']
        self.selected_diag_value = 'Mirnov (Toroidal)'

        self.diag_dropdown = QComboBox()
        self.diag_dropdown.addItems(self.diag_options)
        self.diag_dropdown.setCurrentText(self.selected_diag_value)
        self.diag_dropdown.setEnabled(False)
        self.diag_dropdown.currentTextChanged.connect(self._on_diag_changed)
        grid.addWidget(self.diag_dropdown, 1, 1)
        self.signal_widgets.append(self.diag_dropdown)

        # Channel dropdown (row 2)
        grid.addWidget(QLabel('Channel'), 2, 0)

        self.channel_dropdown = QComboBox()
        self.channel_dropdown.setEnabled(False)
        grid.addWidget(self.channel_dropdown, 2, 1)
        self.signal_widgets.append(self.channel_dropdown)

        parent.layout().addWidget(frame)

    def _create_spectrogram_params(self, parent):
        """Create spectrogram parameter section"""
        frame = QGroupBox("3. Spectrogram Parameters", parent)
        grid = QGridLayout(frame)

        grid.setColumnMinimumWidth(0, self.LABEL_COLUMN_WIDTH)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 1)

        # Time range
        grid.addWidget(QLabel('Time [s]'), 0, 0)

        time_widget = QWidget()
        time_layout = QHBoxLayout(time_widget)
        time_layout.setContentsMargins(0, 0, 0, 0)

        self.time_min_entry = QLineEdit()
        self.time_min_entry.setFixedWidth(60)
        self.time_min_entry.setText('0')
        time_layout.addWidget(self.time_min_entry)

        time_layout.addWidget(QLabel('-'))

        self.time_max_entry = QLineEdit()
        self.time_max_entry.setFixedWidth(60)
        self.time_max_entry.setText('10')
        time_layout.addWidget(self.time_max_entry)
        time_layout.addStretch()

        grid.addWidget(time_widget, 0, 1, 1, 3)

        # Frequency range
        grid.addWidget(QLabel('Freq [kHz]'), 1, 0)

        freq_widget = QWidget()
        freq_layout = QHBoxLayout(freq_widget)
        freq_layout.setContentsMargins(0, 0, 0, 0)

        self.freq_min_entry = QLineEdit()
        self.freq_min_entry.setFixedWidth(60)
        self.freq_min_entry.setText('0')
        freq_layout.addWidget(self.freq_min_entry)

        freq_layout.addWidget(QLabel('-'))

        self.freq_max_entry = QLineEdit()
        self.freq_max_entry.setFixedWidth(60)
        self.freq_max_entry.setText('100')
        freq_layout.addWidget(self.freq_max_entry)
        freq_layout.addStretch()

        grid.addWidget(freq_widget, 1, 1, 1, 3)

        # NFFT
        grid.addWidget(QLabel('NFFT'), 2, 0)

        self.nfft_dropdown = QComboBox()
        self.nfft_dropdown.addItems(self.NFFT_OPTIONS)
        self.nfft_dropdown.setCurrentText(self.DEFAULT_NFFT)
        grid.addWidget(self.nfft_dropdown, 2, 1)

        # Plot button
        plot_btn = QPushButton('Plot Spectrogram')
        plot_btn.clicked.connect(self._plot_spectrogram)
        grid.addWidget(plot_btn, 3, 0, 1, 4)

        # Status label
        self.status_label = QLabel('Ready')
        self.status_label.setStyleSheet('color: gray; font-weight: bold;')
        grid.addWidget(self.status_label, 4, 0, 1, 4)

        parent.layout().addWidget(frame)

    def _create_color_controls(self, parent):
        """Create color range control section"""
        frame = QGroupBox("4. Color Range", parent)
        layout = QVBoxLayout(frame)

        # Dynamic range slider (1 to 11 decades)
        layout.addWidget(QLabel('Dynamic Range (decades):'))

        slider_widget = QWidget()
        slider_layout = QHBoxLayout(slider_widget)
        slider_layout.setContentsMargins(0, 0, 0, 0)

        self._dyn_range_value = 6.0
        self.dyn_range_slider = QSlider(Qt.Horizontal)
        self.dyn_range_slider.setMinimum(10)   # represents 1.0
        self.dyn_range_slider.setMaximum(110)  # represents 11.0
        self.dyn_range_slider.setValue(60)      # represents 6.0
        self.dyn_range_slider.valueChanged.connect(self._on_dyn_range_changed)
        slider_layout.addWidget(self.dyn_range_slider, 1)

        self.dyn_range_label = QLabel('6.0')
        self.dyn_range_label.setFixedWidth(35)
        slider_layout.addWidget(self.dyn_range_label)

        layout.addWidget(slider_widget)

        # Labels for scale reference
        ref_widget = QWidget()
        ref_layout = QHBoxLayout(ref_widget)
        ref_layout.setContentsMargins(0, 0, 0, 0)
        lbl_1 = QLabel('1')
        lbl_1.setStyleSheet('font-size: 9px;')
        ref_layout.addWidget(lbl_1)
        ref_layout.addStretch()
        lbl_11 = QLabel('11')
        lbl_11.setStyleSheet('font-size: 9px;')
        ref_layout.addWidget(lbl_11)
        layout.addWidget(ref_widget)

        parent.layout().addWidget(frame)

    def _create_save_controls(self, parent):
        """Create save data section"""
        frame = QGroupBox("5. Save Data", parent)
        layout = QVBoxLayout(frame)

        btn_widget = QWidget()
        btn_layout = QHBoxLayout(btn_widget)
        btn_layout.setContentsMargins(0, 0, 0, 0)

        self.save_button = QPushButton('Save as .npz')
        self.save_button.clicked.connect(self._save_data)
        self.save_button.setEnabled(False)
        btn_layout.addWidget(self.save_button, 1)

        example_btn = QPushButton('Example Script')
        example_btn.clicked.connect(self._show_example_script)
        btn_layout.addWidget(example_btn, 1)

        layout.addWidget(btn_widget)

        parent.layout().addWidget(frame)

    def _show_example_script(self):
        """Show example script for loading NPZ file with syntax highlighting"""
        script = SPECTROGRAM_EXAMPLE_SCRIPT

        popup = QDialog(self.frame)
        popup.setWindowTitle("Example Script - Spectrogram NPZ")
        popup.resize(650, 500)

        dialog_layout = QVBoxLayout(popup)

        text_widget = QTextEdit()
        text_widget.setReadOnly(True)
        text_widget.setFont(QFont('Courier', 10))
        text_widget.setStyleSheet(
            'background-color: #19232d; color: #ffffff;'
        )
        text_widget.setPlainText(script)

        # Apply syntax highlighting
        self._highlighter = PythonHighlighter(text_widget.document())

        dialog_layout.addWidget(text_widget)

        btn_widget = QWidget()
        btn_layout = QHBoxLayout(btn_widget)
        btn_layout.setContentsMargins(0, 0, 0, 0)

        def copy_to_clipboard():
            clipboard = QGuiApplication.clipboard()
            clipboard.setText(script)
            QMessageBox.information(popup, "Copied", "Script copied to clipboard")

        copy_btn = QPushButton("Copy to Clipboard")
        copy_btn.clicked.connect(copy_to_clipboard)
        btn_layout.addWidget(copy_btn)

        btn_layout.addStretch()

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(popup.accept)
        btn_layout.addWidget(close_btn)

        dialog_layout.addWidget(btn_widget)
        popup.exec()

    def _save_data(self):
        """Save spectrogram data to NPZ file (only selected time/freq range)"""
        if self.spectrogram_data is None or self.signal_data is None:
            QMessageBox.warning(self.frame, "Warning", "No spectrogram data to save")
            return

        # Get parameters for filename
        shot = self.signal_data['shot']
        channel = self.signal_data['channel'].split('(')[0].strip()  # Remove position info
        t_min = float(self.time_min_entry.text())
        t_max = float(self.time_max_entry.text())
        f_min = float(self.freq_min_entry.text())  # kHz
        f_max = float(self.freq_max_entry.text())  # kHz

        # Default filename
        default_name = f"spectrogram_{shot}_{channel}_{t_min:.1f}-{t_max:.1f}s.npz"

        # File dialog
        filepath, _ = QFileDialog.getSaveFileName(
            self.frame,
            'Save Spectrogram Data',
            os.path.join(os.path.expanduser("~"), default_name),
            'NumPy NPZ (*.npz)'
        )

        if not filepath:
            return

        try:
            # Get full data
            f_full = self.spectrogram_data['f']  # Hz
            t_full = self.spectrogram_data['t']  # s
            Sxx_log_full = self.spectrogram_data['Sxx_log']  # (freq, time)

            # Apply frequency mask (convert kHz to Hz for comparison)
            f_mask = (f_full >= f_min * 1e3) & (f_full <= f_max * 1e3)

            # Slice data to selected frequency range
            f_sliced = f_full[f_mask] / 1e3  # Convert to kHz
            Sxx_log_sliced = Sxx_log_full[f_mask, :]  # (freq_sliced, time)

            # Prepare metadata
            metadata = {
                'shot': shot,
                'channel': self.signal_data['channel'],
                'diagnostic': self.diag_dropdown.currentText(),
                'time_range': [t_min, t_max],
                'freq_range': [f_min, f_max],
                'nfft': int(self.nfft_dropdown.currentText()),
            }

            # Save to NPZ
            np.savez(filepath,
                     metadata=np.array(metadata),
                     time=t_full,
                     frequency=f_sliced,  # Already in kHz, sliced
                     power=Sxx_log_sliced.T  # Transpose: (time, freq)
            )

            print(f"[Spectrogram] Data saved to: {filepath}")
            print(f"[Spectrogram]   Time: {len(t_full)} points, Freq: {len(f_sliced)} points")
            QMessageBox.information(self.frame, "Saved", f"Data saved to:\n{filepath}")

        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to save data: {str(e)}")

    def _set_signal_selection_state(self, enabled):
        """Enable or disable signal selection widgets"""
        for widget in self.signal_widgets:
            widget.setEnabled(enabled)

    def _load_shot_info(self):
        """Load shot info and populate channel list"""
        try:
            shot_number = int(self.shot_entry.text())
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
            return

        self.current_shot = shot_number

        # Update shot display label
        self.shot_display_label.setText(f'#{shot_number}')
        self.shot_display_label.setStyleSheet('color: white;')

        # Enable signal selection widgets
        self._set_signal_selection_state(True)

        # Load channels for current diagnostic
        diag_type = self.diag_dropdown.currentText()
        self._update_channel_list(diag_type, shot_number)

    def _update_channel_list(self, diag_type, shot_number):
        """Update channel list based on diagnostic type"""
        if diag_type == 'ECE':
            self._load_ece_channels(shot_number)
        elif diag_type == 'BES':
            self._load_bes_channels(shot_number)
        elif diag_type == 'TCI':
            self._load_tci_channels()
        elif diag_type.startswith('ECEI-'):
            device = diag_type.split('-')[1]  # GT, GR, or HT
            self._load_ecei_channels(shot_number, device)
        else:
            self._load_mirnov_channels(diag_type)

    def _load_ece_channels(self, shot_number):
        """Load ECE channel info for given shot using ECELoader"""
        if shot_number in self.ece_info:
            info = self.ece_info[shot_number]
            self._update_channel_dropdown_ece(info)
            return

        try:
            print(f"[Spectrogram] Loading ECE channel info for #{shot_number}...")

            # Use ECELoader to get channel positions
            positions = self.ece_loader.get_channel_positions(shot_number)

            # Cache the info
            self.ece_info[shot_number] = positions

            print(f"[Spectrogram]   I_TF = {positions['I_TF']:.2f} kA, {len(positions['channels'])} channels")

            self._update_channel_dropdown_ece(positions)

        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load ECE info: {str(e)}")

    def _update_channel_dropdown_ece(self, info):
        """Update channel dropdown with ECE channels (R, Z positions)"""
        channel_list = []
        for ch, R, Z in zip(info['channels'], info['R'], info['Z']):
            channel_list.append(f"ECE{ch:02d} (R={R:.3f}m, Z={Z:.3f}m)")

        self.channel_dropdown.clear()
        self.channel_dropdown.addItems(channel_list)
        if channel_list:
            self.channel_dropdown.setCurrentIndex(0)

    def _load_bes_channels(self, shot_number):
        """Load BES channel info for given shot using BESLoader"""
        if shot_number in self.bes_info:
            info = self.bes_info[shot_number]
            self._update_channel_dropdown_bes(info)
            return

        try:
            print(f"[Spectrogram] Loading BES channel info for #{shot_number}...")

            # Use BESLoader to get channel positions
            positions = self.bes_loader.get_channel_positions(shot_number)

            if not positions['channels']:
                QMessageBox.warning(self.frame, "Warning", "No BES channels available for this shot")
                return

            # Cache the info
            self.bes_info[shot_number] = positions

            print(f"[Spectrogram]   {len(positions['channels'])} BES channels loaded")

            self._update_channel_dropdown_bes(positions)

        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load BES info: {str(e)}")

    def _update_channel_dropdown_bes(self, info):
        """Update channel dropdown with BES channels (R, Z positions)"""
        channel_list = []
        for ch, R, Z in zip(info['channels'], info['R'], info['Z']):
            channel_list.append(f"BES_{ch} (R={R:.3f}m, Z={Z:.3f}m)")

        self.channel_dropdown.clear()
        self.channel_dropdown.addItems(channel_list)
        if channel_list:
            self.channel_dropdown.setCurrentIndex(0)

    def _load_tci_channels(self):
        """Load TCI channel list (processed ne and raw signals)"""
        # Processed ne channels
        channels = [f'TCI{i:02d} (ne)' for i in range(1, 6)]
        # Raw signal channels
        channels += [f'TCI{i:02d} (raw)' for i in range(1, 6)]

        self.channel_dropdown.clear()
        self.channel_dropdown.addItems(channels)
        if channels:
            self.channel_dropdown.setCurrentIndex(0)

    def _load_mirnov_channels(self, diag_type):
        """Load Mirnov channel list"""
        from config.diagnostic_config import DIAGNOSTICS
        mirnov_config = DIAGNOSTICS.get('Mirnov', {})

        if 'Toroidal' in diag_type:
            channels = mirnov_config.get('toroidal_channels', [])
        else:
            channels = mirnov_config.get('poloidal_channels', [])

        self.channel_dropdown.clear()
        self.channel_dropdown.addItems(channels)
        if channels:
            self.channel_dropdown.setCurrentIndex(0)

    def _load_ecei_channels(self, shot_number, device):
        """Load ECEI channel info for given shot and device using ECEILoader"""
        cache_key = (shot_number, device)
        if cache_key in self.ecei_info:
            info = self.ecei_info[cache_key]
            self._update_channel_dropdown_ecei(info, device)
            return

        try:
            print(f"[Spectrogram] Loading ECEI-{device} channel info for #{shot_number}...")

            # Use ECEILoader to get channel positions
            positions = self.ecei_loader.get_channel_positions(shot_number, device)

            # GR excluded vertical channels
            excluded_v = ECEILoader.GR_EXCLUDED_VERTICAL if device == 'GR' else []

            # Build channel list with positions (flattened from 2D arrays)
            channels = []
            R_values = []
            Z_values = []

            for v in range(1, 25):
                if v in excluded_v:
                    continue

                for r in range(1, 9):
                    ch_name = f'{device}{v:02d}{r:02d}'
                    R_pos = positions['R'][v - 1, r - 1]
                    Z_pos = positions['Z'][v - 1, r - 1]

                    channels.append(ch_name)
                    R_values.append(R_pos)
                    Z_values.append(Z_pos)

            # Cache the info
            self.ecei_info[cache_key] = {
                'channels': channels,
                'R': np.array(R_values),
                'Z': np.array(Z_values),
                'Bt': positions['Bt'],
                'mode': positions['mode'],
                'LO': positions['LO']
            }

            print(f"[Spectrogram]   Bt = {positions['Bt']:.1f} T, mode = {positions['mode']}, LO = {positions['LO']} GHz")
            print(f"[Spectrogram]   R range: {min(R_values):.3f} - {max(R_values):.3f} m")
            print(f"[Spectrogram]   {len(channels)} channels loaded")

            self._update_channel_dropdown_ecei(self.ecei_info[cache_key], device)

        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load ECEI-{device} info: {str(e)}")

    def _update_channel_dropdown_ecei(self, info, device):
        """Update channel dropdown with ECEI channels (R, Z)"""
        channel_list = []
        for ch, R, Z in zip(info['channels'], info['R'], info['Z']):
            channel_list.append(f"{ch} (R={R:.3f}m, Z={Z:.3f}m)")

        self.channel_dropdown.clear()
        self.channel_dropdown.addItems(channel_list)
        if channel_list:
            self.channel_dropdown.setCurrentIndex(0)

    def _on_diag_changed(self, text=None):
        """Handle diagnostic type change"""
        if self.current_shot is None:
            return

        diag_type = self.diag_dropdown.currentText()
        self._update_channel_list(diag_type, self.current_shot)

    def _format_frequency(self, freq_hz):
        """Format frequency with appropriate unit (Hz, kHz, MHz, GHz)"""
        if freq_hz >= 1e9:
            return f'{int(freq_hz / 1e9)}GHz'
        elif freq_hz >= 1e6:
            return f'{int(freq_hz / 1e6)}MHz'
        elif freq_hz >= 1e3:
            return f'{int(freq_hz / 1e3)}kHz'
        else:
            return f'{int(freq_hz)}Hz'

    def _on_dyn_range_changed(self, value):
        """Handle dynamic range slider change"""
        val = value / 10.0  # Convert from int slider to float
        self._dyn_range_value = val
        self.dyn_range_label.setText(f'{val:.1f}')

        # Update colorbar if spectrogram exists
        if self.im is not None and self.spectrogram_data is not None:
            vmax = self.spectrogram_data['vmax']
            new_vmin = vmax - val
            self.im.set_clim(vmin=new_vmin, vmax=vmax)
            self.canvas.draw_idle()

    def _load_signal_data(self):
        """Load selected signal data"""
        if self.current_shot is None:
            QMessageBox.critical(self.frame, "Error", "Please fetch a shot first")
            return None

        shot_number = self.current_shot
        channel_str = self.channel_dropdown.currentText()

        if not channel_str:
            QMessageBox.critical(self.frame, "Error", "Please select a channel")
            return None

        diag_type = self.diag_dropdown.currentText()

        try:
            mds = Connection(self.app_config.MDS_IP)
            mds.openTree('kstar', shot_number)

            if diag_type == 'ECE':
                # Parse ECE channel number: "ECE05 (R=2.105m, Z=0.000m)"
                ch_num = int(channel_str.split('(')[0].replace('ECE', '').strip())
                node_name = f'\\ECE{ch_num:02d}'
            elif diag_type == 'BES':
                # Parse BES channel: "BES_0108 (R=...)" or "BES_0108"
                if '(' in channel_str:
                    ch_name = channel_str.split('(')[0].replace('BES_', '').strip()
                else:
                    ch_name = channel_str.replace('BES_', '').strip()
                node_name = f'\\BES_{ch_name}:FOO'
            elif diag_type == 'TCI':
                # Parse TCI channel: "TCI01 (ne)" or "TCI01 (raw)"
                ch_num = int(channel_str[3:5])
                if '(ne)' in channel_str:
                    node_name = f'\\ne_tci{ch_num:02d}'
                else:  # raw
                    node_name = f'\\tci{ch_num:02d}_ls:foo'
            elif diag_type.startswith('ECEI-'):
                # Parse ECEI channel: "GT0101 (R=1.850m, Z=-0.150m)"
                device = diag_type.split('-')[1]
                ch_name = channel_str.split('(')[0].strip()
                node_name = f'\\ECEI_{ch_name}:FOO'
            else:
                # Mirnov channel name directly
                node_name = f'\\{channel_str}'

            print(f"[Spectrogram] Loading {node_name}...")
            data = mds.get(node_name).data()
            time = mds.get(f'dim_of({node_name})').data()

            mds.closeTree('kstar', shot_number)

            # TCI raw signal has one extra time frame at the beginning
            if diag_type == 'TCI' and '(raw)' in channel_str:
                time = time[1:]

            # Check for bad channel (mean == 0)
            if np.mean(data) == 0:
                QMessageBox.warning(self.frame, "Bad Channel",
                    f"{channel_str} appears to be a bad channel (mean=0).\n"
                    "Please select another channel.")
                return None

            # Build result dict
            result = {
                'time': np.array(time),
                'data': np.array(data),
                'shot': shot_number,
                'channel': channel_str,
                'diag_type': diag_type
            }

            # Add position info for ECEI
            if diag_type.startswith('ECEI-'):
                device = diag_type.split('-')[1]
                cache_key = (shot_number, device)
                if cache_key in self.ecei_info:
                    info = self.ecei_info[cache_key]
                    ch_name = channel_str.split('(')[0].strip()
                    try:
                        idx = info['channels'].index(ch_name)
                        result['R'] = info['R'][idx]
                        result['Z'] = info['Z'][idx]
                    except ValueError:
                        pass

            return result

        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load signal: {str(e)}")
            return None

    def _update_status(self, message, color='blue'):
        """Update status label with message and color"""
        self.status_label.setText(message)
        self.status_label.setStyleSheet(f'color: {color}; font-weight: bold;')
        QApplication.processEvents()

    def _plot_spectrogram(self):
        """Calculate and plot spectrogram"""
        t0 = tclock.time()

        # Update status
        self._update_status('Loading signal...', 'blue')

        # Load signal data
        signal = self._load_signal_data()
        if signal is None:
            self._update_status('Failed to load signal', 'red')
            return

        self.signal_data = signal

        # Get parameters
        try:
            t_min = float(self.time_min_entry.text())
            t_max = float(self.time_max_entry.text())
            f_min = float(self.freq_min_entry.text()) * 1e3  # kHz to Hz
            f_max = float(self.freq_max_entry.text()) * 1e3
            nfft = int(self.nfft_dropdown.currentText())
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Invalid parameter values")
            self._update_status('Error', 'red')
            return

        # Time slice
        time = signal['time']
        data = signal['data']

        mask = (time >= t_min) & (time <= t_max)
        if np.sum(mask) < nfft:
            QMessageBox.critical(self.frame, "Error", "Time range too short for selected NFFT")
            self._update_status('Error', 'red')
            return

        time_slice = time[mask]
        data_slice = data[mask]

        # Calculate sampling frequency
        fs = 1.0 / (time[1] - time[0])

        # Calculate spectrogram
        self._update_status('Calculating spectrogram...', 'blue')
        print(f"[Spectrogram] Calculating (NFFT={nfft}, fs={fs/1e6:.2f}MHz)...")
        f, t_spec, Sxx = spectrogram(data_slice, fs, nperseg=nfft)

        # Convert to log scale (a.u.)
        Sxx_log = np.log10(Sxx + 1e-20)  # Avoid log(0)
        vmax = np.max(Sxx_log)

        # Store spectrogram data
        self.spectrogram_data = {
            'f': f,
            't': t_spec + t_min,  # Offset to actual time
            'Sxx_log': Sxx_log,
            'vmax': vmax
        }

        # Clear and plot
        self._update_status('Plotting...', 'blue')
        self.ax.clear()

        # Frequency mask
        f_mask = (f >= f_min) & (f <= f_max)

        # Calculate vmin from dynamic range slider
        dyn_range = self._dyn_range_value
        vmin = vmax - dyn_range

        self.im = self.ax.imshow(
            Sxx_log[f_mask, :],
            aspect='auto',
            origin='lower',
            extent=[t_spec[0] + t_min, t_spec[-1] + t_min,
                    f[f_mask][0]/1e3, f[f_mask][-1]/1e3],
            vmin=vmin,
            vmax=vmax,
            cmap='viridis'
        )

        self.ax.set_xlabel('Time [s]')
        self.ax.set_ylabel('Frequency [kHz]')

        # Build title with position info for ECEI and sampling rate
        if signal.get('diag_type', '').startswith('ECEI-') and 'R' in signal:
            R = signal['R']
            Z = signal['Z']
            ch_name = signal['channel'].split('(')[0].strip()
            title = f"#{signal['shot']} ECEI_{ch_name} (R={R:.3f}m, Z={Z:.3f}m)"
        else:
            title = f"#{signal['shot']} {signal['channel']}"

        # Add sampling rate to title
        fs_str = self._format_frequency(fs)
        title += f" - {fs_str}"
        self.ax.set_title(title)

        # Update or create colorbar
        if self.colorbar is not None:
            self.colorbar.remove()
        self.colorbar = self.figure.colorbar(self.im, ax=self.ax, label='log$_{10}$(Power) [a.u.]')

        self.canvas.draw()

        elapsed = tclock.time() - t0
        self._update_status(f'Done ({elapsed:.1f}s)', 'green')
        print(f"[Spectrogram] Completed in {elapsed:.2f}s")

        # Enable save button
        self.save_button.setEnabled(True)

    def save_settings(self):
        """Save current tab settings"""
        settings = {
            "shot": self.shot_entry.text(),
            "tmin": self.time_min_entry.text(),
            "tmax": self.time_max_entry.text(),
            "fmin": self.freq_min_entry.text(),
            "fmax": self.freq_max_entry.text(),
            "nfft": self.nfft_dropdown.currentText(),
            "dynamic_range": str(self._dyn_range_value)
        }
        set_tab_settings("spectrogram", settings)

    def load_settings(self):
        """Load and apply saved settings"""
        settings = get_tab_settings("spectrogram")

        if settings.get("shot"):
            self.shot_entry.setText(settings["shot"])

        if settings.get("tmin"):
            self.time_min_entry.setText(settings["tmin"])

        if settings.get("tmax"):
            self.time_max_entry.setText(settings["tmax"])

        if settings.get("fmin"):
            self.freq_min_entry.setText(settings["fmin"])

        if settings.get("fmax"):
            self.freq_max_entry.setText(settings["fmax"])

        if settings.get("nfft"):
            self.nfft_dropdown.setCurrentText(settings["nfft"])

        if settings.get("dynamic_range"):
            val = float(settings["dynamic_range"])
            self._dyn_range_value = val
            self.dyn_range_slider.setValue(int(val * 10))
            self.dyn_range_label.setText(settings["dynamic_range"])


class PythonHighlighter(QSyntaxHighlighter):
    """Simple Python syntax highlighter for the example script dialog"""

    def __init__(self, document):
        super().__init__(document)
        import re
        self._rules = []

        # Keywords (purple)
        kw_format = QTextCharFormat()
        kw_format.setForeground(QColor('#c670e0'))
        keywords = r'\b(import|from|as|def|class|if|elif|else|for|while|try|except|with|return|yield|lambda|and|or|not|in|is|None|True|False|print)\b'
        self._rules.append((re.compile(keywords), kw_format))

        # Builtins (orange)
        builtin_format = QTextCharFormat()
        builtin_format.setForeground(QColor('#fab16c'))
        builtins = r'\b(np|plt|data|metadata|time|frequency|power|fig|ax|im)\b'
        self._rules.append((re.compile(builtins), builtin_format))

        # Numbers (yellow)
        num_format = QTextCharFormat()
        num_format.setForeground(QColor('#faed5c'))
        numbers = r'\b(\d+\.?\d*)\b'
        self._rules.append((re.compile(numbers), num_format))

        # Strings (green)
        str_format = QTextCharFormat()
        str_format.setForeground(QColor('#b0e686'))
        strings = r'(\"[^\"]*\"|\'[^\']*\'|f\"[^\"]*\"|f\'[^\']*\')'
        self._rules.append((re.compile(strings), str_format))

        # Comments (gray)
        comment_format = QTextCharFormat()
        comment_format.setForeground(QColor('#999999'))
        comments = r'(#[^\n]*)'
        self._rules.append((re.compile(comments), comment_format))

    def highlightBlock(self, text):
        for pattern, fmt in self._rules:
            for match in pattern.finditer(text):
                start = match.start()
                length = match.end() - match.start()
                self.setFormat(start, length, fmt)
