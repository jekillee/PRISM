"""
User settings management for PRISM
Handles loading, saving, and resetting user preferences
Settings are stored in ~/.config/prism/settings.json
"""

import os
import json
import re

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QCheckBox, QTextEdit, QWidget, QTabWidget,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QPixmap

from config.app_config import VERSION, CONTACT_EMAIL


# Settings file path
SETTINGS_DIR = os.path.expanduser("~/.config/prism")
SETTINGS_FILE = os.path.join(SETTINGS_DIR, "settings.json")

# Default settings
DEFAULT_SETTINGS = {
    "last_seen_version": "0.0.0",
    "theme": "dark",
    "tabs": {
        "tivt_profile": {
            "shot": "",
            "efit_tree": "efitrt1 (RT for PCS)",
            "coord_type": "psi_N",
            "analysis_type": "nn",
            "color_mode": "Gradient(viridis)",
            "show_nodes": False,
            "label_fontsize": 12,
            "legend_fontsize": 8,
            "tick_fontsize": 10,
            "fit_funcs": {"Ti": "mtanh", "vT": "mtanh"}
        },
        "tivt_timetrace": {
            "shot": "",
            "color_mode": "Gradient(viridis)",
            "label_fontsize": 12,
            "legend_fontsize": 8,
            "tick_fontsize": 10
        },
        "nete_profile": {
            "shot": "",
            "efit_tree": "efitrt1 (RT for PCS)",
            "coord_type": "psi_N",
            "diagnostic": "TS+ECE",
            "color_mode": "Gradient(viridis)",
            "show_nodes": False,
            "label_fontsize": 12,
            "legend_fontsize": 8,
            "tick_fontsize": 10,
            "fit_funcs": {"Te": "mtanh", "ne": "mtanh"},
            "tci_validation": False
        },
        "nete_timetrace": {
            "shot": "",
            "diagnostic": "TS",
            "color_mode": "Gradient(viridis)",
            "label_fontsize": 12,
            "legend_fontsize": 8,
            "tick_fontsize": 10
        },
        "mse_profile": {
            "shot": "",
            "efit_tree": "efitrt1 (RT for PCS)",
            "param": "q",
            "color_mode": "Gradient(viridis)",
            "show_nodes": False,
            "label_fontsize": 12,
            "legend_fontsize": 8,
            "tick_fontsize": 10
        },
        "mse_timetrace": {
            "shot": "",
            "color_mode": "Gradient(viridis)",
            "label_fontsize": 12,
            "legend_fontsize": 8,
            "tick_fontsize": 10
        },
        "spectrogram": {
            "shot": "",
            "tmin": "0",
            "tmax": "10",
            "fmin": "0",
            "fmax": "100",
            "nfft": "1024",
            "dynamic_range": "6.0",
            "colormap": "viridis",
            "label_fontsize": 12,
            "tick_fontsize": 10,
            "title_fontsize": 12
        },
        "nmode": {
            "shot": "",
            "tmin": "0.0",
            "tmax": "10.0",
            "tinterval": "0.01",
            "fmin": "0",
            "fmax": "100",
            "nmodes": 5,
            "tolerance": "0.8",
            "fraction": "0.01",
            "sign": 1,
            "integrate": False,
            "detrend": True,
            "plot_type": "contour",
            "color_mode": "Default",
            "label_fontsize": 12,
            "legend_fontsize": 8,
            "tick_fontsize": 10,
            "title_fontsize": 12,
            "contour_levels": "50",
            "contour_linewidth": 0.8,
            "amp_linewidth": 1.5,
            "selected_modes": [1, 2, 3, 4, 5]
        },
        "tv": {
            "shot": ""
        },
        "neutron": {
            "shot": "",
            "color_mode": "Gradient(viridis)",
            "label_fontsize": 12,
            "legend_fontsize": 8,
            "tick_fontsize": 10
        },
        "irvb": {
            "shot": "",
            "efit_tree": "efitrt1 (RT for PCS)",
            "psi_bounds": "0.7, 1.0",
            "trace_color_mode": "Fixed(tab10)",
            "trace_label_fontsize": 12,
            "trace_legend_fontsize": 8,
            "trace_tick_fontsize": 10,
            "plot2d_colormap": "viridis",
            "plot2d_lcfs_color": "black",
            "plot2d_maxis_color": "black",
            "plot2d_limiter_color": "white",
            "plot2d_flux_color": "gray",
            "plot2d_label_fontsize": 12,
            "plot2d_tick_fontsize": 10
        }
    }
}

# Global settings cache
_settings = None


def _ensure_settings_dir():
    """Create settings directory if it doesn't exist"""
    if not os.path.exists(SETTINGS_DIR):
        os.makedirs(SETTINGS_DIR)


def load_settings():
    """Load settings from file, create default if not exists"""
    global _settings

    _ensure_settings_dir()

    if os.path.exists(SETTINGS_FILE):
        try:
            with open(SETTINGS_FILE, 'r') as f:
                _settings = json.load(f)

            # Merge with defaults for any missing keys
            _settings = _merge_defaults(_settings, DEFAULT_SETTINGS)

        except (json.JSONDecodeError, IOError) as e:
            print(f"[Settings] Warning: Failed to load settings, using defaults: {e}")
            _settings = DEFAULT_SETTINGS.copy()
    else:
        _settings = DEFAULT_SETTINGS.copy()

    return _settings


def _merge_defaults(settings, defaults):
    """Recursively merge defaults into settings for missing keys"""
    result = defaults.copy()

    for key, value in settings.items():
        if key in result:
            if isinstance(value, dict) and isinstance(result[key], dict):
                result[key] = _merge_defaults(value, result[key])
            else:
                result[key] = value
        else:
            result[key] = value

    return result


def save_settings():
    """Save current settings to file"""
    global _settings

    if _settings is None:
        return

    _ensure_settings_dir()

    try:
        with open(SETTINGS_FILE, 'w') as f:
            json.dump(_settings, f, indent=4)
        print(f"[Settings] Saved to {SETTINGS_FILE}")
    except IOError as e:
        print(f"[Settings] Warning: Failed to save settings: {e}")


def get_settings():
    """Get current settings (load if not loaded)"""
    global _settings

    if _settings is None:
        load_settings()

    return _settings


def get_tab_settings(tab_name):
    """Get settings for specific tab"""
    settings = get_settings()
    return settings.get("tabs", {}).get(tab_name, {})


def set_tab_settings(tab_name, tab_settings):
    """Set settings for specific tab"""
    global _settings

    if _settings is None:
        load_settings()

    if "tabs" not in _settings:
        _settings["tabs"] = {}

    _settings["tabs"][tab_name] = tab_settings


def get_theme():
    """Get saved theme setting"""
    settings = get_settings()
    return settings.get("theme", "dark")


def set_theme(theme_name):
    """Set theme in settings (call save_settings() to persist)"""
    global _settings

    if _settings is None:
        load_settings()

    _settings["theme"] = theme_name


def reset_settings():
    """Reset all settings to defaults"""
    global _settings
    _settings = DEFAULT_SETTINGS.copy()
    save_settings()
    print("Settings reset to defaults")


def should_show_update_popup():
    """Check if update popup should be shown"""
    settings = get_settings()
    last_seen = settings.get("last_seen_version", "0.0.0")

    # Show popup if version changed
    return VERSION != last_seen


def mark_version_seen():
    """Mark current version as seen (only called when skip checkbox is checked)"""
    global _settings

    if _settings is None:
        load_settings()

    _settings["last_seen_version"] = VERSION


def show_update_popup(parent):
    """Show update notification popup with tabbed changelog (v2, v1)"""
    if not should_show_update_popup():
        return

    # Build changelog for header date display
    _, date = get_changelog_for_version(VERSION)

    # Create popup dialog
    popup = QDialog(parent)
    popup.setWindowTitle("What's New")
    popup.resize(600, 550)
    popup.setWindowModality(Qt.NonModal)

    # Main layout
    layout = QVBoxLayout(popup)
    layout.setContentsMargins(10, 10, 10, 10)

    # Common font
    font_normal = QFont()
    font_normal.setPointSize(10)

    # Logo + version header
    header_widget = QWidget()
    header_layout = QHBoxLayout(header_widget)
    header_layout.setContentsMargins(0, 0, 0, 0)
    header_layout.addStretch()

    icon_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ui', 'icons')
    from ui.theme import ThemeManager
    logo_file = 'prism-logo-dark.svg' if ThemeManager.current_theme == 'dark' else 'prism-logo-light.svg'
    logo_path = os.path.join(icon_dir, logo_file)
    pixmap = QPixmap(logo_path)
    if not pixmap.isNull():
        logo_label = QLabel()
        pixmap = pixmap.scaled(36, 36, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        logo_label.setPixmap(pixmap)
        header_layout.addWidget(logo_label)

    version_label = QLabel(
        '<span style="font-weight:bold;">'
        '<span style="color:#ff4444;">P</span>'
        '<span style="color:#ff8c00;">R</span>'
        '<span style="color:#ffd700;">I</span>'
        '<span style="color:#22c55e;">S</span>'
        '<span style="color:#3b82f6;">M</span>'
        '</span>'
        f' v{VERSION} ({date})'
    )
    version_label.setFont(font_normal)
    header_layout.addWidget(version_label)
    header_layout.addStretch()
    layout.addWidget(header_widget)

    # Tabbed changelog — group by major version
    tab_widget = QTabWidget()
    tab_widget.setFont(font_normal)

    all_changelogs = _get_all_major_changelogs()
    current_major = int(VERSION.split('.')[0])

    for major_ver, cl_text in sorted(all_changelogs.items(), reverse=True):
        text_widget = QTextEdit()
        text_widget.setReadOnly(True)
        text_widget.setFont(font_normal)
        text_widget.setHtml(_render_markdown(cl_text))
        tab_widget.addTab(text_widget, f"v{major_ver}")

    layout.addWidget(tab_widget)

    # Bottom widget for checkbox, contact, and button
    bottom_widget = QWidget()
    bottom_layout = QHBoxLayout(bottom_widget)
    bottom_layout.setContentsMargins(0, 5, 0, 0)

    # Skip checkbox (left side)
    skip_check = QCheckBox("Do not show again")
    skip_check.setFont(font_normal)
    bottom_layout.addWidget(skip_check)

    # Spacer
    bottom_layout.addStretch()

    # Contact info (center-right)
    contact_label = QLabel(f"Bug reports: {CONTACT_EMAIL}")
    contact_label.setFont(font_normal)
    bottom_layout.addWidget(contact_label)

    # OK button (right side)
    def on_close():
        if skip_check.isChecked():
            mark_version_seen()
            save_settings()
        popup.accept()

    ok_button = QPushButton("OK")
    ok_button.setFont(font_normal)
    ok_button.setFixedWidth(80)
    ok_button.clicked.connect(on_close)
    bottom_layout.addWidget(ok_button)

    layout.addWidget(bottom_widget)

    # Handle window close button (X) via finished signal
    popup.finished.connect(lambda result: on_close() if result == 0 else None)

    parent._update_popup = popup  # prevent GC
    popup.show()


def _render_markdown(content):
    """Render markdown content to HTML string (GitHub style)"""
    lines = content.split('\n')
    html_parts = []
    in_list = 0  # nesting level: 0=none, 1=top, 2=sub

    for line in lines:
        stripped = line.strip()

        # Skip empty lines
        if not stripped:
            while in_list > 0:
                html_parts.append('</ul>')
                in_list -= 1
            continue

        # Version headers [1.1.1] - 2026-01-08
        if re.match(r'^\[\d+\.\d+\.\d+\] - \d{4}-\d{2}-\d{2}', stripped):
            while in_list > 0:
                html_parts.append('</ul>')
                in_list -= 1
            html_parts.append(f'<h2 style="color: #0d6efd; margin-top: 15px; margin-bottom: 5px;">{_escape_html(stripped)}</h2>')
            continue

        # Headers (### Added -> Added)
        if stripped.startswith('###'):
            while in_list > 0:
                html_parts.append('</ul>')
                in_list -= 1
            header_text = stripped.replace('###', '').strip()
            html_parts.append(f'<h3 style="margin-top: 10px; margin-bottom: 5px;">{_escape_html(header_text)}</h3>')
            continue

        # Sub-bullet (indented: "  - text")
        if line.startswith('  -') or line.startswith('    -'):
            if in_list == 0:
                html_parts.append('<ul style="margin: 2px 0;">')
                in_list = 1
            if in_list == 1:
                html_parts.append('<ul style="margin: 1px 0 1px 15px;">')
                in_list = 2
            item_text = stripped[1:].strip()
            html_parts.append(f'<li style="margin: 1px 0;">{_format_inline(item_text)}</li>')
            continue

        # Top-level bullet
        if stripped.startswith('-'):
            if in_list == 2:
                html_parts.append('</ul>')
                in_list = 1
            if in_list == 0:
                html_parts.append('<ul style="margin: 2px 0;">')
                in_list = 1
            item_text = stripped[1:].strip()
            html_parts.append(f'<li style="margin: 2px 0;">{_format_inline(item_text)}</li>')
            continue

        # Regular text
        while in_list > 0:
            html_parts.append('</ul>')
            in_list -= 1
        formatted = _format_inline(_escape_html(line))
        html_parts.append(f'<p style="margin: 2px 0;">{formatted}</p>')

    while in_list > 0:
        html_parts.append('</ul>')
        in_list -= 1

    return '\n'.join(html_parts)


def _escape_html(text):
    """Escape HTML special characters"""
    return text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def _format_inline(text):
    """Apply inline formatting: **bold** and `code`"""
    # Replace **bold**
    text = re.sub(r'\*\*([^*]+)\*\*', r'<b>\1</b>', text)
    # Replace `code`
    text = re.sub(r'`([^`]+)`', r'<code style="font-family: Courier;">\1</code>', text)
    return text


def _get_all_major_changelogs():
    """Get changelogs grouped by major version number.

    Returns:
        dict: {major_version_int: changelog_text}
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    changelog_path = os.path.join(script_dir, '..', 'CHANGELOG.md')

    if not os.path.exists(changelog_path):
        changelog_path = os.path.join(script_dir, 'CHANGELOG.md')

    if not os.path.exists(changelog_path):
        return {}

    try:
        with open(changelog_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # Find all version entries: ## [X.Y.Z] - YYYY-MM-DD
        pattern = r'## \[(\d+)\.(\d+\.\d+)\] - (\d{4}-\d{2}-\d{2})\s*\n(.*?)(?=\n## \[|$)'
        matches = list(re.finditer(pattern, content, re.DOTALL))

        result = {}
        for match in matches:
            major = int(match.group(1))
            minor_patch = match.group(2)
            date = match.group(3)
            changelog_text = match.group(4).strip()
            ver = f"{major}.{minor_patch}"

            entry = f"[{ver}] - {date}\n{changelog_text}"
            if major not in result:
                result[major] = []
            result[major].append(entry)

        return {k: "\n\n".join(v) for k, v in result.items()}

    except Exception as e:
        print(f"[Settings] Warning: Could not read CHANGELOG.md: {e}")
        return {}


def get_changelog_for_version(version):
    """Get changelog text for all patch versions in the same minor version from CHANGELOG.md

    For example, if version is 1.1.1, it returns changelog for all 1.1.x versions.
    """
    # Find CHANGELOG.md relative to this file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    changelog_path = os.path.join(script_dir, '..', 'CHANGELOG.md')

    if not os.path.exists(changelog_path):
        # Try alternate location
        changelog_path = os.path.join(script_dir, 'CHANGELOG.md')

    if not os.path.exists(changelog_path):
        return ("Bug fixes and improvements.", "")

    try:
        with open(changelog_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # Parse version to get major.minor prefix
        version_parts = version.split('.')
        if len(version_parts) >= 2:
            minor_prefix = f"{version_parts[0]}.{version_parts[1]}"
        else:
            minor_prefix = version

        # Find all versions matching the minor prefix (e.g., 1.1.*)
        # Pattern: ## [1.1.0] - 2026-01-07
        pattern = rf'## \[({re.escape(minor_prefix)}\.\d+)\] - (\d{{4}}-\d{{2}}-\d{{2}})\s*\n(.*?)(?=\n## \[|$)'
        matches = list(re.finditer(pattern, content, re.DOTALL))

        if not matches:
            return ("Bug fixes and improvements.", "")

        # Build combined changelog (newest first - matches are already in file order)
        combined_changelog = []
        latest_date = ""

        for match in matches:
            ver = match.group(1)
            date = match.group(2)
            changelog_text = match.group(3).strip()

            if not latest_date:
                latest_date = date

            # Add version header
            combined_changelog.append(f"[{ver}] - {date}\n{changelog_text}")

        # Join with separator
        full_changelog = "\n\n".join(combined_changelog)

        return (full_changelog, latest_date)

    except Exception as e:
        print(f"[Settings] Warning: Could not read CHANGELOG.md: {e}")
        return ("Bug fixes and improvements.", "")
