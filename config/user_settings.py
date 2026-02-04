#!/usr/bin/python3.8

"""
User settings management for PRISM
Handles loading, saving, and resetting user preferences
Settings are stored in ~/.config/prism/settings.json
"""

import os
import json
import tkinter as tk
from tkinter import messagebox

from config.app_config import VERSION, CONTACT_EMAIL


# Settings file path
SETTINGS_DIR = os.path.expanduser("~/.config/prism")
SETTINGS_FILE = os.path.join(SETTINGS_DIR, "settings.json")

# Default settings
DEFAULT_SETTINGS = {
    "last_seen_version": "0.0.0",
    "tabs": {
        "tivt_profile": {
            "shot": "",
            "efit_tree": "efitrt1 (RT for PCS)",
            "coord_type": "psi_N",
            "analysis_type": "nn"
        },
        "tivt_timetrace": {
            "shot": ""
        },
        "nete_profile": {
            "shot": "",
            "efit_tree": "efitrt1 (RT for PCS)",
            "coord_type": "psi_N",
            "diagnostic": "TS+ECE"
        },
        "nete_timetrace": {
            "shot": "",
            "diagnostic": "TS"
        },
        "mse_profile": {
            "shot": "",
            "efit_tree": "efitrt1 (RT for PCS)",
            "param": "q"
        },
        "mse_timetrace": {
            "shot": ""
        },
        "spectrogram": {
            "shot": "",
            "tmin": "0",
            "tmax": "10",
            "fmin": "0",
            "fmax": "100",
            "nfft": "1024",
            "dynamic_range": "6.0"
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
            "contour_levels": "50"
        },
        "tv": {
            "shot": ""
        },
        "irvb": {
            "shot": "",
            "efit_tree": "efitrt1 (RT for PCS)",
            "psi_bounds": "0.7, 1.0"
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
    """Show update notification popup"""
    if not should_show_update_popup():
        return
    
    # Build changelog message
    changelog, date = get_changelog_for_version(VERSION)
    
    # Create popup window
    popup = tk.Toplevel(parent)
    popup.title("What's New")
    popup.geometry("600x500")
    popup.resizable(True, True)
    popup.transient(parent)
    popup.grab_set()
    
    # Center the popup
    popup.update_idletasks()
    x = parent.winfo_x() + (parent.winfo_width() - popup.winfo_width()) // 2
    y = parent.winfo_y() + (parent.winfo_height() - popup.winfo_height()) // 2
    popup.geometry(f"+{x}+{y}")
    
    # Common font
    font_normal = ('TkDefaultFont', 10)
    
    # Version info label
    version_label = tk.Label(popup, text=f"Current version: PRISM v{VERSION} ({date})", font=font_normal)
    version_label.pack(pady=(10, 5))
    
    # Bottom frame for checkbox, contact, and button (pack first to reserve space)
    bottom_frame = tk.Frame(popup)
    bottom_frame.pack(side=tk.BOTTOM, fill='x', padx=10, pady=(5, 10))
    
    # Skip checkbox (left side)
    skip_var = tk.BooleanVar(value=False)
    skip_check = tk.Checkbutton(bottom_frame, text="Do not show again", 
                                 variable=skip_var, font=font_normal)
    skip_check.pack(side=tk.LEFT)
    
    # OK button (right side)
    def on_close():
        if skip_var.get():
            mark_version_seen()
            save_settings()
        popup.destroy()
    
    ok_button = tk.Button(bottom_frame, text="OK", command=on_close, width=10, font=font_normal)
    ok_button.pack(side=tk.RIGHT)
    
    # Contact info (center)
    contact_label = tk.Label(bottom_frame, text=f"Bug reports: {CONTACT_EMAIL}", font=font_normal)
    contact_label.pack(side=tk.RIGHT, padx=20)
    
    # Changelog text with scrollbar
    text_frame = tk.Frame(popup)
    text_frame.pack(fill='both', expand=True, padx=10, pady=5)
    
    scrollbar_y = tk.Scrollbar(text_frame, orient='vertical')
    scrollbar_y.pack(side=tk.RIGHT, fill='y')
    
    text_widget = tk.Text(text_frame, wrap='word', font=('TkDefaultFont', 10),
                          yscrollcommand=scrollbar_y.set,
                          bg='white', fg='black',
                          padx=10, pady=10,
                          spacing1=2, spacing3=2)
    
    # Render markdown to text widget
    _render_markdown(text_widget, changelog)
    
    text_widget.config(state='disabled')
    text_widget.pack(side=tk.LEFT, fill='both', expand=True)
    
    scrollbar_y.config(command=text_widget.yview)
    
    # Handle window close button
    popup.protocol("WM_DELETE_WINDOW", on_close)


def _render_markdown(text_widget, content):
    """Render markdown content to text widget (GitHub style)"""
    import re
    
    # Define tags
    text_widget.tag_configure('h2', font=('TkDefaultFont', 11, 'bold'), 
                              foreground='#0366d6', spacing1=15, spacing3=5)
    text_widget.tag_configure('h3', font=('TkDefaultFont', 12, 'bold'), 
                              foreground='#24292e', spacing1=10, spacing3=5)
    text_widget.tag_configure('bold', font=('TkDefaultFont', 10, 'bold'))
    text_widget.tag_configure('code', font=('Courier', 9), background='#f6f8fa', 
                              foreground='#24292e')
    text_widget.tag_configure('bullet', foreground='#6a737d')
    
    lines = content.split('\n')
    
    for line in lines:
        # Skip empty lines but add spacing
        if not line.strip():
            text_widget.insert('end', '\n')
            continue
        
        # Version headers [1.1.1] - 2026-01-08
        if re.match(r'^\[\d+\.\d+\.\d+\] - \d{4}-\d{2}-\d{2}', line.strip()):
            text_widget.insert('end', line + '\n', 'h2')
            continue
        
        # Headers (### Added -> Added)
        if line.startswith('###'):
            header_text = line.replace('###', '').strip()
            text_widget.insert('end', header_text + '\n', 'h3')
            continue
        
        # Process line with inline formatting
        # Replace bullet
        if line.strip().startswith('-'):
            indent = len(line) - len(line.lstrip())
            line = ' ' * indent + '•' + line.strip()[1:]
        
        # Parse inline elements (**bold** and `code`)
        pos = 0
        processed_line = line
        
        # Find all bold and code segments
        segments = []
        
        # Find **bold**
        for match in re.finditer(r'\*\*([^*]+)\*\*', line):
            segments.append((match.start(), match.end(), 'bold', match.group(1)))
        
        # Find `code`
        for match in re.finditer(r'`([^`]+)`', line):
            segments.append((match.start(), match.end(), 'code', match.group(1)))
        
        # Sort by position
        segments.sort(key=lambda x: x[0])
        
        if segments:
            # Build line with tags
            current_pos = 0
            for start, end, tag, text in segments:
                # Add text before this segment
                if start > current_pos:
                    text_widget.insert('end', line[current_pos:start])
                # Add tagged segment
                text_widget.insert('end', text, tag)
                current_pos = end
            # Add remaining text
            if current_pos < len(line):
                text_widget.insert('end', line[current_pos:])
            text_widget.insert('end', '\n')
        else:
            text_widget.insert('end', line + '\n')


def get_changelog_for_version(version):
    """Get changelog text for all patch versions in the same minor version from CHANGELOG.md
    
    For example, if version is 1.1.1, it returns changelog for all 1.1.x versions.
    """
    import re
    
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