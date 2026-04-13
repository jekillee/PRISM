"""
Theme management for PRISM - dark/light theme switching

Usage:
    from ui.theme import ThemeManager
    ThemeManager.apply_theme('dark')   # or 'light'
    ThemeManager.toggle_theme()
    current = ThemeManager.current_theme
"""

import os
import matplotlib as mpl
from PySide6.QtWidgets import QApplication, QScrollArea
from PySide6.QtGui import QPalette, QColor

# Icon paths
_ICON_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'icons')
_CHECK_ICON = os.path.join(_ICON_DIR, 'check-white.svg').replace('\\', '/')
_ARROW_DOWN_LIGHT = os.path.join(_ICON_DIR, 'arrow-down-light.svg').replace('\\', '/')
_ARROW_DOWN_DARK = os.path.join(_ICON_DIR, 'arrow-down-dark.svg').replace('\\', '/')


# ============================================================
# Dark Theme
# ============================================================

_DARK_QSS = """
QMainWindow { background: #2b2b2b; }

/* Matplotlib navigation toolbar */
NavigationToolbar2QT {
    background: #2b2b2b;
    border: none;
    spacing: 2px;
}
NavigationToolbar2QT QToolButton {
    background: transparent;
    border: 1px solid transparent;
    border-radius: 3px;
    padding: 2px;
    margin: 1px;
}
NavigationToolbar2QT QToolButton:hover {
    background: #3c3f41;
    border-color: #555;
}
NavigationToolbar2QT QToolButton:checked {
    background: #0d6efd;
    border-color: #0d6efd;
}
NavigationToolbar2QT QToolButton:pressed {
    background: #0a58ca;
}

/* Sidebar */
QTreeWidget#sidebar {
    background: #1e1e1e;
    color: #cccccc;
    border: none;
    font-size: 13px;
    outline: none;
}
QTreeWidget#sidebar::item {
    padding: 6px 8px;
    border-radius: 4px;
    margin: 1px 4px;
}
QTreeWidget#sidebar::item:selected {
    background: #37373d;
    color: #ffffff;
}
QTreeWidget#sidebar::item:hover:!selected {
    background: #2a2d2e;
}

/* Splitter handle */
QSplitter::handle:horizontal {
    background: #3c3c3c;
    width: 1px;
}

/* Controls */
QGroupBox {
    color: #cccccc;
    border: 1px solid #444;
    border-radius: 6px;
    margin-top: 12px;
    padding-top: 16px;
    font-weight: bold;
    font-size: 12px;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 12px;
    padding: 0 6px;
}
QPushButton {
    background: #0d6efd;
    color: white;
    border: none;
    border-radius: 4px;
    padding: 7px 16px;
    font-size: 12px;
    font-weight: bold;
}
QPushButton:hover { background: #3d8bfd; }
QPushButton:pressed { background: #0a58ca; }
QPushButton:disabled { background: #444; color: #777; }

QComboBox {
    background: #3c3f41;
    color: #ccc;
    border: 1px solid #555;
    border-radius: 4px;
    padding: 4px 6px;
    padding-right: 0px;
    min-height: 20px;
}
QComboBox::drop-down {
    border: none;
    width: 22px;
    subcontrol-origin: padding;
    subcontrol-position: center right;
}
QComboBox::down-arrow {
    image: url(""" + _ARROW_DOWN_LIGHT + """);
    width: 10px;
    height: 10px;
}
QComboBox QAbstractItemView {
    background: #2b2b2b;
    color: #ccc;
    selection-background-color: #37373d;
}
QComboBox:disabled {
    background: #2b2b2b;
    color: #555;
    border: 1px solid #3c3c3c;
}

QLineEdit {
    background: #3c3f41;
    color: #ccc;
    border: 1px solid #555;
    border-radius: 4px;
    padding: 4px 8px;
}
QLineEdit:disabled {
    background: #2b2b2b;
    color: #555;
    border: 1px solid #3c3c3c;
}

QCheckBox { color: #ccc; spacing: 6px; }
QCheckBox::indicator {
    width: 16px; height: 16px;
    border: 1px solid #666;
    border-radius: 3px;
    background: #3c3f41;
}
QCheckBox::indicator:checked {
    background: #0d6efd;
    border-color: #0d6efd;
    image: url(""" + _CHECK_ICON + """);
}

QRadioButton { color: #ccc; spacing: 8px; }
QRadioButton:disabled { color: #555; }
QRadioButton::indicator {
    width: 16px; height: 16px;
    border: 1px solid #666;
    border-radius: 3px;
    background: #3c3f41;
}
QRadioButton::indicator:checked {
    background: #0d6efd;
    border-color: #0d6efd;
}
QRadioButton::indicator:disabled {
    background: #333;
    border-color: #444;
}

QLabel { color: #cccccc; }

QListWidget {
    background: #1e1e1e;
    color: #ccc;
    border: 1px solid #444;
    border-radius: 4px;
    font-size: 12px;
}
QListWidget::item { padding: 3px 6px; }
QListWidget::item:selected { background: #37373d; color: white; }

QSlider::groove:horizontal {
    height: 6px;
    background: #555;
    border-radius: 3px;
}
QSlider::handle:horizontal {
    width: 16px; height: 16px;
    margin: -5px 0;
    background: #0d6efd;
    border-radius: 8px;
}

QStatusBar { background: #1e1e1e; color: #888; font-size: 11px; }
QProgressBar {
    background: #3c3f41;
    border: none;
    border-radius: 4px;
    text-align: center;
    color: white;
    font-size: 10px;
    max-height: 14px;
}
QProgressBar::chunk { background: #0d6efd; border-radius: 4px; }

QScrollBar:vertical {
    background: #2b2b2b;
    width: 12px;
    margin: 0;
}
QScrollBar::handle:vertical {
    background: #555;
    border-radius: 4px;
    min-height: 20px;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }

QScrollBar:horizontal {
    background: #2b2b2b;
    height: 12px;
    margin: 0;
}
QScrollBar::handle:horizontal {
    background: #555;
    border-radius: 4px;
    min-width: 20px;
}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0; }

QDialog { background: #2b2b2b; }

QScrollArea { border: none; background: transparent; }
QScrollArea > QWidget { background: transparent; }

QTextEdit {
    background: #1e1e1e;
    color: #ccc;
    border: 1px solid #444;
    border-radius: 4px;
}
"""

_DARK_PALETTE = {
    'Window':          '#2b2b2b',
    'WindowText':      '#cccccc',
    'Base':            '#1e1e1e',
    'Text':            '#cccccc',
    'Button':          '#3c3f41',
    'ButtonText':      '#cccccc',
    'Highlight':       '#0d6efd',
    'HighlightedText': '#ffffff',
}

_DARK_MPL = {
    'figure.facecolor':  '#2b2b2b',
    'axes.facecolor':    '#1e1e1e',
    'axes.edgecolor':    '#555555',
    'axes.labelcolor':   '#cccccc',
    'xtick.color':       '#cccccc',
    'ytick.color':       '#cccccc',
    'text.color':        '#cccccc',
    'grid.color':        '#444444',
    'legend.facecolor':  '#2b2b2b',
    'legend.edgecolor':  '#555555',
    'legend.labelcolor': '#cccccc',
}

# ============================================================
# Light Theme
# ============================================================

_LIGHT_QSS = """
QMainWindow { background: #f0f0f0; }

/* Matplotlib navigation toolbar */
NavigationToolbar2QT {
    background: #f0f0f0;
    border: none;
    spacing: 2px;
}
NavigationToolbar2QT QToolButton {
    background: transparent;
    border: 1px solid transparent;
    border-radius: 3px;
    padding: 2px;
    margin: 1px;
}
NavigationToolbar2QT QToolButton:hover {
    background: #d0d0d0;
    border-color: #bbb;
}
NavigationToolbar2QT QToolButton:checked {
    background: #0d6efd;
    border-color: #0d6efd;
}
NavigationToolbar2QT QToolButton:pressed {
    background: #0a58ca;
}

/* Sidebar */
QTreeWidget#sidebar {
    background: #e8e8e8;
    color: #333333;
    border: none;
    font-size: 13px;
    outline: none;
}
QTreeWidget#sidebar::item {
    padding: 6px 8px;
    border-radius: 4px;
    margin: 1px 4px;
}
QTreeWidget#sidebar::item:selected {
    background: #d0d0d0;
    color: #000000;
}
QTreeWidget#sidebar::item:hover:!selected {
    background: #dcdcdc;
}

/* Splitter handle */
QSplitter::handle:horizontal {
    background: #cccccc;
    width: 1px;
}

/* Controls */
QGroupBox {
    color: #333333;
    border: 1px solid #cccccc;
    border-radius: 6px;
    margin-top: 12px;
    padding-top: 16px;
    font-weight: bold;
    font-size: 12px;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 12px;
    padding: 0 6px;
}
QPushButton {
    background: #0d6efd;
    color: white;
    border: none;
    border-radius: 4px;
    padding: 7px 16px;
    font-size: 12px;
    font-weight: bold;
}
QPushButton:hover { background: #3d8bfd; }
QPushButton:pressed { background: #0a58ca; }
QPushButton:disabled { background: #cccccc; color: #888; }

QComboBox {
    background: #ffffff;
    color: #333;
    border: 1px solid #cccccc;
    border-radius: 4px;
    padding: 4px 6px;
    padding-right: 0px;
    min-height: 20px;
}
QComboBox::drop-down {
    border: none;
    width: 22px;
    subcontrol-origin: padding;
    subcontrol-position: center right;
}
QComboBox::down-arrow {
    image: url(""" + _ARROW_DOWN_DARK + """);
    width: 10px;
    height: 10px;
}
QComboBox QAbstractItemView {
    background: #ffffff;
    color: #333;
    selection-background-color: #d0d0d0;
}
QComboBox:disabled {
    background: #e8e8e8;
    color: #aaa;
    border: 1px solid #ddd;
}

QLineEdit {
    background: #ffffff;
    color: #333;
    border: 1px solid #cccccc;
    border-radius: 4px;
    padding: 4px 8px;
}
QLineEdit:disabled {
    background: #e8e8e8;
    color: #aaa;
    border: 1px solid #ddd;
}

QCheckBox { color: #333; spacing: 6px; }
QCheckBox::indicator {
    width: 16px; height: 16px;
    border: 1px solid #aaa;
    border-radius: 3px;
    background: #ffffff;
}
QCheckBox::indicator:checked {
    background: #0d6efd;
    border-color: #0d6efd;
    image: url(""" + _CHECK_ICON + """);
}

QRadioButton { color: #333; spacing: 8px; }
QRadioButton:disabled { color: #aaa; }
QRadioButton::indicator {
    width: 16px; height: 16px;
    border: 1px solid #aaa;
    border-radius: 3px;
    background: #ffffff;
}
QRadioButton::indicator:checked {
    background: #0d6efd;
    border-color: #0d6efd;
}
QRadioButton::indicator:disabled {
    background: #e0e0e0;
    border-color: #ccc;
}

QLabel { color: #333333; }

QListWidget {
    background: #ffffff;
    color: #333;
    border: 1px solid #cccccc;
    border-radius: 4px;
    font-size: 12px;
}
QListWidget::item { padding: 3px 6px; }
QListWidget::item:selected { background: #d0d0d0; color: black; }

QSlider::groove:horizontal {
    height: 6px;
    background: #cccccc;
    border-radius: 3px;
}
QSlider::handle:horizontal {
    width: 16px; height: 16px;
    margin: -5px 0;
    background: #0d6efd;
    border-radius: 8px;
}

QStatusBar { background: #e8e8e8; color: #666; font-size: 11px; }
QProgressBar {
    background: #dddddd;
    border: none;
    border-radius: 4px;
    text-align: center;
    color: #333;
    font-size: 10px;
    max-height: 14px;
}
QProgressBar::chunk { background: #0d6efd; border-radius: 4px; }

QScrollBar:vertical {
    background: #f0f0f0;
    width: 12px;
    margin: 0;
}
QScrollBar::handle:vertical {
    background: #cccccc;
    border-radius: 4px;
    min-height: 20px;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }

QScrollBar:horizontal {
    background: #f0f0f0;
    height: 12px;
    margin: 0;
}
QScrollBar::handle:horizontal {
    background: #cccccc;
    border-radius: 4px;
    min-width: 20px;
}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0; }

QDialog { background: #f0f0f0; }

QScrollArea { border: none; background: transparent; }
QScrollArea > QWidget { background: transparent; }

QTextEdit {
    background: #ffffff;
    color: #333;
    border: 1px solid #cccccc;
    border-radius: 4px;
}
"""

_LIGHT_PALETTE = {
    'Window':          '#f0f0f0',
    'WindowText':      '#333333',
    'Base':            '#ffffff',
    'Text':            '#333333',
    'Button':          '#e0e0e0',
    'ButtonText':      '#333333',
    'Highlight':       '#0d6efd',
    'HighlightedText': '#ffffff',
}

_LIGHT_MPL = {
    'figure.facecolor':  '#f0f0f0',
    'axes.facecolor':    '#ffffff',
    'axes.edgecolor':    '#333333',
    'axes.labelcolor':   '#333333',
    'xtick.color':       '#333333',
    'ytick.color':       '#333333',
    'text.color':        '#333333',
    'grid.color':        '#cccccc',
    'legend.facecolor':  '#ffffff',
    'legend.edgecolor':  '#cccccc',
    'legend.labelcolor': '#333333',
}

# ============================================================
# Theme definitions
# ============================================================

THEMES = {
    'dark': {
        'qss': _DARK_QSS,
        'palette': _DARK_PALETTE,
        'mpl': _DARK_MPL,
        'sidebar_category_color': '#888888',
        'title_color': '#0d6efd',
        'dev_label_color': '#888',
        'welcome_color': '#555',
    },
    'light': {
        'qss': _LIGHT_QSS,
        'palette': _LIGHT_PALETTE,
        'mpl': _LIGHT_MPL,
        'sidebar_category_color': '#888888',
        'title_color': '#0d6efd',
        'dev_label_color': '#999',
        'welcome_color': '#aaa',
    },
}


class ThemeManager:
    """Singleton theme manager for PRISM application"""

    current_theme = 'dark'
    _callbacks = []

    @classmethod
    def get_theme(cls):
        """Get the current theme definition"""
        return THEMES[cls.current_theme]

    @classmethod
    def apply_theme(cls, theme_name, save=False):
        """Apply a theme globally (QSS + palette + matplotlib rcParams)

        Args:
            theme_name: 'dark' or 'light'
            save: If True, persist the choice to settings.json
        """
        if theme_name not in THEMES:
            return

        cls.current_theme = theme_name

        # Persist theme choice to settings
        if save:
            try:
                from config.user_settings import set_theme, save_settings
                set_theme(theme_name)
                save_settings()
            except Exception:
                pass
        theme = THEMES[theme_name]

        app = QApplication.instance()
        if app is None:
            return

        # Apply QSS
        app.setStyleSheet(theme['qss'])

        # Apply palette
        palette = app.palette()
        pal = theme['palette']
        palette.setColor(QPalette.Window, QColor(pal['Window']))
        palette.setColor(QPalette.WindowText, QColor(pal['WindowText']))
        palette.setColor(QPalette.Base, QColor(pal['Base']))
        palette.setColor(QPalette.Text, QColor(pal['Text']))
        palette.setColor(QPalette.Button, QColor(pal['Button']))
        palette.setColor(QPalette.ButtonText, QColor(pal['ButtonText']))
        palette.setColor(QPalette.Highlight, QColor(pal['Highlight']))
        palette.setColor(QPalette.HighlightedText, QColor(pal['HighlightedText']))
        app.setPalette(palette)

        # Fix QScrollArea viewport backgrounds for correct theme rendering
        for widget in app.allWidgets():
            if isinstance(widget, QScrollArea):
                widget.viewport().setAutoFillBackground(False)

        # Apply matplotlib rcParams
        mpl.rcParams.update(theme['mpl'])

        # Notify registered callbacks (e.g., to refresh canvases)
        for callback in cls._callbacks:
            try:
                callback(theme_name)
            except Exception:
                pass

    @classmethod
    def toggle_theme(cls):
        """Toggle between dark and light theme (auto-saves to settings)"""
        new_theme = 'light' if cls.current_theme == 'dark' else 'dark'
        cls.apply_theme(new_theme, save=True)

    @classmethod
    def on_theme_changed(cls, callback):
        """Register a callback for theme changes: callback(theme_name)"""
        cls._callbacks.append(callback)

    @classmethod
    def off_theme_changed(cls, callback):
        """Unregister a callback for theme changes"""
        try:
            cls._callbacks.remove(callback)
        except ValueError:
            pass

    @classmethod
    def apply_theme_to_figure(cls, figure):
        """Apply current theme colors to an existing matplotlib figure"""
        theme = THEMES[cls.current_theme]
        mpl_colors = theme['mpl']
        zero_line_color = 'white' if cls.current_theme == 'dark' else 'gray'

        figure.set_facecolor(mpl_colors['figure.facecolor'])
        for ax in figure.get_axes():
            ax.set_facecolor(mpl_colors['axes.facecolor'])
            ax.tick_params(colors=mpl_colors['xtick.color'])
            ax.xaxis.label.set_color(mpl_colors['axes.labelcolor'])
            ax.yaxis.label.set_color(mpl_colors['axes.labelcolor'])
            ax.title.set_color(mpl_colors['text.color'])
            for spine in ax.spines.values():
                spine.set_edgecolor(mpl_colors['axes.edgecolor'])
            # Update grid colors
            for line in ax.xaxis.get_gridlines():
                line.set_color(mpl_colors['grid.color'])
            for line in ax.yaxis.get_gridlines():
                line.set_color(mpl_colors['grid.color'])
            # Update zero reference lines (tagged with gid='zero_ref')
            for line in ax.get_lines():
                if line.get_gid() == 'zero_ref':
                    line.set_color(zero_line_color)
            legend = ax.get_legend()
            if legend:
                legend.get_frame().set_facecolor(mpl_colors['legend.facecolor'])
                legend.get_frame().set_edgecolor(mpl_colors['legend.edgecolor'])
                for text in legend.get_texts():
                    text.set_color(mpl_colors['legend.labelcolor'])
