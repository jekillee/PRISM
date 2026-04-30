"""
Common UI constants for consistent styling across all tabs
"""

from PySide6.QtWidgets import QStyle, QApplication

# Control panel dimensions (pixels)
CONTROL_PANEL_WIDTH = 380


def get_icon(standard_pixmap):
    """Get a QStyle standard icon.

    Usage:
        from ui.ui_constants import get_icon
        button.setIcon(get_icon(QStyle.SP_ArrowUp))
    """
    return QApplication.style().standardIcon(standard_pixmap)


def apply_dark_figure_style(figure):
    """Apply current theme styling to a matplotlib Figure.

    Note: Most styling is handled globally via rcParams.
    This refreshes figure-level colors for dynamically created figures.
    """
    from ui.theme import ThemeManager
    ThemeManager.apply_theme_to_figure(figure)


_STATUS_COLOR_MAP = {
    'red': '#d32f2f', 'green': '#2e7d32', 'blue': '#1976d2',
    'gray': '#888', 'orange': '#e65100',
}


def show_status(widget, prefix, text, color='blue', timeout=0):
    """Push a colored status message to the main window's status bar.

    Args:
        widget: any widget inside the main window (used to find the top-level)
        prefix: short tag shown in brackets, e.g. 'TV', 'EFIT', 'Spectrogram'
        text: status text
        color: 'red' | 'green' | 'blue' | 'gray' | 'orange', or any CSS color
        timeout: ms before auto-clear; 0 = persistent until replaced
    """
    from PySide6.QtWidgets import QApplication
    try:
        w = widget.window()
        if hasattr(w, 'statusBar'):
            sb = w.statusBar()
            css_color = _STATUS_COLOR_MAP.get(color, color)
            sb.setStyleSheet(
                f"QStatusBar {{ color: {css_color}; font-weight: bold; }}")
            sb.showMessage(f"[{prefix}] {text}", timeout)
    except Exception:
        pass
    print(f"[{prefix}] {text}")
    QApplication.processEvents()


def make_or_divider(text="or"):
    """Create a horizontal divider widget with centered text.

    Visual: ───────  or  ───────

    Used to separate alternative input methods (e.g. Fetch from MDS+ vs
    open local file) so it's clear that the user picks one path.

    Returns:
        QWidget ready to be inserted into a layout (typically spanning the
        full width of a QGridLayout via addWidget(w, row, 0, 1, ncols)).
    """
    from PySide6.QtWidgets import QWidget, QHBoxLayout, QFrame, QLabel

    widget = QWidget()
    layout = QHBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(6)

    line_left = QFrame()
    line_left.setFrameShape(QFrame.HLine)
    line_left.setFrameShadow(QFrame.Sunken)
    line_left.setStyleSheet("color: #888;")
    layout.addWidget(line_left, stretch=1)

    label = QLabel(text)
    label.setStyleSheet("color: #888; font-style: italic;")
    layout.addWidget(label)

    line_right = QFrame()
    line_right.setFrameShape(QFrame.HLine)
    line_right.setFrameShadow(QFrame.Sunken)
    line_right.setStyleSheet("color: #888;")
    layout.addWidget(line_right, stretch=1)

    return widget
