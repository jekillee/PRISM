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
