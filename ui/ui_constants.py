"""
Common UI constants for consistent styling across all tabs
"""

import os

from PySide6.QtWidgets import QStyle, QApplication, QFileDialog
from PySide6.QtCore import QSize, Qt
from PySide6.QtGui import QIcon

_ICONS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'icons')


def save_file_async(parent, caption, default_path, name_filter, on_selected):
    """Show a NON-MODAL, non-native "Save file" dialog; call on_selected(path)
    when the user picks a file.

    Unlike the static QFileDialog.getSaveFileName (modal + native), this keeps the
    main window movable and responsive while the dialog is open -- important over
    remote X11, where a modal native dialog locks onto the main window. The dialog
    is parented to `parent` (so it survives until closed) and shown asynchronously,
    so the caller returns immediately and the write runs in on_selected.
    """
    dlg = QFileDialog(parent, caption, default_path, name_filter)
    dlg.setAcceptMode(QFileDialog.AcceptSave)
    dlg.setOption(QFileDialog.DontUseNativeDialog, True)
    dlg.setWindowModality(Qt.NonModal)
    dlg.setAttribute(Qt.WA_DeleteOnClose, True)
    dlg.fileSelected.connect(on_selected)
    dlg.show()
    return dlg
_LISTBOX_ARROW_SIZE = QSize(10, 10)

# Control panel dimensions (pixels)
CONTROL_PANEL_WIDTH = 380


def get_icon(standard_pixmap):
    """Get a QStyle standard icon.

    Usage:
        from ui.ui_constants import get_icon
        button.setIcon(get_icon(QStyle.SP_ArrowUp))
    """
    return QApplication.style().standardIcon(standard_pixmap)


def apply_listbox_arrow_icons(add_btn, remove_btn):
    """Apply white right/left SVG arrow icons to Add/Remove listbox buttons."""
    add_btn.setIcon(QIcon(os.path.join(_ICONS_DIR, 'arrow-right-white.svg')))
    add_btn.setIconSize(_LISTBOX_ARROW_SIZE)
    remove_btn.setIcon(QIcon(os.path.join(_ICONS_DIR, 'arrow-left-white.svg')))
    remove_btn.setIconSize(_LISTBOX_ARROW_SIZE)


def apply_shot_arrow_icons(up_btn, down_btn):
    """Apply white up/down SVG arrow icons to shot-input up/down buttons."""
    up_btn.setIcon(QIcon(os.path.join(_ICONS_DIR, 'arrow-up-white.svg')))
    up_btn.setIconSize(_LISTBOX_ARROW_SIZE)
    down_btn.setIcon(QIcon(os.path.join(_ICONS_DIR, 'arrow-down-white.svg')))
    down_btn.setIconSize(_LISTBOX_ARROW_SIZE)


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
