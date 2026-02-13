"""
Custom navigation toolbar for matplotlib
"""

import os
from PySide6.QtWidgets import QPushButton, QFileDialog, QMessageBox
from PySide6.QtCore import QSize
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT

from ui.theme import ThemeManager


class QuietNavigationToolbar(NavigationToolbar2QT):
    """
    Custom navigation toolbar with 300 DPI PNG/EPS save functionality
    and theme toggle button.
    """

    def __init__(self, canvas, parent=None):
        super().__init__(canvas, parent)
        self._add_theme_toggle()
        # Set initial icon
        self._update_theme_icon(ThemeManager.current_theme)
        # Listen for theme changes
        ThemeManager.on_theme_changed(self._update_theme_icon)
        # Clean up callback when toolbar is destroyed
        self.destroyed.connect(lambda: ThemeManager.off_theme_changed(self._update_theme_icon))

    def set_message(self, s):
        """Suppress coordinate display messages"""
        pass

    def save_figure(self, *args):
        """Save figure as displayed on screen"""
        initial_dir = os.path.expanduser("~")

        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Save figure",
            initial_dir,
            "PNG files (*.png);;All files (*.*)"
        )

        if not filepath:
            return

        base_path = os.path.splitext(filepath)[0]
        png_path = f"{base_path}.png"

        try:
            fig = self.canvas.figure
            fig.savefig(png_path, dpi=300, bbox_inches='tight',
                       facecolor=fig.get_facecolor(), edgecolor='none')

            QMessageBox.information(
                self, "Save Complete", f"Saved:\n\n{png_path}")

        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Failed to save figure:\n{str(e)}")

    def _add_theme_toggle(self):
        """Add theme toggle button (sun/moon) right after save icon in toolbar"""
        self.theme_btn = QPushButton()
        self.theme_btn.setFixedSize(QSize(32, 32))
        self.theme_btn.setToolTip("Toggle dark/light theme")
        self.theme_btn.setStyleSheet(
            "QPushButton { background: transparent; border: none; font-size: 18px; padding: 0; }"
            "QPushButton:hover { background: rgba(128,128,128,0.3); border-radius: 4px; }"
        )
        self.theme_btn.clicked.connect(ThemeManager.toggle_theme)

        # Insert right after the Save action (before coordinate label)
        actions = self.actions()
        # The last action is typically the locLabel (coordinates widget)
        # Insert separator + button before it so it stays next to Save
        if len(actions) >= 2:
            loc_action = actions[-1]
            self.insertSeparator(loc_action)
            self.insertWidget(loc_action, self.theme_btn)
        else:
            self.addSeparator()
            self.addWidget(self.theme_btn)

    def _update_theme_icon(self, theme_name):
        """Update theme toggle button icon and color based on current theme"""
        if theme_name == 'dark':
            self.theme_btn.setText("\u2600\uFE0F")  # Sun
            self.theme_btn.setToolTip("Switch to light theme")
            self.theme_btn.setStyleSheet(
                "QPushButton { background: transparent; border: none; font-size: 18px; "
                "padding: 0; color: #FFD700; }"
                "QPushButton:hover { background: rgba(128,128,128,0.3); border-radius: 4px; }"
            )
        else:
            self.theme_btn.setText("\U0001F319")  # Crescent moon
            self.theme_btn.setToolTip("Switch to dark theme")
            self.theme_btn.setStyleSheet(
                "QPushButton { background: transparent; border: none; font-size: 18px; "
                "padding: 0; color: #4A5568; }"
                "QPushButton:hover { background: rgba(128,128,128,0.3); border-radius: 4px; }"
            )
