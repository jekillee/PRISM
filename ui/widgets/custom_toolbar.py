"""
Custom navigation toolbar for matplotlib
"""

import os
from PySide6.QtWidgets import QPushButton, QToolButton, QFileDialog, QMessageBox
from PySide6.QtCore import QSize
from PySide6.QtGui import QPixmap, QIcon, QPainter, QColor
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT

from ui.theme import ThemeManager


class QuietNavigationToolbar(NavigationToolbar2QT):
    """
    Custom navigation toolbar with 300 DPI PNG save, SVG save,
    and theme toggle button.
    """

    # Override toolitems to remove the default save button
    toolitems = [t for t in NavigationToolbar2QT.toolitems if t[0] != 'Save']

    def __init__(self, canvas, parent=None):
        super().__init__(canvas, parent)
        # Fix icon size consistency across themes
        self.setIconSize(QSize(24, 24))
        # Store original icons from QToolButtons for colorization
        self._original_icons = {}
        for child in self.findChildren(QToolButton):
            if not child.icon().isNull():
                self._original_icons[child] = child.icon()
        self._add_png_save_button()
        self._add_svg_save_button()
        self._add_theme_toggle()
        # Set initial icon and colorize toolbar icons
        self._update_theme_icon(ThemeManager.current_theme)
        # Listen for theme changes
        ThemeManager.on_theme_changed(self._update_theme_icon)
        # Clean up callback when toolbar is destroyed
        self.destroyed.connect(lambda: ThemeManager.off_theme_changed(self._update_theme_icon))

    def set_message(self, s):
        """Suppress coordinate display messages"""
        pass

    def _colorize_icon(self, icon, color):
        """Create a colorized version of an icon"""
        sizes = icon.availableSizes()
        if not sizes:
            sizes = [QSize(24, 24)]
        new_icon = QIcon()
        for size in sizes:
            pixmap = icon.pixmap(size)
            colored = QPixmap(pixmap.size())
            colored.fill(QColor(0, 0, 0, 0))
            painter = QPainter(colored)
            painter.drawPixmap(0, 0, pixmap)
            painter.setCompositionMode(QPainter.CompositionMode_SourceIn)
            painter.fillRect(colored.rect(), QColor(color))
            painter.end()
            new_icon.addPixmap(colored)
        return new_icon

    def _colorize_toolbar_icons(self, color):
        """Colorize all QToolButton icons (except theme toggle)"""
        for btn, original_icon in self._original_icons.items():
            colorized = self._colorize_icon(original_icon, color)
            btn.setIcon(colorized)

    def _add_png_save_button(self):
        """Add PNG save button as text label (replacing default save icon)"""
        self.png_btn = QPushButton("PNG")
        self.png_btn.setFixedSize(QSize(38, 32))
        self.png_btn.setToolTip("Save figure as PNG (300 DPI)")
        self.png_btn.setStyleSheet(
            "QPushButton { background: transparent; border: none; font-size: 11px; padding: 0; }"
            "QPushButton:hover { background: rgba(128,128,128,0.3); border-radius: 4px; }"
        )
        self.png_btn.clicked.connect(self._save_png)

        # Insert before coordinate label
        actions = self.actions()
        if len(actions) >= 2:
            loc_action = actions[-1]
            self.insertWidget(loc_action, self.png_btn)
        else:
            self.addWidget(self.png_btn)

    def _save_png(self):
        """Save figure as PNG with white background"""
        initial_dir = os.path.expanduser("~")

        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Save figure (PNG)",
            initial_dir,
            "PNG files (*.png);;All files (*.*)"
        )

        if not filepath:
            return

        base_path = os.path.splitext(filepath)[0]
        png_path = f"{base_path}.png"

        try:
            fig = self.canvas.figure
            fig.savefig(png_path, dpi=300,
                       facecolor='white', edgecolor='none')

            QMessageBox.information(
                self, "Save Complete", f"Saved:\n\n{png_path}")

        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Failed to save figure:\n{str(e)}")

    def _save_svg(self):
        """Save figure as SVG with white background"""
        initial_dir = os.path.expanduser("~")

        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Save figure (SVG)",
            initial_dir,
            "SVG files (*.svg);;All files (*.*)"
        )

        if not filepath:
            return

        base_path = os.path.splitext(filepath)[0]
        svg_path = f"{base_path}.svg"

        try:
            fig = self.canvas.figure
            fig.savefig(svg_path, format='svg',
                       facecolor='white', edgecolor='none')

            QMessageBox.information(
                self, "Save Complete", f"Saved:\n\n{svg_path}")

        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Failed to save figure:\n{str(e)}")

    def _add_svg_save_button(self):
        """Add SVG save button right after the PNG save icon"""
        self.svg_btn = QPushButton("SVG")
        self.svg_btn.setFixedSize(QSize(38, 32))
        self.svg_btn.setToolTip("Save figure as SVG")
        self.svg_btn.setStyleSheet(
            "QPushButton { background: transparent; border: none; font-size: 11px; padding: 0; }"
            "QPushButton:hover { background: rgba(128,128,128,0.3); border-radius: 4px; }"
        )
        self.svg_btn.clicked.connect(self._save_svg)

        # Insert right after the Save action (before coordinate label)
        actions = self.actions()
        if len(actions) >= 2:
            loc_action = actions[-1]
            self.insertWidget(loc_action, self.svg_btn)
        else:
            self.addWidget(self.svg_btn)

    def _add_theme_toggle(self):
        """Add theme toggle button (sun/moon) right after SVG button in toolbar"""
        self.theme_btn = QPushButton()
        self.theme_btn.setFixedSize(QSize(32, 32))
        self.theme_btn.setToolTip("Toggle dark/light theme")
        self.theme_btn.setStyleSheet(
            "QPushButton { background: transparent; border: none; font-size: 18px; padding: 0; }"
            "QPushButton:hover { background: rgba(128,128,128,0.3); border-radius: 4px; }"
        )
        self.theme_btn.clicked.connect(ThemeManager.toggle_theme)

        # Insert right after the SVG button (before coordinate label)
        actions = self.actions()
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
            dark_btn_style = (
                "QPushButton { background: transparent; border: none; font-size: 11px; "
                "padding: 0; color: #cccccc; }"
                "QPushButton:hover { background: rgba(128,128,128,0.3); border-radius: 4px; }"
            )
            self.png_btn.setStyleSheet(dark_btn_style)
            self.svg_btn.setStyleSheet(dark_btn_style)
            # Colorize toolbar icons to white
            self._colorize_toolbar_icons('#ffffff')
        else:
            self.theme_btn.setText("\U0001F319")  # Crescent moon
            self.theme_btn.setToolTip("Switch to dark theme")
            self.theme_btn.setStyleSheet(
                "QPushButton { background: transparent; border: none; font-size: 18px; "
                "padding: 0; color: #4A5568; }"
                "QPushButton:hover { background: rgba(128,128,128,0.3); border-radius: 4px; }"
            )
            light_btn_style = (
                "QPushButton { background: transparent; border: none; font-size: 11px; "
                "padding: 0; color: #555555; }"
                "QPushButton:hover { background: rgba(128,128,128,0.3); border-radius: 4px; }"
            )
            self.png_btn.setStyleSheet(light_btn_style)
            self.svg_btn.setStyleSheet(light_btn_style)
            # Colorize toolbar icons to black
            self._colorize_toolbar_icons('#000000')
