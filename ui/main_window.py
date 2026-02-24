"""
Main application window with sidebar navigation and on-demand tab initialization
"""

import subprocess
import os
import sys
from functools import partial

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTreeWidget, QTreeWidgetItem, QStackedWidget, QLabel,
    QPushButton, QStatusBar, QProgressBar, QMessageBox, QFrame, QStyle,
    QSizePolicy,
)
from PySide6.QtCore import Qt, Signal, QTimer
from PySide6.QtGui import QColor, QFont, QPixmap, QIcon

from config.app_config import AppConfig, VERSION, UPDATE_DATE, CONTACT_EMAIL, AUTHOR_NAME
from config.diagnostic_config import get_enabled_diagnostics
from config.user_settings import (
    load_settings, save_settings, show_update_popup
)
from ui.tab_factory import TabFactory
from ui.widgets.custom_toolbar import QuietNavigationToolbar
from ui.theme import ThemeManager

_ICON_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'icons')

_CATEGORY_ORDER = ["Profiles", "Time Traces", "Spectral Analysis", "Imaging", "Other"]

_TAB_TYPE_TO_CATEGORY = {
    'profile': "Profiles",
    'timetrace': "Time Traces",
    'spectrogram': "Spectral Analysis",
    'nmode': "Spectral Analysis",
    'tv': "Imaging",
    'tv_startup': "Imaging",
    'irvb': "Imaging",
}


def _categorize_tabs(tab_configs):
    """Categorize tab configs into ordered groups.

    Returns list of (category_name, [(tab_index, tab_name), ...]) in display order.
    """
    categories = {}
    for i, config in enumerate(tab_configs):
        cat = _TAB_TYPE_TO_CATEGORY.get(config['tab_type'], "Other")
        if cat not in categories:
            categories[cat] = []
        categories[cat].append((i, config['tab_name']))

    return [(name, categories[name]) for name in _CATEGORY_ORDER if name in categories]


class SidebarNav(QWidget):
    """Sidebar with categorized tree navigation"""
    tabSelected = Signal(int)

    def __init__(self, tab_configs):
        super().__init__()
        self.tab_configs = tab_configs
        self.setFixedWidth(140)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 8, 0, 0)

        # Logo + PRISM (horizontal, aligned to PRISM font height)
        title_row = QWidget()
        title_layout = QHBoxLayout(title_row)
        title_layout.setContentsMargins(0, 0, 0, 0)
        title_layout.setSpacing(6)
        title_layout.addStretch()

        self.logo_label = QLabel()
        self._update_logo()
        title_layout.addWidget(self.logo_label, alignment=Qt.AlignVCenter)

        title = QLabel(f'<span style="font-size:24px; font-weight:bold; color:#0d6efd;">PRISM</span>')
        title_layout.addWidget(title, alignment=Qt.AlignVCenter)
        title_layout.addStretch()
        layout.addWidget(title_row)

        # Version (below title, tight margin)
        ver_label = QLabel(f'v{VERSION}')
        ver_label.setAlignment(Qt.AlignCenter)
        ver_label.setContentsMargins(0, 0, 0, 0)
        ver_label.setStyleSheet("color: #888; font-size: 11px; margin-top: -4px;")
        layout.addWidget(ver_label)

        # Tree
        self.tree = QTreeWidget()
        self.tree.setObjectName("sidebar")
        self.tree.setHeaderHidden(True)
        self.tree.setRootIsDecorated(False)
        self.tree.setIndentation(0)
        self.tree.setAnimated(True)
        self.tree.setExpandsOnDoubleClick(False)
        layout.addWidget(self.tree)

        self._build_tree()
        self.tree.itemClicked.connect(self._on_click)

    def _build_tree(self):
        """Build categorized tree from tab configs"""
        for cat_name, items in _categorize_tabs(self.tab_configs):
            cat_item = QTreeWidgetItem(self.tree)
            cat_item.setText(0, f"  {cat_name}")
            cat_item.setFlags(Qt.ItemIsEnabled)
            font = cat_item.font(0)
            font.setBold(True)
            font.setPointSize(10)
            cat_item.setFont(0, font)
            cat_item.setForeground(0, QColor("#888888"))

            for tab_index, tab_name in items:
                child = QTreeWidgetItem(cat_item)
                child.setText(0, f"    {tab_name}")
                child.setData(0, Qt.UserRole, tab_index)

            cat_item.setExpanded(True)

    def _on_click(self, item, _col):
        tab_index = item.data(0, Qt.UserRole)
        if tab_index is not None:
            self.tabSelected.emit(tab_index)

    def _update_logo(self):
        """Update logo image based on current theme"""
        theme = ThemeManager.current_theme
        logo_file = 'prism-logo-dark.svg' if theme == 'dark' else 'prism-logo-light.svg'
        logo_path = os.path.join(_ICON_DIR, logo_file)
        pixmap = QPixmap(logo_path)
        if not pixmap.isNull():
            pixmap = pixmap.scaled(32, 32, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            self.logo_label.setPixmap(pixmap)

    def select_first(self):
        """Select the first tab item"""
        first_cat = self.tree.topLevelItem(0)
        if first_cat and first_cat.childCount() > 0:
            first_child = first_cat.child(0)
            self.tree.setCurrentItem(first_child)


class PRISMApp(QMainWindow):
    """Main application window for PRISM with sidebar navigation

    mode='':       Full PRISM with sidebar
    mode='select': Tab selector screen (button grid)
    """

    def __init__(self, mode=''):
        super().__init__()
        self.mode = mode
        self.config = AppConfig()

        # Set window icon
        icon_file = 'prism-logo-dark.svg' if ThemeManager.current_theme == 'dark' else 'prism-logo-light.svg'
        self.setWindowIcon(QIcon(os.path.join(_ICON_DIR, icon_file)))
        self._efit_loader = None
        self._plot_manager = None

        # Load user settings
        load_settings()

        # Tab management
        self.tab_configs = []
        self.tab_cache = {}
        self.tab_widgets = {}

        self._build_tab_configs()

        if self.mode == 'select':
            # Tab selector mode
            self._child_windows = []
            self.setWindowTitle(f'PRISM v{VERSION} - Select Viewer')
            self._create_selector_ui()
            self.setFixedWidth(360)
            self.adjustSize()
        else:
            # Full PRISM mode
            self.setWindowTitle(f'PRISM v{VERSION} - Plasma Research Integrated System for Multi-diagnostics')
            self.resize(1500, 700)
            self._create_widgets()

            # Register theme change callback
            ThemeManager.on_theme_changed(self._on_theme_changed)

            # Print startup message
            self._print_startup_message()

            # Defer first tab creation and update popup to after window is shown
            # This prevents blocking the UI for 3+ seconds during startup
            QTimer.singleShot(0, self._deferred_startup)

    def _deferred_startup(self):
        """Show update popup first (lightweight), then load first tab (heavy)"""
        show_update_popup(self)

        if self.tab_configs:
            self.sidebar.select_first()
            self._switch_tab(0)

    @property
    def efit_loader(self):
        if self._efit_loader is None:
            from data_loaders.efit_loader import EFITLoader
            self._efit_loader = EFITLoader(self.config)
        return self._efit_loader

    @property
    def plot_manager(self):
        if self._plot_manager is None:
            from plotting.plot_manager import PlotManager
            self._plot_manager = PlotManager(self.config)
        return self._plot_manager

    def _build_tab_configs(self):
        """Build full tab configuration list"""
        enabled_diagnostics = get_enabled_diagnostics()
        tab_types = ['profile', 'timetrace']

        for diag in enabled_diagnostics:
            for tab_type in tab_types:
                if not TabFactory.should_create_tab(diag, tab_type):
                    continue
                tab_name = TabFactory.get_tab_name(diag, tab_type)
                self.tab_configs.append({
                    'diagnostic': diag,
                    'tab_type': tab_type,
                    'tab_name': tab_name
                })

        for special_type in ['spectrogram', 'nmode', 'tv', 'tv_startup', 'irvb']:
            self.tab_configs.append({
                'diagnostic': None,
                'tab_type': special_type,
                'tab_name': TabFactory.get_tab_name(None, special_type)
            })

    # ------------------------------------------------------------------
    # Tab Selector Mode (prism select)
    # ------------------------------------------------------------------

    def _create_selector_ui(self):
        """Create tab selector screen with categorized button list"""
        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(30, 15, 30, 15)
        layout.setSpacing(6)

        # Logo + Title
        header = QWidget()
        header_layout = QHBoxLayout(header)
        header_layout.setContentsMargins(0, 0, 0, 0)
        header_layout.setSpacing(6)
        header_layout.addStretch()

        self._selector_logo = QLabel()
        theme = ThemeManager.current_theme
        logo_file = 'prism-logo-dark.svg' if theme == 'dark' else 'prism-logo-light.svg'
        logo_path = os.path.join(_ICON_DIR, logo_file)
        pixmap = QPixmap(logo_path)
        if not pixmap.isNull():
            pixmap = pixmap.scaled(40, 40, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            self._selector_logo.setPixmap(pixmap)
        header_layout.addWidget(self._selector_logo)

        title = QLabel(
            f'<span style="font-size:22px; font-weight:bold; color:#0d6efd;">PRISM</span>'
            f'  <span style="font-size:10px; color:#888;">v{VERSION}</span>'
        )
        header_layout.addWidget(title)
        header_layout.addStretch()
        layout.addWidget(header)

        # Categorized buttons (one row per category, equal-width buttons)
        for cat_name, items in _categorize_tabs(self.tab_configs):
            cat_label = QLabel(cat_name)
            cat_label.setStyleSheet(
                "color: #888; font-size: 10px; font-weight: bold; "
                "padding: 4px 0 1px 0;"
            )
            layout.addWidget(cat_label)

            btn_row = QHBoxLayout()
            btn_row.setSpacing(4)
            for tab_index, tab_name in items:
                btn = QPushButton(tab_name)
                btn.setFixedHeight(26)
                btn.setStyleSheet("font-size: 11px; padding: 2px 8px;")
                btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
                btn.setCursor(Qt.PointingHandCursor)
                btn.clicked.connect(partial(self._launch_single_tab, tab_index))
                btn_row.addWidget(btn, 1)
            layout.addLayout(btn_row)

        self.setCentralWidget(central)

        # Register theme callback so selector follows full PRISM theme
        def _on_selector_theme_changed(theme_name):
            logo_file = 'prism-logo-dark.svg' if theme_name == 'dark' else 'prism-logo-light.svg'
            logo_path = os.path.join(_ICON_DIR, logo_file)
            pixmap = QPixmap(logo_path)
            if not pixmap.isNull():
                pixmap = pixmap.scaled(40, 40, Qt.KeepAspectRatio, Qt.SmoothTransformation)
                self._selector_logo.setPixmap(pixmap)
            icon_file = 'prism-logo-dark.svg' if theme_name == 'dark' else 'prism-logo-light.svg'
            self.setWindowIcon(QIcon(os.path.join(_ICON_DIR, icon_file)))

        self._on_selector_theme_changed = _on_selector_theme_changed
        ThemeManager.on_theme_changed(_on_selector_theme_changed)

    def _launch_single_tab(self, tab_index):
        """Launch a single tab viewer in a new window (selector stays open)"""
        config = self.tab_configs[tab_index]
        tab_name = config['tab_name']
        tab_type = config['tab_type']

        # Build window title
        if tab_type == 'profile':
            title = f"{tab_name} Profile"
        elif tab_type == 'timetrace':
            title = f"{tab_name} Time Trace"
        else:
            title = tab_name

        # Create a new standalone window
        win = _SingleTabWindow(
            config, title, self.config, self.efit_loader, self.plot_manager
        )
        win.show()

        # Keep reference so it doesn't get garbage-collected
        self._child_windows.append(win)

        print(f"[PRISM] Launched: {title}")

    # ------------------------------------------------------------------
    # Full PRISM Mode
    # ------------------------------------------------------------------

    def _create_widgets(self):
        """Create main application layout"""
        central = QWidget()
        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Top area: sidebar + content
        top_area = QWidget()
        top_layout = QHBoxLayout(top_area)
        top_layout.setContentsMargins(0, 0, 0, 0)
        top_layout.setSpacing(0)

        # Sidebar
        self.sidebar = SidebarNav(self.tab_configs)
        self.sidebar.tabSelected.connect(self._switch_tab)

        # Content stack
        self.stack = QStackedWidget()

        top_layout.addWidget(self.sidebar)
        top_layout.addWidget(self.stack)

        main_layout.addWidget(top_area, stretch=1)

        # Bottom bar
        self._create_bottom_bar(main_layout)

        self.setCentralWidget(central)

        # Status bar
        self.progress = QProgressBar()
        self.progress.setFixedWidth(180)
        self.progress.hide()
        self.statusBar().addPermanentWidget(self.progress)
        self.statusBar().showMessage("Ready")

    def _create_bottom_bar(self, parent_layout):
        """Create bottom bar with toolbar and info"""
        self.bottom_frame = QWidget()
        self.bottom_frame.setFixedHeight(42)
        bottom_layout = QHBoxLayout(self.bottom_frame)
        bottom_layout.setContentsMargins(10, 2, 10, 2)

        # Toolbar container (left side)
        self.toolbar_frame = QWidget()
        self.toolbar_layout = QHBoxLayout(self.toolbar_frame)
        self.toolbar_layout.setContentsMargins(0, 0, 0, 0)
        bottom_layout.addWidget(self.toolbar_frame, stretch=1)

        self.toolbar = None

        # Right side: developer info + theme toggle + docs
        self.dev_label = QLabel(f"Developed by {AUTHOR_NAME} ({CONTACT_EMAIL})")
        self.dev_label.setStyleSheet("color: #888; font-size: 11px;")
        bottom_layout.addWidget(self.dev_label)

        docs_button = QPushButton("View Docs")
        docs_button.setFixedWidth(100)
        docs_button.clicked.connect(self._show_docs)
        bottom_layout.addWidget(docs_button)

        parent_layout.addWidget(self.bottom_frame)

    def _switch_tab(self, tab_index):
        """Switch to tab by index (create on demand)"""
        if tab_index not in self.tab_cache:
            self.statusBar().showMessage(f"Loading {self.tab_configs[tab_index]['tab_name']}...")
            self.progress.setRange(0, 0)
            self.progress.show()
            QApplication.processEvents()

            self._create_tab_content(tab_index)

            self.progress.hide()
            self.statusBar().showMessage("Ready")

        if tab_index in self.tab_widgets:
            self.stack.setCurrentWidget(self.tab_widgets[tab_index])
            self._update_toolbar(tab_index)

    def _create_tab_content(self, tab_index):
        """Create actual tab content when first accessed"""
        if tab_index in self.tab_cache:
            return self.tab_cache[tab_index]

        config = self.tab_configs[tab_index]

        tab = TabFactory.create_tab(
            None,
            self.config,
            config['diagnostic'],
            config['tab_type'],
            self.efit_loader,
            self.plot_manager
        )

        if tab is not None:
            self.tab_cache[tab_index] = tab
            self.tab_widgets[tab_index] = tab.frame
            self.stack.addWidget(tab.frame)

        return tab

    def _update_toolbar(self, tab_index):
        """Update toolbar for the current tab"""
        # Remove existing toolbar
        if self.toolbar is not None:
            self.toolbar.setParent(None)
            self.toolbar.deleteLater()
            self.toolbar = None

        current_tab = self.tab_cache.get(tab_index)
        if not current_tab:
            return

        if hasattr(current_tab, 'canvas') and current_tab.canvas:
            self.toolbar = QuietNavigationToolbar(current_tab.canvas, self.toolbar_frame)
            self.toolbar_layout.addWidget(self.toolbar)
            current_tab.toolbar = self.toolbar

    def _on_theme_changed(self, theme_name):
        """Handle theme change: refresh canvases and UI elements"""
        # Update sidebar logo
        if hasattr(self, 'sidebar'):
            self.sidebar._update_logo()

        # Refresh all cached matplotlib figures
        for tab in self.tab_cache.values():
            if hasattr(tab, 'canvas') and tab.canvas:
                ThemeManager.apply_theme_to_figure(tab.canvas.figure)
                tab.canvas.draw_idle()

    def _show_docs(self):
        """Open docs folder"""
        try:
            if not os.path.exists(self.config.DOCS_PATH):
                QMessageBox.information(
                    self, "Docs Not Found",
                    f"Docs folder not found.\n\nPlease contact {AUTHOR_NAME} ({CONTACT_EMAIL})."
                )
                return

            subprocess.run(["xdg-open", self.config.DOCS_PATH], check=True)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening docs: {str(e)}")

    def _print_startup_message(self):
        """Print startup message to console"""
        print("\n")
        print("+" + "=" * 62 + "+")
        print("|" + " " * 62 + "|")
        print("|" + f"PRISM v{VERSION} ({UPDATE_DATE})".center(62) + "|")
        print("|" + "Plasma Research Integrated System for Multi-diagnostics".center(62) + "|")
        print("|" + " " * 62 + "|")
        print("+" + "=" * 62 + "+")
        print("|" + " " * 62 + "|")
        print("|" + f"Developed by {AUTHOR_NAME}".center(62) + "|")
        print("|" + CONTACT_EMAIL.center(62) + "|")
        print("|" + " " * 62 + "|")
        print("+" + "=" * 62 + "+")
        print()
        print("  * prism       Launch full PRISM with all diagnostics")
        print("  * prism -s    Select and launch individual diagnostic viewers")
        print("  * prism -h    Show help")
        print()

    def closeEvent(self, event):
        """Handle window close event"""
        for tab_index, tab in self.tab_cache.items():
            if hasattr(tab, 'save_settings'):
                tab.save_settings()

        save_settings()
        event.accept()

    def run(self):
        """Start the application main loop (compatibility with old interface)"""
        self.show()


class _SingleTabWindow(QMainWindow):
    """Standalone window for a single tab launched from the selector"""

    def __init__(self, tab_config, title, app_config, efit_loader, plot_manager):
        super().__init__()
        self.setWindowTitle(f'PRISM v{VERSION} - {title}')
        self.resize(1400, 700)

        # Set window icon
        icon_file = 'prism-logo-dark.svg' if ThemeManager.current_theme == 'dark' else 'prism-logo-light.svg'
        self.setWindowIcon(QIcon(os.path.join(_ICON_DIR, icon_file)))

        self.tab = None
        self.toolbar = None

        # Central widget
        central = QWidget()
        main_layout = QVBoxLayout(central)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Content area (no sidebar)
        self.stack = QStackedWidget()
        main_layout.addWidget(self.stack, stretch=1)

        # Bottom bar
        self.bottom_frame = QWidget()
        self.bottom_frame.setFixedHeight(42)
        bottom_layout = QHBoxLayout(self.bottom_frame)
        bottom_layout.setContentsMargins(10, 2, 10, 2)

        self.toolbar_frame = QWidget()
        self.toolbar_layout = QHBoxLayout(self.toolbar_frame)
        self.toolbar_layout.setContentsMargins(0, 0, 0, 0)
        bottom_layout.addWidget(self.toolbar_frame, stretch=1)

        dev_label = QLabel(f"Developed by {AUTHOR_NAME} ({CONTACT_EMAIL})")
        dev_label.setStyleSheet("color: #888; font-size: 11px;")
        bottom_layout.addWidget(dev_label)

        main_layout.addWidget(self.bottom_frame)
        self.setCentralWidget(central)

        # Create the tab
        self.tab = TabFactory.create_tab(
            None, app_config, tab_config['diagnostic'],
            tab_config['tab_type'], efit_loader, plot_manager
        )

        if self.tab:
            self.stack.addWidget(self.tab.frame)

            # Setup toolbar
            if hasattr(self.tab, 'canvas') and self.tab.canvas:
                self.toolbar = QuietNavigationToolbar(
                    self.tab.canvas, self.toolbar_frame
                )
                self.toolbar_layout.addWidget(self.toolbar)
                self.tab.toolbar = self.toolbar

        # Register theme callback
        ThemeManager.on_theme_changed(self._on_theme_changed)

    def _on_theme_changed(self, theme_name):
        if self.tab and hasattr(self.tab, 'canvas') and self.tab.canvas:
            ThemeManager.apply_theme_to_figure(self.tab.canvas.figure)
            self.tab.canvas.draw_idle()

    def closeEvent(self, event):
        ThemeManager.off_theme_changed(self._on_theme_changed)
        if self.tab and hasattr(self.tab, 'save_settings'):
            self.tab.save_settings()
        save_settings()
        event.accept()
