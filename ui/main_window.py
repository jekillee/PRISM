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
    QGridLayout, QSizePolicy,
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QFont

from config.app_config import AppConfig, VERSION, UPDATE_DATE, CONTACT_EMAIL, AUTHOR_NAME
from config.diagnostic_config import get_enabled_diagnostics
from config.user_settings import (
    load_settings, save_settings, show_update_popup
)
from data_loaders.efit_loader import EFITLoader
from plotting.plot_manager import PlotManager
from ui.tab_factory import TabFactory
from ui.widgets.custom_toolbar import QuietNavigationToolbar
from ui.theme import ThemeManager


class SidebarNav(QWidget):
    """Sidebar with categorized tree navigation"""
    tabSelected = Signal(int)

    def __init__(self, tab_configs):
        super().__init__()
        self.tab_configs = tab_configs
        self.setFixedWidth(140)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 8, 0, 0)

        # Title with version
        title = QLabel(f'<p style="line-height:70%; margin:10;"><span style="font-size:35px; font-weight:bold;">PRISM</span></p><p style="line-height:100%; margin:0;"><span style="font-size:15px;">v{VERSION}</span></p>')
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("color: #0d6efd; font-weight: bold; padding: 0 0 0 0;")
        layout.addWidget(title)

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
        # Categorize tabs
        categories = {}
        for i, config in enumerate(self.tab_configs):
            tab_type = config['tab_type']
            if tab_type == 'profile':
                cat = "Profiles"
            elif tab_type == 'timetrace':
                cat = "Time Traces"
            elif tab_type in ('spectrogram', 'nmode'):
                cat = "Spectral Analysis"
            elif tab_type in ('tv', 'tv_startup', 'irvb'):
                cat = "Imaging"
            else:
                cat = "Other"

            if cat not in categories:
                categories[cat] = []
            categories[cat].append((i, config['tab_name']))

        # Maintain order
        cat_order = ["Profiles", "Time Traces", "Spectral Analysis", "Imaging", "Other"]
        for cat_name in cat_order:
            if cat_name not in categories:
                continue

            cat_item = QTreeWidgetItem(self.tree)
            cat_item.setText(0, f"  {cat_name}")
            cat_item.setFlags(Qt.ItemIsEnabled)
            font = cat_item.font(0)
            font.setBold(True)
            font.setPointSize(10)
            cat_item.setFont(0, font)
            cat_item.setForeground(0, QColor("#888888"))

            for tab_index, tab_name in categories[cat_name]:
                child = QTreeWidgetItem(cat_item)
                child.setText(0, f"    {tab_name}")
                child.setData(0, Qt.UserRole, tab_index)

            cat_item.setExpanded(True)

    def _on_click(self, item, _col):
        tab_index = item.data(0, Qt.UserRole)
        if tab_index is not None:
            self.tabSelected.emit(tab_index)

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
        self.efit_loader = EFITLoader(self.config)
        self.plot_manager = PlotManager(self.config)

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
            self.resize(500, 400)
            self._create_selector_ui()
        else:
            # Full PRISM mode
            self.setWindowTitle(f'PRISM v{VERSION} - Plasma Research Integrated System for Multi-diagnostics')
            self.resize(1500, 700)
            self._create_widgets()

            # Register theme change callback
            ThemeManager.on_theme_changed(self._on_theme_changed)

            # Print startup message
            self._print_startup_message()

            # Create first tab immediately
            if self.tab_configs:
                self.sidebar.select_first()
                self._switch_tab(0)

            # Show update popup if new version
            show_update_popup(self)

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
        """Create tab selector screen with button grid"""
        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(30, 20, 30, 20)
        layout.setSpacing(15)

        # Title
        title = QLabel(
            f'<span style="font-size:28px; font-weight:bold; color:#0d6efd;">PRISM</span>'
            f'  <span style="font-size:12px; color:#888;">v{VERSION}</span>'
        )
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        subtitle = QLabel("Select a viewer to launch")
        subtitle.setAlignment(Qt.AlignCenter)
        subtitle.setStyleSheet("color: #888; font-size: 13px; margin-bottom: 10px;")
        layout.addWidget(subtitle)

        # Button grid
        grid = QGridLayout()
        grid.setSpacing(8)

        # Define buttons: (label, tab_indices)
        buttons = []
        for i, config in enumerate(self.tab_configs):
            tab_type = config['tab_type']
            tab_name = config['tab_name']

            # Build display label
            if tab_type == 'profile':
                label = f"{tab_name} Profile"
            elif tab_type == 'timetrace':
                label = f"{tab_name} Time Trace"
            else:
                label = tab_name

            buttons.append((label, i))

        # Layout: 2 columns
        for idx, (label, tab_index) in enumerate(buttons):
            btn = QPushButton(label)
            btn.setMinimumHeight(38)
            btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            btn.setCursor(Qt.PointingHandCursor)
            btn.clicked.connect(partial(self._launch_single_tab, tab_index))
            grid.addWidget(btn, idx // 2, idx % 2)

        layout.addLayout(grid)
        layout.addStretch()

        # Developer info
        dev_label = QLabel(f"Developed by {AUTHOR_NAME} ({CONTACT_EMAIL})")
        dev_label.setAlignment(Qt.AlignCenter)
        dev_label.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(dev_label)

        self.setCentralWidget(central)

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

        # Right side: developer info + theme toggle + manual
        self.dev_label = QLabel(f"Developed by {AUTHOR_NAME} ({CONTACT_EMAIL})")
        self.dev_label.setStyleSheet("color: #888; font-size: 11px;")
        bottom_layout.addWidget(self.dev_label)

        manual_button = QPushButton("View Manual")
        manual_button.setFixedWidth(100)
        manual_button.clicked.connect(self._show_manual)
        bottom_layout.addWidget(manual_button)

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
        # Refresh all cached matplotlib figures
        for tab in self.tab_cache.values():
            if hasattr(tab, 'canvas') and tab.canvas:
                ThemeManager.apply_theme_to_figure(tab.canvas.figure)
                tab.canvas.draw_idle()

    def _show_manual(self):
        """Open user manual PDF"""
        try:
            if not os.path.exists(self.config.MANUAL_PATH):
                QMessageBox.information(
                    self, "Manual Not Found",
                    f"User manual not found.\n\nPlease contact {AUTHOR_NAME} ({CONTACT_EMAIL}) for the manual."
                )
                return

            subprocess.run(["xdg-open", self.config.MANUAL_PATH], check=True)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening manual: {str(e)}")

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
