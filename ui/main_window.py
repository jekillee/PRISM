"""
Main application window with sidebar navigation and on-demand tab initialization
"""

import subprocess
import os
import sys
from functools import partial

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTreeWidget, QTreeWidgetItem, QStackedWidget, QLabel, QLineEdit,
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

_CATEGORY_ORDER = ["Profiles", "Time Traces", "Spectral", "Imaging", "TRANSP", "Other"]

# Tab types highlighted with accent color in sidebar
_NEW_TABS = set()  # no accent-colored tabs

_TAB_TYPE_TO_CATEGORY = {
    'profile': "Profiles",
    'timetrace': "Time Traces",
    'spectrogram': "Spectral",
    'nmode': "Spectral",
    'tv': "Imaging",
    'tv_startup': "Imaging",
    'irvb': "Imaging",
    'neutron': "Time Traces",
    'bi_profile': "Profiles",
    'bi_timetrace': "Time Traces",
    'transp_profile': "TRANSP",
    'transp_timetrace': "TRANSP",
}


def _categorize_tabs(tab_configs):
    """Categorize tab configs into ordered groups.

    Returns list of (category_name, [(tab_index, tab_name, tab_type), ...]) in display order.
    """
    categories = {}
    for i, config in enumerate(tab_configs):
        cat = _TAB_TYPE_TO_CATEGORY.get(config['tab_type'], "Other")
        if cat not in categories:
            categories[cat] = []
        categories[cat].append((i, config['tab_name'], config['tab_type']))

    return [(name, categories[name]) for name in _CATEGORY_ORDER if name in categories]


class SidebarNav(QWidget):
    """Sidebar with categorized tree navigation"""
    tabSelected = Signal(int)

    def __init__(self, tab_configs, mode=''):
        super().__init__()
        self.tab_configs = tab_configs
        self.mode = mode
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

        # Version + date (below title, tight margin)
        ver_label = QLabel(f'v{VERSION} ({UPDATE_DATE})')
        ver_label.setAlignment(Qt.AlignCenter)
        ver_label.setContentsMargins(0, 0, 0, 0)
        ver_label.setStyleSheet("color: #888; font-size: 11px; margin-top: -8px;")
        layout.addWidget(ver_label)

        # Global shot input (single row: apply button + entry)
        shot_widget = QWidget()
        shot_row = QHBoxLayout(shot_widget)
        shot_row.setContentsMargins(8, 6, 8, 2)
        shot_row.setSpacing(3)

        self.global_shot_entry = QLineEdit()
        self.global_shot_entry.setAlignment(Qt.AlignCenter)
        self.global_shot_entry.setFixedHeight(24)
        self.global_shot_entry.setPlaceholderText("Enter Shot")
        self.global_shot_entry.setFocusPolicy(Qt.ClickFocus)
        shot_row.addWidget(self.global_shot_entry, stretch=1)

        self.apply_all_btn = QPushButton()
        self.apply_all_btn.setIcon(self.style().standardIcon(QStyle.SP_DialogOkButton))
        self.apply_all_btn.setFixedSize(24, 24)
        self.apply_all_btn.setToolTip("Apply shot to all tabs")
        shot_row.addWidget(self.apply_all_btn)

        layout.addWidget(shot_widget)

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
        self.tree.itemExpanded.connect(self._on_category_toggled)
        self.tree.itemCollapsed.connect(self._on_category_toggled)

    def _build_tree(self):
        """Build categorized tree from tab configs"""
        if self.mode == 'transp':
            self._build_bi_tree()
            return

        # Load saved expand states
        from config.user_settings import get_tab_settings
        saved = get_tab_settings('sidebar')
        collapsed_cats = saved.get('collapsed_categories', [])

        self._cat_items = {}  # {cat_name: QTreeWidgetItem}

        for cat_name, items in _categorize_tabs(self.tab_configs):
            expanded = cat_name not in collapsed_cats
            arrow = '\u25bc' if expanded else '\u25b6'  # ▼ or ▶

            cat_item = QTreeWidgetItem(self.tree)
            cat_item.setText(0, f" {arrow} {cat_name}")
            cat_item.setFlags(Qt.ItemIsEnabled)
            font = cat_item.font(0)
            font.setBold(True)
            font.setPointSize(10)
            cat_item.setFont(0, font)
            cat_item.setForeground(0, QColor("#888888"))
            cat_item.setData(0, Qt.UserRole + 1, cat_name)

            for tab_index, tab_name, tab_type in items:
                child = QTreeWidgetItem(cat_item)
                child.setText(0, f"    {tab_name}")
                child.setData(0, Qt.UserRole, tab_index)

            cat_item.setExpanded(expanded)
            self._cat_items[cat_name] = cat_item

    def _on_category_toggled(self, item):
        """Save category expand/collapse state and update arrow"""
        cat_key = item.data(0, Qt.UserRole + 1)
        if cat_key is None:
            return
        # Display label: use part after '/' for sub-categories, full key for roots
        display = cat_key.split('/')[-1] if '/' in cat_key else cat_key
        # Indent: sub-categories get extra space
        indent = "  " if '/' in cat_key else " "
        arrow = '\u25bc' if item.isExpanded() else '\u25b6'
        item.setText(0, f"{indent}{arrow} {display}")
        # Save state
        from config.user_settings import get_tab_settings, set_tab_settings
        saved = get_tab_settings('sidebar')
        collapsed = set(saved.get('collapsed_categories', []))
        if item.isExpanded():
            collapsed.discard(cat_key)
        else:
            collapsed.add(cat_key)
        saved['collapsed_categories'] = list(collapsed)
        set_tab_settings('sidebar', saved)

    def _build_bi_tree(self):
        """Build TRANSP viewer sidebar with collapsible root nodes."""
        from config.user_settings import get_tab_settings
        saved = get_tab_settings('sidebar')
        collapsed_cats = saved.get('collapsed_categories', [])

        self._cat_items = {}

        _BI_CATS = {"Profiles", "Time Traces"}
        categorized = _categorize_tabs(self.tab_configs)
        bi_cats = [(n, items) for n, items in categorized if n in _BI_CATS]
        transp_cats = [(n, items) for n, items in categorized if n not in _BI_CATS]

        def _add_root(label, cats, flat=False):
            expanded = label not in collapsed_cats
            arrow = '\u25bc' if expanded else '\u25b6'

            root = QTreeWidgetItem(self.tree)
            root.setText(0, f" {arrow} {label}")
            root.setFlags(Qt.ItemIsEnabled)
            root.setData(0, Qt.UserRole + 1, label)
            font = root.font(0)
            font.setBold(True)
            font.setPointSize(10)
            root.setFont(0, font)
            root.setForeground(0, QColor("#888888"))
            self._cat_items[label] = root

            if flat:
                # Tabs directly under root (no sub-categories)
                for _, items in cats:
                    for tab_index, tab_name, tab_type in items:
                        child = QTreeWidgetItem(root)
                        child.setText(0, f"    {tab_name}")
                        child.setData(0, Qt.UserRole, tab_index)
            else:
                for cat_name, items in cats:
                    cat_key = f"{label}/{cat_name}"
                    cat_expanded = cat_key not in collapsed_cats
                    cat_arrow = '\u25bc' if cat_expanded else '\u25b6'

                    cat_item = QTreeWidgetItem(root)
                    cat_item.setText(0, f"  {cat_arrow} {cat_name}")
                    cat_item.setFlags(Qt.ItemIsEnabled)
                    cat_item.setData(0, Qt.UserRole + 1, cat_key)
                    cat_font = cat_item.font(0)
                    cat_font.setBold(True)
                    cat_item.setFont(0, cat_font)
                    cat_item.setForeground(0, QColor("#888888"))
                    self._cat_items[cat_key] = cat_item

                    for tab_index, tab_name, tab_type in items:
                        child = QTreeWidgetItem(cat_item)
                        child.setText(0, f"      {tab_name}")
                        child.setData(0, Qt.UserRole, tab_index)

                    cat_item.setExpanded(cat_expanded)

            root.setExpanded(expanded)

        _add_root("BiProfile", bi_cats)
        if transp_cats:
            _add_root("TRANSP", transp_cats, flat=True)

    def _on_click(self, item, _col):
        tab_index = item.data(0, Qt.UserRole)
        if tab_index is not None:
            self.tabSelected.emit(tab_index)
        else:
            # Category item clicked — toggle expand/collapse
            cat_name = item.data(0, Qt.UserRole + 1)
            if cat_name is not None:
                item.setExpanded(not item.isExpanded())

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
        if self.mode == 'transp':
            # 3-level: BiProfile > category > tab
            bi = self.tree.topLevelItem(0)
            if bi and bi.childCount() > 0:
                cat = bi.child(0)
                if cat and cat.childCount() > 0:
                    self.tree.setCurrentItem(cat.child(0))
            return
        first_cat = self.tree.topLevelItem(0)
        if first_cat and first_cat.childCount() > 0:
            first_child = first_cat.child(0)
            self.tree.setCurrentItem(first_child)


class PRISMApp(QMainWindow):
    """Main application window for PRISM with sidebar navigation

    mode='':       Full PRISM with sidebar
    mode='select': Tab selector screen (button grid)
    mode='transp':     TRANSP viewer (transport analysis profiles)
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

        if self.mode == 'transp':
            self._build_bi_tab_configs()
        else:
            self._build_tab_configs()

        # Print startup message for all modes
        self._print_startup_message()

        if self.mode == 'select':
            # Tab selector mode
            self._child_windows = []
            self.setWindowTitle(f'PRISM v{VERSION} - Select Viewer')
            self._create_selector_ui()
            self.setFixedWidth(360)
            self.adjustSize()
            QTimer.singleShot(0, lambda: show_update_popup(self))
        elif self.mode == 'transp':
            # TRANSP viewer mode
            self._bi_shot_data = {}
            self.setWindowTitle(f'PRISM v{VERSION} - TRANSP Viewer')
            self.resize(1500, 700)
            self._create_widgets()
            ThemeManager.on_theme_changed(self._on_theme_changed)
            QTimer.singleShot(0, self._deferred_startup_bi)
        else:
            # Full PRISM mode
            self.setWindowTitle(f'PRISM v{VERSION} - Plasma Research Integrated System for Multi-diagnostics')
            self.resize(1500, 700)
            self._create_widgets()

            # Register theme change callback
            ThemeManager.on_theme_changed(self._on_theme_changed)

            # Defer first tab creation and update popup to after window is shown
            QTimer.singleShot(0, self._deferred_startup)

    def _deferred_startup(self):
        """Show update popup first (lightweight), then load first tab (heavy)"""
        show_update_popup(self)

        if self.tab_configs:
            self.sidebar.select_first()
            self._switch_tab(0)

    def _deferred_startup_bi(self):
        """TRANSP viewer mode startup"""
        show_update_popup(self)
        if self.tab_configs:
            self.sidebar.select_first()
            self._switch_tab(0)

    # ------------------------------------------------------------------
    # BiProfile mode
    # ------------------------------------------------------------------

    def _build_bi_tab_configs(self):
        """Build tab configs for TRANSP viewer mode"""
        self.tab_configs = [
            {'diagnostic': None, 'tab_type': 'bi_profile',
             'tab_name': 'Ti, vT', 'bi_params': ('Ti', 'vT')},
            {'diagnostic': None, 'tab_type': 'bi_profile',
             'tab_name': 'ne, Te', 'bi_params': ('ne', 'Te')},
            {'diagnostic': None, 'tab_type': 'bi_timetrace',
             'tab_name': 'Ti, vT', 'bi_params': ('Ti', 'vT')},
            {'diagnostic': None, 'tab_type': 'bi_timetrace',
             'tab_name': 'ne, Te', 'bi_params': ('ne', 'Te')},
            {'diagnostic': None, 'tab_type': 'transp_profile',
             'tab_name': 'Profile'},
            {'diagnostic': None, 'tab_type': 'transp_timetrace',
             'tab_name': 'Time Trace'},
        ]

    def fetch_biprofile_shot(self, shot):
        """Load all biprofile data for a shot (parallel). Cached.

        Two-phase: BIPROFILE + DIAG_PARAMS + CES + Thomson in parallel,
        then EFIT using tree from BIPROFILE metadata.
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed
        from data_loaders.biprofile_loader import (
            load_biprofile, load_diag_params,
            load_efit_psin, load_ces_raw, load_thomson_raw,
        )

        self.statusBar().showMessage(f"Loading #{shot}...")
        self.progress.setRange(0, 0)
        self.progress.show()
        QApplication.processEvents()

        bi_data = {}
        diag = None
        ces = None
        thomson = None
        with ThreadPoolExecutor(max_workers=8) as executor:
            futures = {}
            for param in ['Ti', 'vT', 'Te', 'ne']:
                futures[executor.submit(load_biprofile, shot, param)] = param
            futures[executor.submit(load_diag_params, shot)] = 'DIAG_PARAMS'
            futures[executor.submit(load_ces_raw, shot)] = 'CES'
            futures[executor.submit(load_thomson_raw, shot)] = 'Thomson'

            for future in as_completed(futures):
                key = futures[future]
                try:
                    result = future.result()
                    if key == 'DIAG_PARAMS':
                        diag = result
                    elif key == 'CES':
                        ces = result
                    elif key == 'Thomson':
                        thomson = result
                    else:
                        bi_data[key] = result
                except Exception as e:
                    print(f"[Loader] BiProfile {key} failed: {e}")

        if not bi_data:
            self.progress.hide()
            self.statusBar().showMessage(f"#{shot}: No data loaded")
            return None

        # Phase 2: EFIT
        efit = None
        efit_tree = 'EFIT01'
        for param_data in bi_data.values():
            if param_data.get('efit_used'):
                efit_tree = param_data['efit_used']
                break
        try:
            efit = load_efit_psin(shot, efit_tree)
        except Exception as e:
            print(f"[EFIT] {efit_tree} failed: {e}")

        self.progress.hide()
        self._bi_shot_data[shot] = {
            'bi': bi_data, 'diag': diag,
            'efit': efit, 'ces': ces, 'thomson': thomson,
        }
        loaded = list(bi_data.keys())
        self.statusBar().showMessage(f"#{shot}: {', '.join(loaded)} loaded")
        return self._bi_shot_data[shot]

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

        # Tabs that require nkstar-local resources (TV images, IRVB HTTP server)
        import socket as _socket
        _host = _socket.gethostname()
        _nkstar_only = {'tv', 'tv_startup', 'irvb'}

        for special_type in ['spectrogram', 'nmode', 'tv', 'tv_startup', 'irvb', 'neutron']:
            if special_type in _nkstar_only and not _host.startswith('nkstar'):
                continue
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
            for tab_index, tab_name, tab_type in items:
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
        self.sidebar = SidebarNav(self.tab_configs, mode=self.mode)
        self.sidebar.tabSelected.connect(self._switch_tab)
        self.sidebar.apply_all_btn.clicked.connect(self._apply_shot_to_all_tabs)

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

    def _apply_shot_to_all_tabs(self):
        """Apply global shot number to all tabs (creating uncached ones first)"""
        shot_text = self.sidebar.global_shot_entry.text().strip()
        if not shot_text:
            return

        self.statusBar().showMessage(f"Applying shot #{shot_text} to all tabs...")
        self.progress.setRange(0, 0)
        self.progress.show()
        QApplication.processEvents()

        count = 0
        for i in range(len(self.tab_configs)):
            if i not in self.tab_cache:
                self._create_tab_content(i)
            tab = self.tab_cache.get(i)
            if tab and hasattr(tab, 'shot_entry'):
                tab.shot_entry.setText(shot_text)
                count += 1

        self.progress.hide()
        self.statusBar().showMessage(f"Shot #{shot_text} applied to {count} tab(s)")

    def _create_tab_content(self, tab_index):
        """Create actual tab content when first accessed"""
        if tab_index in self.tab_cache:
            return self.tab_cache[tab_index]

        config = self.tab_configs[tab_index]

        # BiProfile / TRANSP tabs (custom, not via TabFactory)
        if config['tab_type'] in ('bi_profile', 'bi_timetrace'):
            tab = self._create_bi_tab(config)
        elif config['tab_type'] in ('transp_profile', 'transp_timetrace'):
            tab = self._create_transp_tab(config)
        else:
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

    def _create_bi_tab(self, config):
        """Create a BiProfile tab instance"""
        if config['tab_type'] == 'bi_profile':
            from ui.biprofile_profile_tab import BiProfileTab
            return BiProfileTab(self, config['bi_params'])
        else:
            from ui.biprofile_timetrace_tab import BiTimeTraceTab
            return BiTimeTraceTab(self, config['bi_params'])

    def _create_transp_tab(self, config):
        """Create a TRANSP CDF tab instance"""
        if config['tab_type'] == 'transp_profile':
            from ui.transp_profile_tab import TranspProfileTab
            return TranspProfileTab(self)
        else:
            from ui.transp_timetrace_tab import TranspTimeTraceTab
            return TranspTimeTraceTab(self)

    def _update_toolbar(self, tab_index):
        """Update toolbar for the current tab"""
        # Disconnect and remove existing toolbar
        if self.toolbar is not None:
            # Disconnect matplotlib callbacks before destroying
            if hasattr(self.toolbar, '_id_press'):
                self.toolbar.canvas.mpl_disconnect(self.toolbar._id_press)
            if hasattr(self.toolbar, '_id_release'):
                self.toolbar.canvas.mpl_disconnect(self.toolbar._id_release)
            if hasattr(self.toolbar, '_id_drag'):
                self.toolbar.canvas.mpl_disconnect(self.toolbar._id_drag)
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
        print("  * prism -t    Launch TRANSP viewer")
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
