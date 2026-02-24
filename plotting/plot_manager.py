"""
Plot management and styling
"""


# Legend item limits
PROFILE_LEGEND_LIMIT = 30
TIMETRACE_LEGEND_LIMIT = 15


class ColorManager:
    """Manages color assignment using matplotlib colormaps"""

    # Discrete colormaps: use fixed index (color N is always the same)
    DISCRETE_CMAPS = {'tab10', 'tab20', 'tab20b', 'tab20c',
                      'Set1', 'Set2', 'Set3', 'Paired', 'Accent',
                      'Dark2', 'Pastel1', 'Pastel2'}

    def get_colors_for_entries(self, entries, colormap='tab10'):
        """Generate colors for entries using specified colormap.

        For discrete colormaps (tab10, tab20, etc.): uses fixed indices
        so color assignment is stable regardless of entry count.
        For continuous colormaps (viridis, etc.): distributes colors evenly.
        """
        import matplotlib.pyplot as plt

        n_colors = len(entries)
        if n_colors == 0:
            return []

        cmap = plt.get_cmap(colormap)

        if colormap in self.DISCRETE_CMAPS:
            # Fixed index: entry 0 → color 0, entry 1 → color 1, ...
            # Use .colors list for discrete colormaps (avoids float [0,1] mapping)
            palette = cmap.colors  # list of RGB tuples
            n_total = len(palette)
            colors = [palette[i % n_total] for i in range(n_colors)]
        else:
            # Continuous: distribute evenly across the colormap
            if n_colors == 1:
                return [cmap(0.5)]
            colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]

        return colors


def apply_legend_with_limit(ax, max_items, **kwargs):
    """Apply legend with item limit, showing '... +N more' if exceeded"""
    from matplotlib.lines import Line2D

    handles, labels = ax.get_legend_handles_labels()

    if len(handles) == 0:
        return

    if len(handles) <= max_items:
        ax.legend(handles, labels, **kwargs)
    else:
        # Show only first max_items
        truncated_handles = handles[:max_items]
        truncated_labels = labels[:max_items]

        # Add "... +N more" indicator
        remaining = len(handles) - max_items
        dummy_handle = Line2D([], [], color='none', marker='', linestyle='')
        truncated_handles.append(dummy_handle)
        truncated_labels.append(f'... +{remaining} more')

        ax.legend(truncated_handles, truncated_labels, **kwargs)


class PlotManager:
    """Manages matplotlib plotting operations"""
    
    def __init__(self, config):
        self.config = config
        self.color_manager = ColorManager()
    
    def setup_profile_plot(self, figure, y1_label, y2_label):
        """Setup profile plot axes (side-by-side)"""
        figure.subplots_adjust(left=0.10, right=0.97, top=0.93, bottom=0.10, wspace=0.20)
        ax1 = figure.add_subplot(121)
        ax1.set_xlabel('x')
        ax1.set_ylabel(y1_label)

        ax2 = figure.add_subplot(122, sharex=ax1)
        ax2.set_xlabel('x')
        ax2.set_ylabel(y2_label)

        return ax1, ax2

    def setup_timetrace_plot(self, figure, y1_label, y2_label):
        """Setup time trace plot axes (stacked)"""
        figure.subplots_adjust(left=0.10, right=0.97, top=0.95, bottom=0.10, hspace=0.15)
        ax1 = figure.add_subplot(211)
        ax1.set_ylabel(y1_label)

        ax2 = figure.add_subplot(212, sharex=ax1)
        ax2.set_xlabel('Time [s]')
        ax2.set_ylabel(y2_label)

        # Add zero line for velocity-like parameters
        if 'v' in y2_label.lower():
            from ui.theme import ThemeManager
            zc = 'white' if ThemeManager.current_theme == 'dark' else 'gray'
            ax2.axhline(y=0, c=zc, ls='--', gid='zero_ref')

        return ax1, ax2
    
    def apply_common_styling(self, ax1, ax2, plot_type='profile', skip_legend=False,
                             legend_fontsize=8, label_fontsize=None, tick_fontsize=None):
        """Apply common styling to axes with legend limit"""
        if plot_type == 'timetrace':
            max_items = TIMETRACE_LEGEND_LIMIT
        else:
            max_items = PROFILE_LEGEND_LIMIT

        for ax in [ax1, ax2]:
            if not skip_legend:
                apply_legend_with_limit(ax, max_items, frameon=False, fontsize=legend_fontsize)
            if label_fontsize:
                ax.xaxis.label.set_fontsize(label_fontsize)
                ax.yaxis.label.set_fontsize(label_fontsize)
            if tick_fontsize:
                ax.tick_params(labelsize=tick_fontsize)
            import matplotlib as mpl
            ax.grid(ls='--', lw=0.3, c=mpl.rcParams.get('grid.color', '#444444'))