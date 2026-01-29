#!/usr/bin/python3.8

"""
Plot management and styling
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Legend item limits
PROFILE_LEGEND_LIMIT = 30
TIMETRACE_LEGEND_LIMIT = 15


class ColorManager:
    """Manages color assignment using viridis colormap"""
    
    def get_colors_for_entries(self, entries):
        """Generate colors for entries using viridis colormap"""
        n_colors = len(entries)
        if n_colors == 0:
            return []
        elif n_colors == 1:
            return [plt.cm.viridis(0.5)]
        
        colors = [plt.cm.viridis(i / (n_colors - 1)) for i in range(n_colors)]
        return colors


def apply_legend_with_limit(ax, max_items, **kwargs):
    """Apply legend with item limit, showing '... +N more' if exceeded"""
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
        """Setup profile plot axes"""
        ax1 = figure.add_subplot(121)
        ax1.set_xlabel('x')
        ax1.set_ylabel(y1_label)
        
        ax2 = figure.add_subplot(122, sharex=ax1)
        ax2.set_xlabel('x')
        ax2.set_ylabel(y2_label)
        
        return ax1, ax2
    
    def setup_timetrace_plot(self, figure, y1_label, y2_label):
        """Setup time trace plot axes"""
        ax1 = figure.add_subplot(211)
        ax1.set_ylabel(y1_label)
        
        ax2 = figure.add_subplot(212, sharex=ax1)
        ax2.set_xlabel('Time [s]')
        ax2.set_ylabel(y2_label)
        
        # Add zero line for velocity-like parameters
        if 'v' in y2_label.lower():
            ax2.axhline(y=0, c='silver', ls='--')
        
        return ax1, ax2
    
    def apply_common_styling(self, ax1, ax2, plot_type='profile', skip_legend=False):
        """Apply common styling to axes with legend limit"""
        if plot_type == 'timetrace':
            max_items = TIMETRACE_LEGEND_LIMIT
        else:
            max_items = PROFILE_LEGEND_LIMIT
        
        for ax in [ax1, ax2]:
            if not skip_legend:
                apply_legend_with_limit(ax, max_items, frameon=False, fontsize=8)
            ax.grid(ls='--', lw=0.3, c='lightgray')