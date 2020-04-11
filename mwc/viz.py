import bokeh.io
import bokeh.plotting
import bokeh.layouts
import bokeh.palettes
import bokeh.themes.theme
import itertools
import holoviews as hv
import skimage.io
import skimage.measure
import seaborn as sns
import numpy as np
import pandas as pd
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.signal import gaussian, convolve
from matplotlib.path import Path
from matplotlib.patches import BoxStyle
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import make_axes_locatable


def despine(ax, offset=2):
    if (type(ax) != np.ndarray) & (type(ax) != list):
        ax = [ax]
    for a in ax:
         a.spines['bottom'].set_position(('outward', 2))
         a.spines['left'].set_position(('outward', 2))

def color_palette():
    """ 
    Returns my preferred color scheme 
    """
    colors = {
        "dark_purple": "#5F2E88",
        "purple": "#7E59A2",
        "light_purple": "#A17DB8",
        "pale_purple": "#dfd6e5",
        "dark_orange": "#F38227",
        "orange": "#E39943",
        "light_orange": "#EEBA7F",
        "pale_orange": "#f2d4b6",
        "dark_blue": "#3F60AC",
        "blue": "#7292C7",
        "light_blue": "#A5B3CC",
        "pale_blue": "#dae4f1",
        "dark_red": "#9C372F",
        "red": "#C76A6A",
        "light_red": "#E39C9D",
        "pale_red": "#edcccc",
        "dark_green": "#395A34",
        "green": "#688A2F",
        "light_green": "#B3CD86",
        "pale_green": "#d8e2c3",
        "dark_brown": "#764f2a",
        "brown": "#c2996c",
        "light_brown": "#e1bb96",
        "pale_brown": "#efccaf",
        "black": "#444147",
        "grey": "#EFEFEF",
        "gray": "#EFEFEF",
        "light_grey": "#6D6F72",
        "light_gray": "#6D6F72",
    }
    palette = [
        v
        for k, v in colors.items()
        if (v not in ["grey", "gray", "dark_purple", "light_grey"])
        and ("pale" not in [v])
    ]
    return colors, palette


def plotting_style():
    """
    Sets the plotting style to my preference
    """
    rc = {
        "axes.facecolor": "#EFEFEF",
        "font.family": "sans-serif",
        "font.family": "Myriad Pro",
        "font.style": "normal",
        "pdf.fonttype": 42,
        "axes.edgecolor": "#444147",
        "axes.labelcolor": "#444147",
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        "axes.axisbelow": True,
        "axes.linewidth": 0.25,
        "axes.titlesize": 8,
        "text.color": "#444147",
        "axes.grid": False,
        "lines.linewidth": 0.75,
        "lines.dash_capstyle": "round",
        "patch.linewidth": 0.25,
        "axes.labelsize": 8,
        "xtick.color": "#444147",
        "ytick.color": "#444147",
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.25,
        "ytick.major.width": 0.25,
        "xtick.major.pad": 6,
        "ytick.major.pad": 6,
        "xtick.minor.size": 0,
        "ytick.minor.size": 0,
        "legend.fontsize": 6,
        "legend.frameon": True,
        "legend.edgecolor": "#444147",
        "axes.xmargin": 0.03,
        "axes.ymargin": 0.03,
        "figure.facecolor": "white",
        "figure.dpi": 300,
        "errorbar.capsize": 1,
        "savefig.bbox": "tight",
    }

    plt.rc("text.latex", preamble=r"\usepackage{mathpazo}")
    matplotlib.style.use(rc)
    colors, palette = color_palette()
    sns.set_palette(palette)
    return colors, palette


def titlebox(ax, text, color='white',  bgcolor=None, size=8, boxsize="12%", pad=0.05,
            loc=10, **kwargs):
    """Sets a colored box about the title with the width of the plot"""
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size=boxsize, pad=pad)
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.spines['top'].set_visible(True)
    cax.spines['right'].set_visible(True)
    plt.setp(cax.spines.values(), color=color)
    if bgcolor != None:
        cax.set_facecolor(bgcolor) 
    else:
        cax.set_facecolor('white')
    at = AnchoredText(text, loc=loc, frameon=False, prop=dict(size=size, color=color))
    cax.add_artist(at)

def ylabelbox(ax, text, color,  bgcolor=None, size=6, boxsize="15%", pad=0.02, **kwargs):
    """Sets a colored box about the title with the width of the plot"""
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("left", size=boxsize, pad=pad)
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.spines['top'].set_visible(True)
    cax.spines['right'].set_visible(True)
    plt.setp(cax.spines.values(), color=color)
    if bgcolor != None:
        cax.set_facecolor(bgcolor) 
    else:
        cax.set_facecolor('white')

    at = AnchoredText(text, loc=10, frameon=False, prop=dict(fontweight='bold', rotation='vertical', 
                    size=size, color=color))
    cax.add_artist(at)


def bokeh_theme(return_color_list=True):
    theme_json =  {
    'attrs' : {
        'Figure' : {
            'background_fill_color': '#EEEEEE',
        },
        'Axis': {
            'axis_line_color': 'slategray',
            'major_tick_line_color': None,
            'minor_tick_line_color': None,
        },
        'Legend': {
            'border_line_color': 'slategray',
            'background_fill_color': '#EEEEEE',
            'border_line_width': 0.75,
            'background_fill_alpha': 0.75,
        },
        'Grid': {
            'grid_line_color': '#FFFFFF',
            'grid_line_width': 0.75,
        },
        'Text': {
            'text_font_style': 'italic',
            'text_font': 'Arial', 
            'text_font_size':10,
        },
        'Title': {
            'background_fill_color': '#EEEEEE',
            'text_color': '#3c3c3c',
            'align': 'center',
            'text_font': 'Arial',
            'offset': 2,
         }
    }
    }

    colors = {'dark_purple': '#5F2E88', 'dark_orange':'#F38227', 'black':'#444147',
        'dark_blue': '#3F60AC', 'dark_red':'#9C372F', 'dark_green':'#395A34',
        'purple': '#7E59A2', 'orange':'#E39943', 'blue': '#7292C7', 'red':'#C76A6A',
        'green':'#688A2F', 'light_purple':'#A17DB8', 'light_orange':'#EEBA7F',
        'light_blue':'#A5B3CC', 'light_red':'#E39C9D', 'light_green':'#B3CD86', 
        'grey': '#EFEFEF', 'gray': '#EFEFEF', 'light_grey':'#6D6F72'}
        
    color_items = [v for v in colors.values()]
    theme = bokeh.themes.Theme(json=theme_json)
    bokeh.io.curdoc().theme = theme
    if return_color_list:
        return [colors, color_items]
    else:
        return colors

def growth_animation(images, fname, contours=None, descriptors={'bar_length':10,
    'ip_dist':0.07, 'temperature':'', 'carbon':'', 'time_interval':0}):
    """
    Generates and saves a growth movie from supplied images and segmentation information.

    Parameters
    ----------
    images : list
        Image sequence to be animated
    fname: str
        Desired filename of animation.
    contours: list of arrays
        Detected contours in segemented image to be plotted over the images. If None,
        no contours will be generated.
    descriptors: dict
        List of sample descriptors to be printed on the image. Allowed descriptors are
            bar_length : float [µm]
                Length of scale bar in units of µm.
            ip_dist : float [µm per pixel]
                Interpixel distance of the image.
            temperature: float or str [C]
                Temperature of the experiment in degrees C.
            carbon: string
                Carbon source of the experiment.
            time_interval: float
                Interval between frames in units of minutes.
    """

    # Instantiate the figure.
    fig, ax = plt.subplots(1,1)
    def update(it):
        ax.clear()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

        # Set the appropriate labels.
        bar = ax.hlines(-30, 10, 10 + descriptors['bar_length']/descriptors['ip_dist'],
                        lw=3, color='k')
        bar_label = ax.text(10, -60, '{} µm'.format(descriptors['bar_length']), color='k',
                        fontsize=8)
        title = ax.set_title('{}° C {}'.format(descriptors['temperature'], descriptors['carbon']),
                        color='k', fontsize=8, y=0.97)
        time = ax.text(0.8, 0.99, '{} min'.format(it * descriptors['time_interval']),
                    transform=ax.transAxes, fontsize=8)
        # Plot the image.
        _ = ax.imshow(images[it], cmap=plt.cm.magma)

        # Determine if contours should be plotted.
        if contours != None:
            for c in contours:
                _ = ax.plot(c[:, 1], c[:, 0], color='lightskyblue', lw=1)
        plt.tight_layout()

    # Set the animator and save.
    anim = FuncAnimation(fig, update, frames=np.arange(0, len(images), 1),
                        interval=100)
    anim.save(fname, dpi=200, writer='imagemagick')


def fill_between(p, domain, range_1, range_2, **kwargs):
    """
    Creates a shaded area between two curves on a bokeh figure
    similar to matplotlibs `fill_between`.
    Parameters
    ----------
    p : Bokeh Figure
        Figure canvas on which to plot shaded region.
    domain : numpy array
        Range of x-values over which to shade
    range_1, range_2 : numpy array
        The two curves to fill between.
    """

    # Ensure that domain, range_1, and range_2 are the same length.
    lengths = np.unique([len(domain), len(range_1), len(range_2)])
    if len(lengths) > 1:
        raise RuntimeError(
            "domain, range_1, and range_2 must be the same length.")

    # Set up the patch objects.
    patch_x = np.append(domain, domain[::-1])
    patch_y = np.append(range_1, range_2[::-1])

    # Add a patch to the bokeh Figure.
    p.patch(patch_x, patch_y, **kwargs)
