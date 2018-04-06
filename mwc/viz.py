import bokeh.io
import bokeh.plotting
import bokeh.layouts
import bokeh.palettes
import skimage.io
import skimage.measure
import seaborn as sns
import numpy as np

import os
import matplotlib.pyplot as plt
import matplotlib
from scipy.signal import gaussian, convolve


def pub_style(return_colors=True):
    """
    Sets the style to the publication style

    Parameters
    ----------
    return_colors: Bool
        If True, this will also return a dictionary

    """
    rc = {'axes.facecolor': '#E3DCD0',
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': '-',
          'grid.linewidth': 0.5,
          'grid.alpha': 0.75,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.facecolor': '#FFEDCE',
          'figure.dpi': 150}

    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)
    colors = ['#D56C55', '#738FC1', '#7AA974', '#D4C2D9', '#F1D4C9', '#C9D7EE',
              '#DCECCB']
    if return_colors:
        return colors


def personal_style():
    """
    Sets the plotting style to my personal preference.
    """
    rc = {'axes.facecolor': '#EAECEE',
          'font.family': 'Arial',
          'font.style': 'italic',
          'axes.grid': True,
          'axes.edgecolor': 'slategray',
          'axes.spines.right': False,
          'axes.spines.top': False,
          'axes.axisbelow': True,
          'axes.linewidth': 0.75,
          'axes.titlesize': 8,
          'axes.grid': True,
          'lines.linewidth': 1.2,
          'lines.dash_capstyle': 'round',
          'grid.linestyle': '-',
          'grid.linewidth': 0.75,
          'grid.alpha': 0.5,
          'grid.color': '#B9BBBC',
          'axes.labelsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'legend.fontsize': 8,
          'legend.frameon': False,
          'xtick.color': '#4b4b4b',
          'ytick.color': '#4b4b4b',
          'axes.xmargin': 0.01,
          'axes.ymargin': 0.01,
          'figure.dpi': 150}

    # plt.rc('mathtext', fontset='dejavuserif', sf='sans')
    plt.rc('text.latex', preamble=r'\usepackage{mathpazo}')
    matplotlib.style.use(rc)
    flat = ['#34495e', '#c0392b', '#3498db', '#27ae60', '#7B1FA2', '#d35400']
    sns.set_palette(flat)
    return flat


def format_axes(pub_style=False):
    """
    Executes a seaborn despining function with my prefered offset and trimming.
    """
    if pub_style == False:
        sns.despine(offset=7)


# Data specific plotting functions.
def generate_flow_summary(df, file, savedir='./output', pub_style='False'):
    """
    Generates a summary plot for flow cytometry experiments.

    Parameters
    ----------
    df : pandas DataFrame
        DataFrame containing all measured events after gating.
    file : str
        Converted .csv file of raw flow cytometry measurement with gating
        identifier. This will be used to make the example flow cloud.
    savedir : str
        Target directory for saved figure.
    pub_style : bool
        If true, the axes will not be formatted to match the publication
        style.

    Returns
    -------
    fig : matplotlib Figure canvas
        The generated figure canvas.
    """
