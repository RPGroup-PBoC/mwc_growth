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
from matplotlib.animation import FuncAnimation

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


# Generate the movie.
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
