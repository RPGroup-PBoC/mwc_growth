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

# Thematic settings
def set_style(holoviews=False):
    theme = bokeh.themes.theme.Theme(
        json={
        'attrs': {
            'Figure': {
                'background_fill_color': '#EEEEEE'},
            'Grid': {
                'grid_line_width': 1,
                'grid_line_color': '#ffffff'},
        'Text' : {
            'text_font': "Cabin"},
        'Axis': {
            'axis_label_text_font': "Open Sans",
            'axis_label_text_font_style': "normal",
            'major_label_text_font': "Open Sans"}
        }
        })
    if holoviews == True:
        hv.renderer('bokeh').theme = theme

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
    Sets the plotting style to my preference
    """ 
    rc = {'axes.facecolor': '#f1f2f6',  #'#EAECEE',
          'font.family': 'sans-serif',
          'font.style': 'italic',
          'axes.grid': True,
          'axes.edgecolor': 'slategray',
          'axes.spines.right': False,
          'axes.spines.top': False,
          'axes.axisbelow': True,
          'axes.linewidth': 1,
          'axes.titlesize': 8,
          'axes.grid': True,
          'lines.linewidth': 2,
          'lines.dash_capstyle': 'round',
          'grid.linestyle': '-',
          'grid.linewidth': 0.75,
          'grid.alpha': 0.5,
          'grid.color': '#B9BBBC',
          'axes.labelsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'legend.fontsize': 8,
          'legend.frameon': True,
          'xtick.color': '#4b4b4b',
          'ytick.color': '#4b4b4b',
          'axes.xmargin': 0.01,
          'axes.ymargin': 0.01,
          'figure.dpi': 100}

    # plt.rc('mathtext', fontset='dejavuserif', sf='sans')
    plt.rc('text.latex', preamble=r'\usepackage{mathpazo}')
    matplotlib.style.use(rc)
    flat = ['#64767C', '#484B3E','#95B7D8', '#699FCE','#6B5E86','#8389B4',  '#A6DCE8', '#72A2B6',
            '#6D7960']
    sns.set_palette(flat)
    return flat      


def pboc_style(grid=True):
    """
    Sets the style to the publication style
    """
    rc = {'axes.facecolor': '#E3DCD0',
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': '-',
          'grid.linewidth': 0.5,
          'grid.alpha': 0.75,
          'grid.color': '#ffffff',
          'axes.grid': grid,
          'ytick.direction': 'in',
          'xtick.direction': 'in',
          'xtick.gridOn': True,
          'ytick.gridOn': True,
          'ytick.major.width':5,
          'xtick.major.width':5,
          'ytick.major.size': 5,
          'xtick.major.size': 5,
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.facecolor': '#FFEDCE',
          'figure.dpi': 150,
           'xtick.color': 'k',
           'ytick.color': 'k'}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)


def phd_style():
    """
    Sets the plotting style to my preference
    """   
    rc = {'axes.facecolor': '#EFEFEF', # '#F5F9FC', #E5E8EA', #DFE8EF', #EAEAEA', #E0E1E2', 
          'font.family': 'sans-serif',
          'font.style': 'italic',
          'font.weight': 400,
          'font.family': 'Arial',
          'axes.edgecolor': 'slategray',
          'axes.spines.right': False,
          'axes.spines.top': False,
          'axes.axisbelow': True,
          'axes.linewidth': 0.75,
          'axes.titlesize': 6,
          'axes.grid': True,
          'lines.linewidth': 0.75,
          'lines.dash_capstyle': 'round',
          'grid.linestyle': '-',
          'grid.linewidth': 0.35,
          'grid.color': '#ffffff',
          'axes.labelsize': 6,
          'xtick.labelsize': 5,
          'ytick.labelsize': 5,
          'legend.fontsize': 6,
          'legend.frameon': True,
          'axes.xmargin': 0.02,
          'axes.ymargin': 0.02,
          'figure.dpi': 200}

    plt.rc('text.latex', preamble=r'\usepackage{mathpazo}')
    matplotlib.style.use(rc)
    colors = {'dark_purple': '#5F2E88', 'dark_orange':'#F38227', 'black':'#444147',
              'dark_blue': '#3F60AC', 'dark_red':'#9C372F', 'dark_green':'#395A34',
              'purple': '#7E59A2', 'orange':'#E39943', 'blue': '#7292C7', 'red':'#C76A6A',
               'green':'#688A2F', 'light_purple':'#A17DB8', 'light_orange':'#EEBA7F',
               'light_blue':'#A5B3CC', 'light_red':'#E39C9D', 'light_green':'#B3CD86', 
               'grey': '#EFEFEF', 'gray': '#EFEFEF', 'light_grey':'#6D6F72'}
    palette = [v for k, v in colors.items() if v not in ['grey', 'gray', 'dark_purple', 'light_grey']]
    sns.set_palette(palette)
    return colors     

def bokeh_theme(return_color_list=True):
    theme_json =  {
    'attrs' : {
        'Figure' : {
            'background_fill_color': '#EEEEEE',
        },
        'Axis': {
            'axis_line_color': 'black',
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
            'text_font': 'Stone Sans', 
            'text_font_size':10,
        },
        'Title': {
            'text_color': 'black',
            'align': 'left',
            'text_font_style': 'italic',
            'text_font': 'Helvetica',
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

def format_axes(pub_style=False):
    """
    Executes a seaborn despining function with my prefered offset and trimming.
    """
    if pub_style == False:
        sns.despine(offset=7)


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