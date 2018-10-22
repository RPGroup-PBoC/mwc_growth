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
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.signal import gaussian, convolve

# Thematic settings
def set_style(theme=''):
    theme = bokeh.themes.theme.Theme(
        json={
        'attrs': {
            'Figure': {
                'background_fill_color': '#EEEEEE'},
            'Grid': {
                'grid_line_width': 1,
                'grid_line_color': '#ffffff'},
        'Text' : {
            'text_font': "Open Sans"},
        'Axis': {
            'axis_label_text_font': "Helvetica",
            'axis_label_text_font_style': "normal",
            'major_label_text_font': "Helvetica"}
        }
        })
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

# Specialized viewing functions.

def bokeh_traceplot(samples, varnames=None):
    """
    Generate a Bokey traceplot of a series of parameters and their sampling 
    history. As of now, the only shown distribution is 'ecdf'. 

    Parameters
    ----------
    samples: StanFit4Model object
        Sampling output from Stan.
    varnames: list
        List of variable names to include in the plot.
    """
    params = samples.model_pars
    sample_values = samples.extract()
    palette = bokeh.palettes.Category10_10
    
    # Define the standard plot sizes.   
    pairs = []
    if varnames != None:
        iterator = varnames
    else:
        iterator = params
    for p in iterator:
        colors = itertools.cycle(palette)
        if len(np.shape(sample_values[p])) == 1:
            _samples = np.array([sample_values[p]]).T
        else:
            _samples = sample_values[p]
      
        dfs = []
        trace = bokeh.plotting.figure(plot_width=400, plot_height=200, 
                                      x_axis_label='sample number', y_axis_label=f'{p}',
                                     title=f'sampling history for {p}', background_fill_color='#ecedef') 

        dist = bokeh.plotting.figure(plot_width=400, plot_height=200, 
                                     x_axis_label=f'{p}', y_axis_label='ECDF',
                                    title=f'posterior distribution for {p}', background_fill_color='#ecedef')   
        
        # Add visual formatting
        trace.xgrid.grid_line_color = '#FFFFFF'
        trace.ygrid.grid_line_color = '#FFFFFF'
        dist.xgrid.grid_line_color = '#FFFFFF'
        dist.ygrid.grid_line_color = '#FFFFFF'    
        
        for i, color in zip(range(np.shape(_samples)[1]), colors): 
            # Extract the specific values. 
            _values = _samples[:, i]
            x, y = np.sort(_values), np.arange(0, len(_values), 1) / len(_values)
            _df = pd.DataFrame(np.array([x, y, 
                                         _values, np.arange(0, len(_values), 1)]).T, 
                               columns=['x', 'ecdf', 'samples', 'step_no'])            
            dist.line(_df['x'], _df['ecdf'], line_width=2, color=color)
            trace.line(_df['step_no'], _df['samples'], line_width=1, color=color) 
        pairs.append([dist, trace])  
    return bokeh.io.show(bokeh.layouts.gridplot(pairs))

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
