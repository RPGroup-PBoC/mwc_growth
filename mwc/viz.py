import bokeh.io
import bokeh.plotting
import bokeh.layouts
import bokeh.palettes
import skimage.io
import skimage.measure
import seaborn as sns
import numpy as np
import cairosvg
import os
import matplotlib.pyplot as plt
from scipy.signal import gaussian, convolve

# ---------------------------------------------------------------------------
# Plotting styles
# ---------------------------------------------------------------------------


def set_plotting_style(return_colors=True):
    """
    Sets the plotting style.

    Parameters
    ----------
    return_colors: Bool
        If True, this will also return a palette of eight color-blind safe
        colors with the hideous yellow replaced by 'dusty purple.'
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
    # colors = sns.color_palette('colorblind', n_colors=8)
    # colors[4] = sns.xkcd_palette(['dusty purple'])[0]
    colors = {'green': '#7AA974', 'light_green': '#BFD598',
              'pale_green': '#DCECCB', 'yellow': '#EAC264',
              'light_yellow': '#F3DAA9', 'pale_yellow': '#FFEDCE',
              'blue': '#738FC1', 'light_blue': '#A9BFE3',
              'pale_blue': '#C9D7EE', 'red': '#D56C55', 'light_red': '#E8B19D',
              'pale_red': '#F1D4C9', 'purple': '#AB85AC',
              'light_purple': '#D4C2D9'}
    if return_colors:
        return colors


def boilerplate(**kwargs):
    """
    A boiler-plate for a bokeh plotting figure.
    """
    # Make a bokeh figure axis.
    if kwargs is not None:
        p = bokeh.plotting.figure(**kwargs)
    else:
        p = bokeh.plotting.figure()

    # Apply the styling to the figure axis.
    p.background_fill_color = '#E3DCD0'
    p.grid.grid_line_color = '#FFFFFF'
    p.grid.grid_line_width = 0.75
    p.axis.minor_tick_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.axis_line_color = None
    p.axis.axis_label_text_font = 'Lucida Sans Unicode'
    p.axis.major_label_text_font = 'Lucida Sans Unicode'
    p.axis.axis_label_text_font_style = 'normal'
    p.axis.axis_label_text_font_size = '1em'
    p.axis.major_label_text_font_size = '0.75em'
    p.axis.axis_label_text_color = '#3c3c3c'
    p.axis.axis_label_standoff = 3
    return p


def bokeh_to_pdf(p, fname):
    """
    Save a bokeh figure to a pdf using cairosvg.

    Parameters
    ----------
    p : bokeh figure object
        The handle of the bokeh figure to be saved.
    fname: str
        The file name of the resulting object.
    """
    # Check if temp file exists.
    if os.path.exists('._tmp.svg'):
        os.remove('._tmp.svg')

    # Export temporarily to svg.
    bokeh.io.export_svgs(p, '._tmp.svg')

    # Convert svg object to pdf
    cairosvg.svg2pdf(url='./._tmp.svg', write_to=fname)

    # Remove temporary file.
    os.remove('./._tmp.svg')
    print('Saved bokeh figure as {0}'.format(fname))


def imshow(im, color_mapper=None, plot_height=400, length_units='pixels',
           interpixel_distance=1.0, return_glyph=False):
    """
    Display an image in a Bokeh figure.

    Parameters
    ----------
    im : 2-dimensional Numpy array
        Intensity image to be displayed.
    color_mapper : bokeh.models.LinearColorMapper instance, default None
        Mapping of intensity to color. Default is 256-level Viridis.
    plot_height : int
        Height of the plot in pixels. The width is scaled so that the
        x and y distance between pixels is the same.
    length_units : str, default 'pixels'
        The units of length in the image.
    interpixel_distance : float, default 1.0
        Interpixel distance in units of `length_units`.
    return_glyph : book, default False
        If True, the image GlyphRenderer is also returned.

    Returns
    -------
    output : bokeh.plotting.figure instance
        Bokeh plot with image displayed.

    Notes
    -----
    We thank Justin Bois for writing this function. http://bebi103.caltech.edu
    """
    # Get shape, dimensions
    n, m = im.shape
    dw = m * interpixel_distance
    dh = n * interpixel_distance

    # Set up figure with appropriate dimensions
    plot_width = int(m / n * plot_height)
    kwargs = {'plot_height': plot_height, 'plot_width': plot_width,
              'x_range': [0, dw], 'y_range': [0, dh],
              'x_axis_label': length_units, 'y_axis_label': length_units,
              'tools': 'pan, box_zoom, wheel_zoom, reset'}
    p = boilerplate(**kwargs)

    # Set color mapper; we'll do Viridis with 256 levels by default
    if color_mapper is None or color_mapper is 'viridis':
        color_mapper = bokeh.models.LinearColorMapper(
            bokeh.palettes.viridis(256))
    if color_mapper is 'Greys_r':
        color_mapper = bokeh.models.LinearColorMapper(
            bokeh.palettes.Greys_r(256))

    # Display the image
    im_bokeh = p.image(image=[im[::-1, :]], x=0, y=0, dw=dw, dh=dh,
                       color_mapper=color_mapper)

    if return_glyph is True:
        return p, im_bokeh
    else:
        return p


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
