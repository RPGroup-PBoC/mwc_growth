# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import mwc.viz
import bokeh.io
import bokeh.layouts
import bokeh.plotting
from bokeh.models import Dropdown, RadioButtonGroup, ColumnDataSource, Slider
bokeh.io.output_notebook()
colors, colors_list = mwc.viz.bokeh_theme()


# Load the fluctuation data
flucts = pd.read_csv('../../data/analyzed_fluctuations.csv')

def dilution_doc(doc):
    # Define the interactions
    button = RadioButtonGroup(labels=['carbon source variation', 'temperature variation'], active=0)
    drop1 = Dropdown(label='select carbon source', menu=[('glucose', 'glucose'), ('glycerol', 'glycerol'), ('acetate', 'acetate')], value='glucose')
    drop2 = Dropdown(label='select temperature', menu=[('42° C', '42'), ('37° C', '37'), ('32° C', '32')], value='37')
    bins = Slider(value=50, start=5, end=1000, step=1, title='events per bin', bar_color=colors['purple'])

    # Define the figure canvas and layouts
    p = bokeh.plotting.figure(x_axis_type='log',y_axis_type='log', 
                             x_axis_label='I\u2081 + I\u2082',
                             y_axis_label='(I\u2081 - I\u2082)\u00b2',
                             width=600, height=400)
    row = bokeh.layouts.row(drop1, drop2)
    lay = bokeh.layouts.column(button, row, bins, p)

    # Define the source
    fluct_source = ColumnDataSource(dict(summed=[], fluct=[]))
    bin_source = ColumnDataSource(dict(I_range=[], fit=[]))
    fit_source = ColumnDataSource(dict(binned_summed=[], binned_fluct=[]))
    I_tot_range = np.logspace(2, 6, 200)
    def update():
        if button.active == 0:
            temp = 37
            carbon = drop1.value

        elif button.active == 1:
            temp = int(drop2.value)
            carbon = 'glucose'

        # Isolate the dataset 
        _df = flucts[(flucts['carbon']==carbon) & (flucts['temp']==temp)]
        alpha = _df['alpha_mu'].values[0]

        # Get the number of bins
        binning = mwc.stats.bin_by_events(_df, bins.value)
        p.title.text = f'{carbon}, {temp}°C, no antibiotic'
        fluct_source.data = dict(summed=_df['summed'], fluct=_df['fluct'])
        fit_source.data = dict(I_range=I_tot_range, fit=alpha * I_tot_range)
        bin_source.data = dict(binned_summed=binning['summed'], binned_fluct=binning['fluct'])
        bokeh.io.push_notebook()

    p.circle(x='summed', y='fluct', color='black', size=1, alpha=0.5, legend='division', source=fluct_source)
    p.circle(x='binned_summed', y='binned_fluct', line_color=colors['purple'], fill_color=colors['light_purple'], legend='binned divisions', source=bin_source)
    p.line(x='I_range', y='fit', color=colors['orange'], line_width=2, legend='best fit', source=fit_source)
    p.legend.location = 'top_left'
    for v in [drop1, drop2, bins]:
        v.on_change('value', lambda attr, old, new: update())
    update()
    doc.add_root(lay)

bokeh.io.show(dilution_doc)



#%%


#%%
