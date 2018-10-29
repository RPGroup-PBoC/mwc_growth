# -*- coding: utf-8 -*-    
import sys
import mwc.io
from dilution_experiments import  dilution_forager
import pandas as pd
import bokeh.io
import bokeh.plotting
from bokeh.models import ColumnDataSource, Div
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup
import glob

 # Load all microscopy data and include determination of accepted or rejected.
experiments = glob.glob('code/processing/microscopy/*dilution')
print(experiments)
info_df = pd.DataFrame([])
for exp in experiments:
    date, run_number, temp, carbon, operator, _ = exp.split('/')[-1].split('_')
    status = mwc.io.scrape_frontmatter(exp)
    data = {'date':date, 'run_number':int(run_number.split('r')[1]), 
            'temp': int(temp.split('C')[0]), 'carbon':carbon, 'operator':operator,
            'status':status['status'].lower(), 'reason':status['reason']} 
    info_df = info_df.append(data, ignore_index=True)

tab1 = dilution_forager(info_df)
tabs = bokeh.models.widgets.Tabs(tabs=[tab1])
bokeh.io.curdoc().add_root(tabs)


