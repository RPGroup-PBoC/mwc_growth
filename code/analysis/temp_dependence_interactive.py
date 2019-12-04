#%%
import numpy as np
import pandas as pd
import phd.viz
import bokeh.io
import numpy as np
import bokeh.io
import bokeh.plotting
import bokeh.layouts
import bokeh.palettes
from bokeh.models import *
colors, palette = phd.viz.bokeh_theme()
bokeh.io.output_file('temp_dependence.html')
# Load the data
data = pd.read_csv('../../data/analyzed_foldchange.csv')
data = data[(data['repressors'] > 0) & (data['fold_change'] >= 0) &
            (data['size']=='large') & (data['temp']  != 37)]
data = data.groupby(['date', 'run_number', 'temp', 'atc_ngml']).mean().reset_index()
data = data.groupby(['temp', 'atc_ngml']).mean().reset_index()


# Define the sliders. 
slope_epRA = Slider(title='∆S_R log10 slope', start=-10, end=-2, step=0.0001, value=0)
slope_epAI = Slider(title='∆S_AI log10 slope', start=-10, end=-2, step=0.0001, value=0)


cb_code = """
   var temp = source.data['temp'];
   var reps = source.data['reps'];
   var epRA = -13.9;
   var epAI = 4.5;
   var Nns = 4600000;
   var temp_ref = 37 + 273.15

   for (var i = 1; i < temp.length; i++) {
       var delS_RA = Math.pow(10, slope_epRA.value) * temp[i];
       var delS_AI = Math.pow(10, slope_epAI.value) * temp[i];
       var epRA_star = delS_RA * (temp_ref - temp[i]) + epRA;
       var epAI_star = delS_AI * (temp_ref - temp[i]) + epAI;
       var pact = Math.pow(1 + Math.exp(-epAI_star), -1);
       var fc = Math.pow(1 + pact * (reps[i] / Nns) * Math.exp(-epRA_star), -1);
       source.data['fc'][i] = fc;
    }
   source.change.emit()
"""

# define the data source. 
df1 = pd.DataFrame([], columns=['reps', 'temp', 'fc'])
df2 = pd.DataFrame([], columns=['reps', 'temp', 'fc'])
rep_range = np.logspace(0, 4, 800)
df1['reps'] = rep_range
df1['temp'] = 32 + 273.15
df1['fc'] = np.ones(len(rep_range))
df1['color'] = colors['blue']
df2['reps'] = rep_range
df2['temp'] = 42 + 273.15
df2['fc'] = np.ones(len(rep_range))
df2['color'] = colors['red']
df = pd.concat([df1, df2])
source = ColumnDataSource(df)
                          

# Define the callback. 
cb = CustomJS(args={'source':source, 'slope_epRA':slope_epRA, 
            'slope_epAI':slope_epAI}, code=cb_code)

# Set up the figure canvas
p = bokeh.plotting.figure(width=500, height=500, x_axis_type='log', y_axis_type='log',
                        x_range=[1, 1E3], y_range=[1E-3, 1.1])
_colors = {32:colors['blue'], 42:colors['red']}
for g, d in data.groupby(['temp']):
    p.circle(d['repressors'], d['fold_change'], color=_colors[g],
            legend_label=f'{g} °C', size=9)

slope_epRA.js_on_change('value', cb)
slope_epAI.js_on_change('value', cb)

# Plot the lines. 
col = bokeh.layouts.column(slope_epRA, slope_epAI)
p.circle(x='reps', y='fc', color='color', source=source, size=2)
lay = bokeh.layouts.row(p, col)
bokeh.io.save(lay)


# %%
