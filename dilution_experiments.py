import sys
import mwc.viz
import mwc.stats
import mwc.io
import mwc.model
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import itertools
from bokeh.models import ColumnDataSource, Div
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup
from datetime import datetime
import glob
constants = mwc.model.load_constants()
import os

def dilution_forager(info_df):
    print(os.getcwd())
    # Compute the thoery line. 
    rep_range = np.logspace(0, 4, 200)
    theo = mwc.model.SimpleRepression(rep_range, ep_r=constants['O2'], ka=constants['ka'],
                                     ki=constants['ki'], ep_ai=constants['ep_ai'], 
                                      n_ns=constants['n_ns'], n_sites=constants['n_sites'],
                                     effector_conc=0).fold_change()
    # Define the data sources  
    fluct_source = ColumnDataSource(dict(x=[], y=[]))
    alpha_source = ColumnDataSource(dict(x=[], y=[]))
    fit_source = ColumnDataSource(dict(x=[], y1=[], y2=[]))
    binned_source = ColumnDataSource(dict(x=[], y=[]))
    mch_expression_source = ColumnDataSource(dict(x=[], y=[], c=[]))
    yfp_expression_source = ColumnDataSource(dict(x=[], y=[], c=[]))
    rep_atc_source = ColumnDataSource(dict(x=[], y=[], xs=[], ys=[]))
    fc_source = ColumnDataSource(dict(x=[], y=[]))
    
    
    
    # Instantiate the figure cavases
    p_fluct = bokeh.plotting.figure(width=400, height=400, 
                              x_axis_type='log', y_axis_type='log',
                              x_axis_label='summed intensity [a.u.]',
                              y_axis_label='squared fluctuations [a.u.^2]')
    p_alpha= bokeh.plotting.figure(width=400, height=400,
                                  x_axis_label='calibration factor [a.u. / molecule]',
                                  y_axis_label='sampling frequency')
    
    p_mch_expression = bokeh.plotting.figure(width=310, height=200,
                                  x_axis_label='mean mCherry pixel intensity [a.u.]',
                                  y_axis_label='ECDF')    
    p_yfp_expression = bokeh.plotting.figure(width=310, height=200,
                                  x_axis_label='mean YFP pixel intensity [a.u.]',
                                  y_axis_label='ECDF')     
    p_rep_atc = bokeh.plotting.figure(width=400, height=400,
                                  x_axis_label='ATC concentration [ng / mL]',
                                  y_axis_label='repressors per cell')
    
    p_fc = bokeh.plotting.figure(width=400, height=400,
                                 x_axis_label='repressors per cell',
                                 y_axis_label='fold-change',
                                 x_axis_type='log', y_axis_type='log')
    # Add starting glyphs
    p_fluct.circle(x='x', y='y', source=fluct_source, size=1, color='black')
    p_fluct.circle(x='x', y='y', source=binned_source, color='tomato', size=5)
    p_alpha.step(x='x', y='y', source=alpha_source, line_width=2, color='black')
    p_fluct.line(x='x', y='y1', source=fit_source, line_width=2, color='dodgerblue')
    p_fluct.line(x='x', y='y2', source=fit_source, line_width=2, color='dodgerblue')
    p_rep_atc.circle(x='x', y='y', source=rep_atc_source, color='blue')
    p_rep_atc.line(x='x', y='y', source=rep_atc_source, color='blue')
    p_rep_atc.multi_line(xs='xs', ys='ys', source=rep_atc_source, color='blue')
    p_mch_expression.circle(x='x', y='y', color='c', source=mch_expression_source, size=1)
    p_yfp_expression.circle(x='x', y='y', color='c', source=yfp_expression_source, size=1)
    
    _span_low = bokeh.models.Span(location=0, dimension='height', line_color='dodgerblue',
                                 line_width=3)
    _span_high = bokeh.models.Span(location=1, dimension='height', line_color='dodgerblue',
                                line_width=3)
    p_alpha.add_layout(_span_low)
    p_alpha.add_layout(_span_high)
    
    # Fold-change plotting. 
    p_fc.line(rep_range, theo, color='black', line_width=2)
    p_fc.circle(x='x', y='y', color='tomato', size=5, source=fc_source)
    
    # Define the selector callbacks
    def restrict_date(attr, old, new): 
        active = {0: 'accepted', 1:'rejected'} 
        date.options = list(info_df[(info_df['carbon']==new) &\
                                    (info_df['status']==active[status.active])]['date'].unique())
    def restrict_run(attr, old, new): 
        active = {0: 'accepted', 1:'rejected'} 
        run.options = list(info_df[(info_df['carbon']==carbon.value) &
                                   (info_df['date']==new)] &
                                   (info_df['status']==active[status.active])['run_number'].unique())

    # Define the selector buttons
    status = RadioButtonGroup(labels=['Accepted', 'Rejected'], active=0)
    carbon = Select(title='Carbon Source', value='glucose', options=list(info_df['carbon'].unique()))
    carbon.on_change('value', restrict_date)   
    date = Select(title='Date', value='20181002', options=list(info_df['date'].unique()))   
    carbon.on_change('value', restrict_date)   
    run = Select(title='Run Number', value='1.0', options=list(info_df['run_number'].unique().astype(str)))
    
    
    # Define the adjustment buttons
    agg_stat = RadioButtonGroup(labels=['Mean', 'Median', 'Mode'], active=0)
    bins = Slider(title='Events per Bin', value=50, start=5, end=500, step=5) 
    cred_range = Slider(title='% Credible Region', value=95, start=5, end=99, step=5)
    
    # Function to choose the experiment based off the selector values
    def select_experiment():
        # Get slider state
        carbon_val = carbon.value 
        date_val = date.value
        run_val = run.value
        
        # Define relative paths
        pref = f'{date_val}_r{int(float(run_val))}_37C_{carbon_val}'
        path = f'{pref}_O2_dilution/output/{pref}_O2_cal_factor_samples.csv'
        
        # Load data files
        samples = pd.read_csv(f'code/processing/microscopy/{path}')
        fluct_df = pd.read_csv(f'code/processing/microscopy/{pref}_O2_dilution/output/{pref}_O2_fluctuations.csv')
        fc_df = pd.read_csv(f'code/processing/microscopy/{pref}_O2_dilution/output/{pref}_O2_foldchange.csv')
        
        # Compute the initial hpd, assume 95%
        hpd_min, hpd_max = mwc.stats.compute_hpd(samples['alpha'], mass_frac=0.95)
        
        # Perform necessary data cleaning. 
        mean_auto = fc_df[fc_df['strain']=='auto']['mean_mCherry'].mean()
        fluct_df['I_1_sub'] = (fluct_df['I_1'] - mean_auto) * fluct_df['area_1']
        fluct_df['I_2_sub'] = (fluct_df['I_2'] - mean_auto) * fluct_df['area_2']
        fluct_df['summed'] = fluct_df['I_1_sub'] + fluct_df['I_2_sub']
        fluct_df['sq_fluct'] = (fluct_df['I_1_sub'] - fluct_df['I_2_sub'])**2
        
        hist, bins = np.histogram(samples['alpha'], bins=75)
        return [fluct_df, samples, {'bins':bins[:-1], 'hist':hist}, fc_df]
    
    # Adjust the source data to show the binnned data
    def _update_binned_events(fluct_df):
        binned = mwc.stats.bin_by_events(fluct_df, average=['summed', 'sq_fluct'], bin_size=bins.value)
        binned_source.data = dict(x=binned['summed'], y=binned['sq_fluct'])
                
    def _update_cred_plot(hpd_min, hpd_max, flucts):
        summed_range = np.logspace(np.log10(flucts['summed'].min() - 1), np.log10(flucts['summed'].max() + 1), 200)
        min_fit = hpd_min * summed_range
        max_fit = hpd_max * summed_range
        fit_source.data = dict(x=summed_range, y1=min_fit, y2=max_fit)
        return hpd_min, hpd_max
        
    def _update_cred_region(alpha_samples, flucts):
        hpd_min, hpd_max = mwc.stats.compute_hpd(alpha_samples, mass_frac=cred_range.value/100)
        _span_low.location = hpd_min
        _span_high.location = hpd_max
        _update_cred_plot(hpd_min, hpd_max, flucts)
        return hpd_min, hpd_max
    
    def _update_expression_distribution(fc_data): 
        fc_data = fc_data[fc_data['strain']=='dilution']
        mch_dfs = []
        yfp_dfs = []
        m_colors = bokeh.palettes.inferno(len(fc_data['atc_ngml'].unique()) + 3)
        y_colors = bokeh.palettes.viridis(len(fc_data['atc_ngml'].unique()) + 3) 
        i = 0
        for g, d in fc_data.groupby(['atc_ngml']):
            
            mch_x = np.sort(d['mean_mCherry'])
            yfp_x = np.sort(d['mean_yfp'])
            _y = np.arange(0, len(d), 1) / len(d) 
            _mdf = pd.DataFrame(np.array([mch_x,_y]).T, columns=['x', 'y'])
            _ydf = pd.DataFrame(np.array([yfp_x,_y]).T, columns=['x', 'y'])
            _mdf['c'] = m_colors[i]
            _ydf['c'] = y_colors[i]
            
            mch_dfs.append(_mdf)
            yfp_dfs.append(_ydf)
            i += 1
        mch_df = pd.concat(mch_dfs)
        yfp_df = pd.concat(yfp_dfs)
        mch_expression_source.data = dict(x=mch_df['x'], y=mch_df['y'], c=mch_df['c'])
        yfp_expression_source.data = dict(x=yfp_df['x'], y=yfp_df['y'], c=yfp_df['c'])
         
    def _update_atc_titration(fc_data, alpha_samps, hpd_min, hpd_max):
        _agg_stat = {0:np.median(alpha_samps['alpha']), 
                     1:np.mean(alpha_samps['alpha']),
                    2:alpha_samps.iloc[np.argmax(alpha_samps['log_prob'].values)]['alpha']}

        # Compute the repressors
        mean_auto = np.mean(fc_data[fc_data['strain']=='auto']['mean_mCherry'])
        fc_data = fc_data[fc_data['strain']=='dilution'].copy()
        fc_data['rep'] = ((fc_data['mean_mCherry'].values - mean_auto) * fc_data['area_pix'].values) / _agg_stat[agg_stat.active]
        fc_data['rep_min'] = (fc_data['mean_mCherry'].values - mean_auto) * fc_data['area_pix'].values / hpd_max
        fc_data['rep_max'] = (fc_data['mean_mCherry'].values - mean_auto) * fc_data['area_pix'].values / hpd_min
        _grouped = fc_data.groupby('atc_ngml').mean().reset_index()
        rep_atc_source.data = dict(x=_grouped['atc_ngml'], y=_grouped['rep'], xs=[(x, x) for x in _grouped['atc_ngml']],
                                   ys=[(ymin, ymax) for ymin, ymax in zip(_grouped['rep_min'], _grouped['rep_max'])])
        return _grouped
        
    def _update_foldchange(fc_data):
        fc_source.data = dict(x=fc_data['rep'], y=fc_data['fold_change'])
        
    # Main callback for updating the entire app
    def update():
        fluct_df, samples, step_dict, fc_df = select_experiment()
        _update_binned_events(fluct_df)
        _update_expression_distribution(fc_df)
        hpd_min, hpd_max = _update_cred_region(samples['alpha'], fluct_df)
        summarized_fc = _update_atc_titration(fc_df, samples, hpd_min, hpd_max)
        _update_foldchange(summarized_fc)
        fluct_source.data = dict(x=fluct_df['summed'], y=fluct_df['sq_fluct'])
        alpha_source.data = dict(x=step_dict['bins'], y=step_dict['hist'])
         
    controls = [status, carbon, date, run, bins, cred_range, agg_stat]
    for control in controls[1:-1]:
        control.on_change('value', lambda attr, old, new: update()) 
    for control in [status, agg_stat]:
        control.on_change('active', lambda attr, old, new: update())
    inputs = widgetbox(*controls, sizing_mode='scale_width')
    _lay_expression = bokeh.layouts.Column(p_mch_expression, p_yfp_expression)
    lay = layout([[inputs, p_fluct, p_alpha], [_lay_expression, p_rep_atc, p_fc]], sizing_mode='fixed')
    tab = bokeh.models.widgets.Panel(child=lay, title='Dilution Experiments')
    return tab

