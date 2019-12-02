# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
from .stats import compute_hpd, bin_by_events
from .model import SimpleRepression
import bokeh.io
import bokeh.plotting
import bokeh.layouts

def dilution_summary(fluct_df, mean_auto, alpha_samples, fname=None, ip_dist=0.065,
                    title=None):
    """Generates a Bokeh plot summarzing the dilution data."""
    # Compute necessary properties.
    fluct_df['I_1_sub'] = fluct_df['I_1'] - mean_auto
    fluct_df['I_2_sub'] = fluct_df['I_2'] - mean_auto
    
    # Determine if the fluctuations have been calculated
    if 'summed' not in fluct_df.keys():
        fluct_df['summed'] = fluct_df['I_1_sub'] * fluct_df['area_1'] +\
                        fluct_df['I_2_sub'] * fluct_df['area_2']
    if 'fluct' not in fluct_df.keys():
        fluct_df['fluct'] = (fluct_df['I_1_sub'] * fluct_df['area_1'] -\
                       fluct_df['I_2_sub'] * fluct_df['area_2'])**2
    
    # Compute the eCDF of fractional intensity. 
    frac_int = (fluct_df['I_1_sub'] * fluct_df['area_1']) / fluct_df['summed']
    x, ecdf = np.sort(frac_int), np.arange(0, len(frac_int), 1) / len(frac_int)
    
    # Bin the data by events. 
    binned = bin_by_events(fluct_df, 50)
    
    # Generate the histogram of alpha samples
    hist, bins = np.histogram(alpha_samples, bins=50, density=True)
 
    # Compute the fit line
    hpd = compute_hpd(alpha_samples, 0.95)
    summed_range = np.logspace(np.log10(np.min(fluct_df['summed'])) - 1,
                              np.log10(np.max(fluct_df['summed'])) + 1, 200)
    fit_lower = hpd[0] * summed_range
    fit_upper = hpd[1] * summed_range
    
    # Instantiate the figure canvas
    p1 = bokeh.plotting.figure(width=700, height=300,
                               x_axis_type='log', y_axis_type='log',
                               x_axis_label='summed intensity [a.u.]',
                               y_axis_label='fluctuations [a.u.^2]')
    p2 = bokeh.plotting.figure(width=300, height=200,
                               x_axis_label='cell area [µm^2]', 
                               y_axis_label='mean intensity [a.u. / pix]') 
    p3 = bokeh.plotting.figure(width=200, height=200, 
                               x_axis_label='partitioned intensity', 
                               y_axis_label='ECDF')
    p4 = bokeh.plotting.figure(width=200, height=200,
                              x_axis_label='α [a.u. / molecule]',
                              y_axis_label='~ probability')
    
    # Determine if the title should be added
    if title != None:
        p1.title.text = title

    # Plot the fluctuation glyphs
    p1.circle(fluct_df['summed'], fluct_df['fluct'], size=1, color='black', legend='divisions')
    p1.circle(binned['summed'], binned['fluct'], size=4, color='tomato', legend='binned data (50 events)')
    p1.line(summed_range, fit_upper, color='dodgerblue')
    p1.line(summed_range, fit_lower, color='dodgerblue')
   
    # Plot the intensity versus area
    p2.circle(fluct_df['area_1'] * ip_dist**2, fluct_df['I_1'], size=1.5, color='black', legend='I_1')
    p2.circle(fluct_df['area_2'] * ip_dist**2, fluct_df['I_2'], size=1.5, color='tomato', legend='I_2')
    
    # Plot the ecdf of partitioned intensity
    p3.line(x, ecdf, color='black', line_width=2)
    
    # Plot the cal factor distribution. 
    p4.step(bins[:-1], hist, color='black', line_width=2)
    
    # Create the layout and return
    _layout = bokeh.layouts.gridplot([[p1], [p2, p3, p4]])
     
    # Save the figure as html and png. 
    if fname != None:
        bokeh.io.export_png(_layout, f"output/{fname.split('.')[0]}.png")
        bokeh.io.save(_layout, f"output/{fname.split('.')[0]}.html", 
                    title=fname.split('.')[0])
    return _layout

def fc_summary_microscopy(fc_data, alpha_samples, fname=None, constants=None, 
                          operator='O2', title=None):
    """Generates a summary plot of measured fold-change"""    
    # Compute the credible regions and median for repressors per cell
    mean_auto = fc_data[fc_data['strain']=='auto']['mean_mCherry'].mean()
    hpd = compute_hpd(alpha_samples['alpha'], 0.95)
    int_mCherry = (fc_data['mean_mCherry'] - mean_auto) * fc_data['area_pix'] 
    fc_data['repressors_min'] = int_mCherry / hpd[0]
    fc_data['repressors_max'] = int_mCherry / hpd[1]
    fc_data['repressors_median'] = int_mCherry / np.median(alpha_samples['alpha'])

    # Compute the theory line
    rep_range = np.logspace(0, 4, 200)
    arch = SimpleRepression(rep_range, ep_r=constants[operator], ep_ai=constants['ep_ai'],
                                 n_ns=constants['n_ns'], ka=constants['ka'], ki=constants['ki'],
                                 n_sites=constants['n_sites'], effector_conc=0).fold_change()

    # Define the color palettes
    mch_colors = {conc: color for conc, color in zip(np.sort(fc_data['atc_ngml'].unique()), bokeh.palettes.Reds9)}
    yfp_colors = {conc: color for conc, color in zip(np.sort(fc_data['atc_ngml'].unique()), bokeh.palettes.Greens9)}

    # Instantiate the figure canvases
    p1 = bokeh.plotting.figure(width=350, height=250, 
                      x_axis_label='mean mCherry intensity [a.u. / pix]',
                      y_axis_label='mean YFP intensity [a.u. / pix]', x_axis_type='log', y_axis_type='log')
    p2 = bokeh.plotting.figure(width=350, height=250, 
                      x_axis_label='mean mCherry intensity [a.u. / pix]',
                      y_axis_label='ECDF', x_axis_type='log')
    p3 = bokeh.plotting.figure(width=350, height=250, x_axis_label='ATC [ng / mL]',
                          y_axis_label='rep. per cell')
    p4 = bokeh.plotting.figure(width=350, height=250, x_axis_type='log', y_axis_type='log',
                           x_axis_label='rep. per cell', y_axis_label='fold-change')

    # Add title if desired
    if title != None:
        p1.title.text = title

    for g, d in fc_data.groupby(['atc_ngml','strain']):
        # Compute the ECDF of mCherry intensity
        x_mch, y_mch = np.sort(d['mean_mCherry']), np.arange(0, len(d), 1) /len(d)
        if g[0] == 0: legend = g[1]
        else: legend = str(g[0])

        # Plot the various glyphs
        p1.circle(d['mean_mCherry'], d['mean_yfp'], color=mch_colors[g[0]], size=0.75, legend=legend)
        p2.line(x_mch, y_mch, color=mch_colors[g[0]], line_width=2, legend=legend)
        p3.circle(g[0], np.mean(d['repressors_median']))
        p3.line([g[0],g[0]],  [np.mean(d['repressors_min']), np.mean(d['repressors_max'])])
        p4.circle(np.mean(d['repressors_median']), np.mean(d['fold_change']))
        p4.line(rep_range, arch, line_width=2, color='black')
    
    p1.legend.label_text_font_size = '8pt'
    p2.legend.label_text_font_size = '8pt'
    p1.legend.spacing = -4
    p2.legend.spacing = -4
    p1.legend.click_policy = 'hide'
    p2.legend.click_policy = 'hide'
    
    # Assemble the layout and determine if the plot should be saved
    _layout = bokeh.layouts.gridplot([[p1, p2], [p3, p4]])
    if fname != None:
        bokeh.io.export_png(_layout, f"output/{fname.split('.')[0]}.png")
        bokeh.io.save(_layout, f"output/{fname.split('.')[0]}.html", title=fname)
    return _layout 

def _generate_summary(self):
    """Generates a validation summary string"""
    failures = self.failures
    summary = """
Test: `{}` ({} of {} passed)
=============================================================================
The following files are in error:\n
""".format(self.name, self.n_passed, self.len)

    if len(failures) == 0:
        failures = ['None\n']
    for f in failures:
        summary += f

    summary += """
-----------------------------------------------------------------------------
"""
    return summary


class MCMCValidation():
    def __init__(self):
        return True

    def test_convergence():
        return True

    def test_real_values():
        return True

    def run_suite():
        return True


class FlowValidation():
    """
    Validation test suite for flow cytometry data sets.
    """

    def __init__(self, fnames=None):
        """
        """
        if type(fnames) is not list:
            fnames = list(fnames)
        self.fnames = fnames
        self.dfs = [pd.read_csv(f, comment='#') for f in fnames]

    def test_flow_cols(self):
        """Ensures that csv file has appropriate columns."""
        self.name = 'test_flow_cols'
        fnames = self.fnames
        dfs = self.dfs
        cols = ['FSC-A', 'FSC-H', 'SSC-A', 'SSC-H', 'FITC-A', 'FITC-H', 'gate']
        failed = []
        n_passed = 0
        for d, f in zip(dfs, fnames):
            if list(d.columns) == cols:
                n_passed += 1
            else:
                failed.append(f.split('/')[-1] + '\n')

        self.n_passed = n_passed
        self.len = len(dfs)
        self.failures = failed
        return _generate_summary(self)

    def test_gate(self):
        """Ensures gated column is either 1 or 0."""
        self.name = 'test_gate'
        dfs = self.dfs
        fnames = self.fnames
        failed = []
        n_passed = 0
        for f, d in zip(fnames, dfs):
            ids = len(d[(d['gate'] == 0) | (d['gate'] == 1)])
            if ids == len(d):
                n_passed += 1
            else:
                failed.append(f.split('/')[-1] + '\n')
        self.n_passed = n_passed
        self.len = len(dfs)
        self.failures = failed
        return _generate_summary(self)

    def test_positivity(self):
        """Ensures < 1 percent of data is negative."""
        self.name = 'test_positivity'
        dfs = self.dfs
        fnames = self.fnames
        failed = []
        n_passed = 0
        for i, d in enumerate(dfs):
            perc_neg = (d[d['gate'] == 1]['FITC-A'] < 0).sum() / len(d)
            if perc_neg >= 0.01:
                failed.append(fnames[i].split('/')[-1] +
                              ' -> {} %\n'.format(perc_neg * 100))
            else:
                n_passed += 1

        self.n_passed = n_passed
        self.len = len(dfs)
        self.failures = failed
        return _generate_summary(self)

    def test_event_count(self, n_events=1E5):
        """Ensures there are at least a given number of events in the sample."""
        self.name = 'test_event_count'
        dfs = self.dfs
        fnames = self.fnames
        failed = []
        n_passed = 0
        for f, d in zip(fnames, dfs):
            if len(d) < n_events:
                failed.append(f.split('/')[-1] +
                              ' -> {0:.1f}% of desired\n'.format(len(d) * 100 / n_events))
            else:
                n_passed += 0
        self.n_passed = n_passed
        self.len = len(dfs)
        self.failures = failed
        return _generate_summary(self)

    def run_suite(self, save_output=True, output_dir='./validation/',
                  fname='flow_cytometry_validation.txt'):
        """
        Executes the entire test suite and saves the output to disk.

        Parameters
        ----------
        save_output : bool
            If True, the output will be saved as a `.txt` file.
        output_dir : str
            Path to the output directory. Default is in the current working
            directory in a folder called 'validation'. If this directory
            does not exist, it will be created.
        fname : str
            Name of the saved file. Default name is
            `flow_cytometry_validation.txt`

        Returns
        -------
        output : str
            Compiled string of the data validation results.
        """
        fns = [self.test_flow_cols(), self.test_gate(), self.test_positivity(),
               self.test_event_count()]
        output = """
*****************************************************************************
                Test suite for flow cytometry measurements
*****************************************************************************
"""
        for f in fns:
            output += f

        if save_output:
            if os.path.isdir('./validation') == False:
                os.mkdir('./validation')
            with open('{}{}'.format(output_dir, fname), 'w') as f:
                f.write(output)
        return output


class ImageValidation():
    def __init__(self):
        return True

    def test_fname():
        return True

    def test_dtype():
        return True

    def test_duration():
        return True

    def run_suite():
        return True


class DilutionDataValidator():
    """
    Data validation suite for the various data sets involved in the dilution
    experiment.

    """
    def __init__(self, dirs=None, clist_file=None, fluct_file=None):
        return True

    def test_positivity():
        return True

    def test_complete_clist():
        return True

    def test_error_fraction():
        return True

    def test_num_divisions():
        return True

    def run_suite():
        return True


class FoldChangeValidation():
    def __init__(self):
        return True

    def check_bounds():
        return True

    def run_suite():
        return True
