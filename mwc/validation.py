# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import os


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

    """Validation test suite for flow cytometry data sets.

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

    def test_temporal_stability():
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
