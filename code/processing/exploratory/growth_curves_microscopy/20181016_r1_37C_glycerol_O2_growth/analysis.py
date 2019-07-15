#-*- coding: utf-8 -*-
import sys 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pystan 
sys.path.insert(0, '../../../../')
import mwc.viz
import mwc.stats
colors = mwc.viz.personal_style()
import altair as alt

# Experimental constants
DATE = 20181016
RUN_NO = 1
TEMP = 37
CARBON = 'glycerol'
OPERATOR = 'O2'

# Load the data. 
data = pd.read_csv(f'output/{DATE}_r{RUN_NO}_{TEMP}C_{CARBON}_{OPERATOR}_growth.csv')
data

alt.Chart(data).mark_point().encode(
x='time_min:Q',
y='area:Q',
color='colony_idx:O')