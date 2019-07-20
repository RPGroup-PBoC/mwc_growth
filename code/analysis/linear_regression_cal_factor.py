#%%[markdown]
# # Estimation of $\alpha$ Via Bayesian linear Regression
# ---
#%%
import numpy as np 
import pandas as pd
import mwc.bayes
import mwc.viz
import mwc.stats
import bokeh.io
import bokeh.plotting

# Load the lineages
lineages = pd.read_csv('../../data/analyzed_lineages.csv')

# Choose a singlee carbon source.
gluc = lineages[(lineages['temp']==37) & (lineages['carbon']=='glucose')].copy()

#%%[markdown]
# I'm getting the feeling that my 'proper' Bayesian approach won't really work.
# There's too much uncertainty and it's proving to bee too sensitive to my
# choice of prior. Another way that I can go about this is by returning to the
# Rosenfeld approach and trying to use a binning scheme. Rather than just
# choosing an arbitrary bin and hoping for the best, I can try grouping by
# events ranging from 10 to the total number of measurements in, say, steps of 10,
# fitting $\alpha$ for each step. This is similar to how Brewster went about the
# problem, but without having to worry about bootstrapping to get a proper
# measurement error. 
#
# To start, I'm just going to take one replicate, set the bins, and perform the
# inference. Unfortunately, I'm worried this will be slow, but it should be
# parallelizeable. 
#

#%%
# Choose a date and replicate. We'll sort by count to get the most. 
date, run = gluc.groupby(['date', 'run_number'])['I_1'].count().idxmax()
samp = gluc[(gluc['date']==date) & (gluc['run_number']==run)].copy()

# Compute the summed values and the fluctuations. 
samp['summed'] = samp['I_1_sub'] + samp['I_2_sub']
samp['fluct'] = (samp['I_1_sub'] - samp['I_2_sub'])**2

# Determine how many binning groups I will have. 
step = 10
n_groups = int(len(samp) / step)

#
#%%[markdown]
# Okay, so for this particular replicate, I have around 50 inferences to perform
# (!). Hopefully it's not too slow. I think that assigning bins should be
# relatively easy to do 

#%%
def assign_bin_idx(df, n_events, sortby='summed'):
    _df = df.copy()
    # Sort the dataframe.
    _df.sort_values(sortby, inplace=True)

    # Set the bins.
    bin_idx = np.array([[i+1] * n_events for i in range(
                        int((len(_df) + n_events) / n_events))]).flatten()

    # Assign bins, trimming the left over
    _df['bin_idx'] = bin_idx[:len(_df)]

    # If the last bin doesn't have at least n_events, add them to the previous bin
    if (len(_df)%n_events) != 0:
        _df.loc[_df['bin_idx']==_df['bin_idx'].max(), 'bin_idx'] = _df['bin_idx'].max() - 1
    return _df
#%%[markdown]
# Tested and this works, although I absolutely should write a test suite for
# this.  Let's look at the fluctuations and means. 
#%%
#%%kI wrote a stan model for doing this called 'cal_factor_bin_simple.stan'.

#%%
import tqdm
model = mwc.bayes.StanModel('../stan/cal_factor_bin_simple.stan', force_compile=True) 
bin_size = 10 
bins = np.arange(0, len(samp)+bin_size, bin_size)
alpha_means = []
for b in tqdm.tqdm(bins):
    binned = assign_bin_idx(samp, 10)
    # Assemble the data dictionary. 
    data_dict = {'N':len(binned), 'J':binned['bin_idx'].max(),
            'idx':binned['bin_idx'],
            'summed':binned['summed'], 
            'fluct':binned['fluct']}

    fit, samples = model.sample(data_dict, iter=2000) #, control=dict(adapt_delta=0.95))
    alpha_means.append(np.mean(samples['alpha']))

#%%
p = bokeh.plotting.figure()
p.circle(bins, alpha_means)
p.line(bins, alpha_means)
bokeh.io.show(p)
#%%
p = bokeh.plotting.figure(x_axis_type='log', y_axis_type='log', width=500, height=300)

p.circle(binned['summed'], binned['fluct'], size=1, alpha=1, color=colors['black'])

summary = binned.groupby(['bin_idx']).agg(('mean', 'sem'))
summary['summed_low'] = summary['summed']['mean'] - summary['summed']['sem']
summary['summed_high'] = summary['summed']['mean'] + summary['summed']['sem']
summary['fluct_low'] = summary['fluct']['mean'] - summary['fluct']['sem']
summary['fluct_high'] = summary['fluct']['mean'] + summary['fluct']['sem']
p.segment(x0=summary['summed_low'], x1=summary['summed_high'], y0=summary['fluct']['mean'],
        y1=summary['fluct']['mean'], line_width=1, color=colors['purple'])
p.segment(x0=summary['summed']['mean'], x1=summary['summed']['mean'], y0=summary['fluct_low'],
        y1=summary['fluct_high'], line_width=1, color=colors['orange'])
p.circle(summary['summed']['mean'], summary['fluct']['mean'], line_color=colors['orange'],
                fill_color=colors['light_orange'], size=5, line_width=1)
I_tot = np.logspace(2, 5.5)
p.line(I_tot, np.mean(samples['alpha']) * I_tot, color='tomato')
bokeh.io.show(p) 
#


#%%
