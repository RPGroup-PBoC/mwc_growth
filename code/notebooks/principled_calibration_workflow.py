#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'code/notebooks'))
	print(os.getcwd())
except:
	pass
#%% [markdown]
# # A Principled Bayesian Inference of a Fluorescence Calibration Factor
# 
# © 2019 Griffin Chure. This work is licensed under a [Creative Commons Attribution License CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). All code contained herein is licensed under an [MIT license](https://opensource.org/licenses/MIT).
# 
# --- 

#%%
import numpy as np
import pandas as pd
import pystan 
import mwc.viz
import mwc.stats
import mwc.bayes
import scipy.stats
import matplotlib.pyplot as plt
import bokeh.io
import bokeh.plotting
import bokeh.palettes
import multiprocessing as mp
import joblib
import mwc.viz
import tqdm
seed = np.random.seed(666)
colors, colors_list = mwc.viz.bokeh_theme() 
bokeh.io.output_notebook()
get_ipython().run_line_magic('load_ext', 'stanmagic')

#%% [markdown]
# In this notebook, we lay out a principled workflow for the parameter estimation of a fluorescence calibration factor. This notebook covers the inference of the calibratoin factor from a single measurement. A principled analysis of a hierarchical approach will be written at a later date. 
#%% [markdown]
# ## Writing A Physical Model
# 
# It has become trivial in the era of molecular biology to label your favorite protein with a reporter that fluoresces at your favorite wavelength. Similarly, single-cell and single-molecule microscopy has become commonplace, making precise measurement of total cell fluorescence and/or localization a relatively painless procedure. However, it is still remarkably difficult to translate that precise measurement of fluorescence into the absolute copy number of your protein. 
# 
# By measuring the fluctuations in intensity between a mother/daughter pair after cell division, we can back calculate how bright a single molecule of interest is in arbitrary unitls, permitting the relatively easy calculation of protein copy number. We operate under the assumption that protein degradation is negligible and that protein production is ceased immediately before division. Using these assumptions, we can say that the total intensity of a given cell $I_\text{tot}$ is related to the number of fluorescent proteins and their relative brightness,
# 
# $$
# I_\text{tot} = \alpha N_\text{tot}, \tag{1}
# $$
# 
# where $N_\text{tot}$ is the total number of proteins and $\alpha$ is the brightness of a single fluorophore. Assuming there is no more production or deradation, we can relate the intensity of the mother cell to the daughter cells by dictating that the fluorescence must be conserved,
# 
# $$
# I_\text{tot} = I_1 + I_2 = \alpha (N_1 + N_2), \tag{2}
# $$
# 
# where we have used $I_1$ and $I_2$ to represent the total intensity of daughter cells 1 and 2. By looking at the fluctuations in intensity between the two daughter cells, followed by invocation of the mean and variance of the Binomial distribution, we arrive at the simple relation that
# 
# $$
# \langle (I_1 - I_2)^2 \rangle = \alpha I_\text{tot}. \tag{3}
# $$
# 
# Of course, lumped in with $\alpha$ is all of the features of the detector, fluorophore quantum efficiency, and other minutae of measurement. While incorporating these details into the generative model building is the proper thing to do, I know from experience in these experiments and quantitative biological. microscopy in general that the noise in these measurements is much smaller than the noise of the biological system. To this end, we will neglect them for simplicity. 
#%% [markdown]
# ## Building a Generative Statistical Model
#%% [markdown]
# The posterior probability distribution for our parameter of interest $\alpha$ is given by Bayes' theorem as
# 
# $$
# g(\alpha\,\vert\,[I_1, I_2]) \propto f([I_1, I_2]\,\vert\, \alpha)g(\alpha)
# \tag{4},
# $$
# 
# where I have used $g$ and $f$ to denote probability density functions over parameters and data, respectively. As $I_1$ and $I_2$ are related to each other (through Eq. 2), the likelihood $f([I_1, I_2]\,\vert\, \alpha)$ can be rewritten in the form,
# 
# $$
# f([I_1, I_2]\,\vert\,\alpha) = f([I_1]\,\vert\, \alpha, [I_2]).\tag{5}
# $$
# 
# Ultimately, we would like to translate the observed intensity into the nubmer of fluorescent proteins in the cell. Using Eq. 1 and the change-of-variables formula, Eq. 5 can be written in terms of $N_1$ as
# 
# $$
# f([I_1]\,\vert\, \alpha, [I_2]) = f([N_1]\,\vert\,\alpha, [I_2]) {d N_1 \over d I_1} = {1 \over \alpha^k}f([N_1]\,\vert\,\alpha,[I_2]),\tag{6}
# $$
# 
# where $k$ is the number of observed division events in the set of observed intensities $[I_1, I_2]$. 
# 
# The number of proteins present in the arbitrarily chosen daughter cell 1 ($N_1$) should be binomial distributed with $N_{tot}$ proteins in the mother cell. As such, the functional form of our likelihood can be written as
# 
# $$
# f([N_1]\,\vert\,\alpha,[I_2]) = {1 \over \alpha^k}\prod\limits_{i=1}^k {\Gamma\left({I_{1,i} + I_{2,i} \over \alpha} + 1\right) \over \Gamma\left({I_1 \over \alpha} + 1\right) \Gamma\left({I_2 \over \alpha} + 1\right)} 2^{-{I_1 + I_2 \over \alpha}},\tag{7}
# $$
# 
# where Gamma functions are used to avoid factorials of continuous variables. Note here we have used the simplifying assumption that the proteins are partitioned between the two proteins with an equal probability, $p = 1/2$. This is a fair approximation as the repressors to which the fluorescent proteins are covalently bound are always bound to the DNA and therefore are partitioned via segregation of the chromosomes. 
# 
# 
# With the likelihood in hand, we are now assigned the task of setting a prior for $\alpha$, which will be verified by prior predictive checks. This is certainly not an obvious endeavor as we do not have a good intuitoin for how bright a single fluorophore should be. However, we do know some limits. For example, it is not possible for the fluoropohore to have a negative value as the detector returns unsigned integeters per pixel as the readout. We also know that our resolution is diffraction limited meaning in the extreme case in which the repressors are frozen in a specific point, their fluorescence will be detected across several pixels. In our system, the interpixel distance of the camera is on the order of 50 nm / pixel meaning that a single fluorophore will have a diameter of 4 - 5 pixels (assuming the extreme scenario in which diffusion of the fluorophore is much slower than the frame rate of our camera). With a bit-depth of 12, the maximum value that a pixel can be assigned is $2^{12} - 1 = 4095$ arbitrary units. Thus, the maximum brightness of a single fluorophore that can be detected on our system is $ 4 \times 4095 \approx 1.5\times 10^4$ a.u. /molecule. This, of course is a large upper bound and in our system the fluorescence will be below this value. However, we can use these bounds to  define a prior. As flat priors are difficult to sample and a uniform prior is not appropriate, we can assign a log-normal prior centered at $\approx 5$ (~200 [a.u.]) with a scale parameter of 3.
# 
# ## Prior Predictive Checks
# The best way to check if this choice of prior is reasonable is to draw $N$ values from this prior distribution and use $\alpha$ to generate a data set of $M$ values for $I_1$ based on a set of $M$ protein copy numbers. To start, we will define a set protein copy numbers $N_1$ that cover physiologically reasonable values of intenisities.

#%%
# Define the set of copy numbers. 
M = 500 
N1 = np.random.gamma(10, 10, M).astype(int) # Gamma distributed to mimic bursty gene expression

#%% [markdown]
# We can now make $N$ draws out of the prior distribution for $\alpha$ and compute the observed intensities $I_1$. We'll make this as a data frame for use later. 

#%%
# Make N draws out of a lognormal distribution for alpha. 
N = 500
alpha_draws = np.random.lognormal(5, 2, N)

# Compute the intensities. 
n, a = np.meshgrid(N1, alpha_draws)
I1 = n * a

# Make the dataframe. 
ppc_df = pd.DataFrame(np.array([I1.flatten(), n.flatten(), a. flatten()]).T, 
                     columns=['I1', 'N1', 'alpha_true'])



#%%



#%%
# Compute percentiles of the data.
percentiles = [0.99, 0.95, 0.75, 0.5, 0.25, 0.1, 0.05, 0.01]
cmap = bokeh.palettes.Purples9
ax = bokeh.plotting.figure(width=600, height=400, 
                          x_axis_label='cell intensity [a.u.]',
                         y_axis_label='cumulative distribution',
                         x_axis_type='log')

percentile_df = pd.DataFrame([])
for g, d in ppc_df.groupby(['N1']):
    x = np.sort(d['I1'].values) 
    y = np.linspace(0, 1, len(d))
    ax.line(x, y, line_width=0.5, color=colors_list[0])

bokeh.io.show(ax)

#%% [markdown]
# These distributions look pretty good -- most of the generated intensities are between $10^2$ and $10^7$ a.u. per cell which is in line with what I would see in the extremes of my experiments. There s also barely any density below $1$ a.u. / cell, which is what we would want to see during the prior predictive checks. 
# 
# Satisfied with our choice of prior, we can now define the model using Stan, as given in the cell below. 

#%%
get_ipython().run_cell_magic('stan', '-v simple', 'functions{\n    /** \n    * Approximate the Binomial distirubution for continuous variables \n    * as a ratio of Gamma functions \n    * \n    * @param I1: Observed fluorescence of daughter cell 1. \n    * @param I2: Observed fluorescence of daughter cell 2.\n    * @param alpha: Fluorescenc calibration factor in units of a.u. / molecule\n    * @param N: Total number of measurements \n    **/\n    real GammaApproxBinom_lpdf(vector I1, vector I2, real alpha) {\n        return sum(-log(alpha) + lgamma(((I1 + I2) ./ alpha) + 1) - lgamma((I1 ./ alpha) + 1) - lgamma((I2 ./ alpha) + 1) - ((I1 + I2) ./ alpha) * log(2));\n    } \n\n}\n     \ndata {\n    int<lower=0> N; // Number of data points\n    vector<lower=0>[N] I1; // Observed fluorescence of daughter cell 1\n    vector<lower=0>[N] I2; // Observed fluorescence of daughter cell 2\n}\n\n\nparameters {\n    // Generate entered parameters\n    real<lower=0> alpha;\n}\n\n\nmodel {     \n    alpha ~ lognormal(5, 2);\n    I1 ~ GammaApproxBinom(I2, alpha);  \n}')


#%%
get_ipython().run_cell_magic('stan', '-v proper', 'data {\n    int<lower=1> N;\n    vector<lower=0>[N] I1;\n    vector<lower=0>[N] I2;\n}\n\nparameters {\n    vector<lower=0>[N] n_tot;\n    vector<lower=0>[N] n1;\n    real<lower=0> alpha;\n    real<lower=0> sigma;\n}\n\nmodel {\n    vector[N] n2;\n    n_tot ~ lognormal(3, 3);\n    sigma ~ normal(0, 1);\n    alpha ~ lognormal(3, 3);\n    n1 ~ normal(n_tot/2, n_tot/4);\n    n2 = n_tot - n1;\n    I1 ~ normal(alpha * n1, sigma);\n    I2 ~ normal(alpha * n2, sigma);\n}')


#%%
# Compile the model 
model = pystan.StanModel(model_code=simple.model_code)

#%% [markdown]
# ## Computational Tractability
#%% [markdown]
# Our next step in the workflow is to ensure that our model is well specified and computationally accessible. To do so, we perform [Simulation Based Calibration]() as put forward by Talts. et al. in 2018. To do so, given the draws from the prior, we generate $M$ simulated data sets for the partitioning with $N$ divisions each, where the intensity is dicted by the draw from the $\alpha$ prior distribution. Using this data set, we then perform inference on these data and compare the data-averaged posterior distribution for $\alpha$ to the true prior distribution. If everything works as expected, we should recover the prior distribution and have uniformly distributed rank statistics. In addition, we can check the degree to which the prior distribution shrinks (shrinkage) as well as the $z-$score, which quantifies how accurately the posterior captures the model configuration.
# 
# To permit direct comparison between inferential runs, we will use the same randomly generated protein copy numbers and only change the measured intensity $I$. 

#%%
# Define the number of divisions
M = 800
n_tot = np.random.gamma(2, 100, M).astype(int)
n1 = np.random.binomial(n_tot, p=0.5)
n2 = n_tot - n1




# %%
model = mwc.bayes.StanModel('../stan/calibration_factor.stan')

#%%
# Define a function to evaluate a single run. 
def sbc(alpha):
        df = pd.DataFrame([])
        # Compute intensities
        I1 = n1 * alpha
        I2 = n2 * alpha
    
        # Perform inference. 
        fit, samples = model.sample(dict(N=M, I1=I1, I2=I2), n_jobs=1) 

    
        # Compute the statistics. 
        samples['alpha'] = np.exp(samples['alpha'])
        mean_alpha = np.mean(samples['alpha'])
        shrinkage = 1 - (np.var(samples['alpha']) / np.var(alpha_draws))
        z_score = (mean_alpha - alpha) / np.std(samples['alpha'])
        df = df.append(dict(alpha_true=alpha, mean_alpha=mean_alpha, 
                           shrinkage=shrinkage, z_score=z_score), ignore_index=True)
        return df


#%%
result = joblib.Parallel(n_jobs=-1)(joblib.delayed(sbc)(a) for a in tqdm.tqdm(alpha_draws))


#%%
sbc_df = pd.concat(result)
# Compute the rank statistic. 
rank = [np.sum(sbc_df['mean_alpha'] < a) for a in alpha_draws]
sbc_df['rank'] = rank
sbc_df.to_csv('../../data/sbc_calibration_factor.csv', index=False)

# %%
sbc_df = pd.read_csv('../../data/sbc_calibration_factor.csv')
#%%
ax = bokeh.plotting.figure(width=600, height=400, x_axis_label='rank statistic', y_axis_label='cumulative distribution')
ax.step(np.sort(sbc_df['rank']), np.linspace(0, 1, len(sbc_df)), color=colors['dark_purple'])
bokeh.io.show(ax)


#%%
p = bokeh.plotting.figure(width=600, height=400, 
                          x_axis_label='shrinkage', y_axis_label='z-score',
                          x_range=[0, 1], y_range=[0, 5])
p.circle(sbc_df['shrinkage'], sbc_df['z_score'])
bokeh.io.show(p)


#%%
sbc_df


#%%

plt.hist(sbc_df['mean_alpha'], bins=100)
# plt.xscale('log')


#%%
output = joblib.Parallel(n_jobs=mp.cpu_count())(joblib.delayed(sbc)(a) for a in alpha_draws)


