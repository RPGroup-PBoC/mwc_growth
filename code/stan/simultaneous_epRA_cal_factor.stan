/* 
* Fitting DNA Binding Energy
* -------------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description 
* ------------------------------------------------
* This model infers the calibration factor, repressor
* copy number, as well as the DNA binding energy 
* for a given data set. 
*/
#include functions.stan

data {
    // Define dimensional parameters.
    int<lower=1> J_conc; //Number of unique ATC concentrations.
    int<lower=1> J_runs; // Number of unique experimental runs.
    int<lower=1> N_fluct; // Number of fluctuation measurements.
    int<lower=1> N_fc; // Number of fold-change measurements (means).
    
    // Define identification vectors. 
    int<lower=1, upper=J_runs> fluct_index[N_fluct];
    int<lower=1, upper=J_runs> fc_index_run[N_fc];
    int<lower=1, upper=J_conc> fc_index_conc[N_fc];
    int<lower=1, upper=N_fc> fc_index_replicates[N_fc];
     
    // Define measured parameters.
    vector<lower=0>[N_fluct] I_1; // Intensity of daughter 1.
    vector<lower=0>[N_fluct] I_2; // Intensity of daughter 2.
    vector<lower=0>[N_fc] mean_mCherry; // Mean mCherry integrated cell intensity.
    vector[N_fc] fc; // Measured mean fold change per measurement. 
    
    // Define known parameters.
    real<lower=0> Nns; // Number of non-specific  binding sites.
    real ep_AI; // Energetic difference between allosteric states in kBT.
}

parameters {
    // Define hyper parameters.
    real ep_RA;
    real<lower=0> sigma;    
    vector<lower=0>[J_conc] R_mu;
    vector<lower=0>[J_conc] R_sigma;
    real<lower=0> alpha_mu;
    real<lower=0> alpha_sigma;
    
    // Define low-level parameters.
    vector[J_runs] alpha_run_raw;
}

transformed parameters {    
    
    // Uncenter the alpha parameter
    vector[J_runs] alpha_run = alpha_mu + alpha_run_raw * alpha_sigma;
    
    // Compute the number of repressors per cell based off of inferred alpha. 
    vector[N_fc] inferred_R = mean_mCherry ./ alpha_run[fc_index_run];
}

model {
    // Define the top-level priors.
    alpha_mu ~ normal(0, 10);
    alpha_sigma ~ normal(0, 1);
    alpha_run_raw ~ normal(0, 1);
    ep_RA ~ normal(0, 10);
    sigma ~ normal(0, 1);
    R_mu ~ lognormal(0, 3);
    R_sigma ~ lognormal(0, 2); 
    
    
    // Compute the calibration factor likelihood.
    I_1 ~ GammaApproxBinom(I_2, alpha_run[fluct_index]);
    
    // Compute the mean repressor copy numbers per experimental run.
    inferred_R ~ normal(R_mu[fc_index_conc], R_sigma[fc_index_conc]);
    
    // Compute the likelihood for ep_RA.
    fc ~ normal(fold_change(inferred_R, Nns, ep_RA, ep_AI), sigma);    
}