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
functions{
    /** 
    * Approximate the Binomial distirubution for continuous variables 
    * as a ratio of Gamma functions 
    * 
    * @param I1: Observed fluorescence of daughter cell 1. 
    * @param I2: Observed fluorescence of daughter cell 2.
    * @param alpha: Fluorescence calibration factor in units of a.u. / molecule
    * @param N: Total number of measurements 
    **/
    real GammaApproxBinom_lpdf(vector I1, vector I2, vector alpha ) { 
            return sum(-log(alpha))  + sum(lgamma(((I1 + I2) ./ alpha) + 1) - lgamma((I1 ./ alpha) + 1)
                        - lgamma((I2 ./ alpha) + 1) - ((I1 + I2) ./ alpha) * log(2));
        }
        
    /**
    * Calculate the fold-change in gene expression.
    *
    * @param R The number of repressors per cell
    * @param Nns The number of nonspecific repressor binding sites.
    * @param ep_r The  binding energy of the repressor to the DNA in kBT.
    * @param ep_ai The energetic difference between the active and inactive
    *        states of the repressor in kBT.
    **/
    vector fold_change(vector R, real Nns, real ep_r, real ep_ai) {
        // Compute the various componenets piecewise for simplicity.
        real pact = 1 / (1 + exp(-ep_ai));
        return 1 ./ (1 + pact * (R / Nns) * exp(-ep_r));
        } 
}

data {
    // Define dimensional parameters.
    int<lower=1> J_conc; //Number of unique ATC concentrations.
    int<lower=1> J_runs; // Number of unique experimental runs.
    int<lower=1> J_conc_runs; // Number of unique experimental runs per ATC concentration. 
    int<lower=1> N_fluct; // Number of fluctuation measurements.
    int<lower=1> N_fc; // Number of fold-change measurements.
    
    // Define identification vectors. 
    int<lower=1, upper=J_runs> fluct_index[N_fluct];
    int<lower=1, upper=J_runs> fc_index_run[N_fc];
    int<lower=1, upper=J_conc_runs> fc_index_replicates[N_fc];
    int<lower=1, upper=J_conc> fc_index_conc[J_conc_runs];
    
    
    // Define measured parameters
    vector<lower=0>[N_fluct] I_1;
    vector<lower=0>[N_fluct] I_2;
    vector<lower=0>[N_fc] mCherry;
    vector[N_fc] fc; // Measured foldchange
}

parameters {
    // Define hyper parameters
    real ep_RA;
    real sigma;
    vector<lower=0>[J_conc] R_mu;
    vector<lower=0>[J_conc] R_sigma;
    vector<lower=0>[J_conc] R_sigma_mu;
    vector<lower=0>[J_conc] R_sigma_sigma;
    real fc_mu[J_conc];
    real<lower=0> fc_sigma[J_conc];
    real<lower=0> fc_sigma_mu[J_conc];
    real<lower=0> fc_sigma_sigma[J_conc];
     
    // Define low-level parameters
    vector<lower=0>[J_runs] alpha;
    vector<lower=0>[J_conc_runs] R_mu_run;
    vector<lower=0>[J_conc_runs] R_sigma_run;
    real fc_mu_run[J_conc_runs];
    real<lower=0> fc_sigma_run[J_conc_runs];
}

model {
    // Define the top-level priors
    alpha ~ lognormal(0, 3);
    ep_RA ~ normal(0, 10);
    sigma ~ normal(0, 1);
    R_mu ~ lognormal(0, 3);
    R_sigma ~ lognormal(0, 2);
    R_sigma_mu ~ normal(0, 10);
    R_sigma_sigma ~ normal(0, 1);
    fc_mu ~ normal(0, 1);
    fc_sigma ~ normal(0, 1);
    fc_sigma_mu ~ normal(0, 1);
    fc_sigma_sigma ~ normal(0, 1);
    
    
    // Define the low-level priors
    R_mu_run ~ normal(R_mu[fc_index_conc], R_sigma[fc_index_conc]);
    R_sigma_run ~ normal(R_sigma_mu[fc_index_conc] ,R_sigma_sigma[fc_index_conc]);
    fc_mu_run ~ normal(fc_mu[fc_index_conc], fc_sigma[fc_index_conc]);   
    fc_sigma_run ~ normal(fc_sigma_mu[fc_index_conc], fc_sigma_sigma[fc_index_conc]);   

    
    // Compute the calibration factor likelihood
    I_1 ~ GammaApproxBinom(I_2, alpha[fluct_index]);
    
    // Compute the mean repressor copy numbers per experimental run   
    mCherry ~ normal(R_mu_run[fc_index_replicates] .* alpha[fc_index_run], R_sigma_run[fc_index_replicates] .* alpha[fc_index_run]);
    
    // Compute the mean fold-change per concentration
    fc ~ normal(fc_mu_run[fc_index_replicates], fc_sigma_run[fc_index_replicates]);
        
    // Compute the likelihood for ep_RA
    fc_mu ~ normal(fold_change(R_mu, 4.6E6, ep_RA, 4.5), sigma);    
}