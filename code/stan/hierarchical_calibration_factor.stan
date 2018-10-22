/* 
* Hierarchical Model for Calibration Factor Inference
* ---------------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
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
    
}
data {
    // Dimensional parameters for calibration factor determination
    int<lower=1> J_1; // Number of unique growth media
    int<lower=1> J_2; // Number of unique experiments across entire data set
    int<lower=1> N; // Total number of measurements for fluctuations
    int<lower=1> N_auto;
    int<lower=1> N_delta;
    int<lower=1> M; // Total number of repressor copy number measurements
    int<lower=1, upper=J_1> index_1[J_2];
    int<lower=1, upper=J_2> index_2[N];
    int<lower=1, upper=J_1> auto_index_1[J_2];
    int<lower=1, upper=J_2> auto_index_2[N_auto];
    int<lower=1, upper=J_1> delta_index_1[J_2];
    int<lower=1, upper=J_2> delta_index_2[N_delta];
    int<lower=1, upper=J_1> fc_index_1[J_2];
    int<lower=1, upper=J_2> fc_index_2[M];
    
    // Experimental data for calibration factor determination
    vector<lower=0>[N] I_1; // Observed mean pixel intensity of daughter cell 1
    vector<lower=0>[N] I_2; // Observed mean pixel intensity of daughter cell 2 

    // Experimental data for fold-change and repressor copy number calculation
    vector<lower=0>[N_auto] autofluo;
    vector<lower=0>[N_delta] deltafluo;
    vector<lower=0>[M] mCherry_exp;
    vector<lower=0>[M] yfp_exp;
    vector<lower=1>[M] area; // Cell area in µm^2
}
   
parameters {
    // Define how hyperparameters vary
    real<lower=0> tau_alpha;
    real<lower=0> tau_auto;
    real<lower=0> tau_auto_sigma;
    real<lower=0> tau_delta;
    real<lower=0> tau_delta_sigma;
    
    // Level-1 parameters 
    vector<lower=0>[J_1] alpha_1; 
    vector<lower=0>[J_1] autofluo_1;
    vector<lower=0>[J_1] log_autofluo_sigma;
    vector<lower=0>[J_1] deltafluo_1;
    vector<lower=0>[J_1] log_deltafluo_sigma;
    
    // Level-2 parameters
    vector[J_2] alpha_2_raw;
    vector[J_2] autofluo_2_raw;
    vector[J_2] log_autofluo_sigma_2_raw;
    vector[J_2] deltafluo_2_raw;
    vector[J_2] log_deltafluo_sigma_2_raw;
}

transformed parameters {
    // Non-centered parameterization
    vector<lower=0>[J_2] alpha_2 = alpha_1[index_1] + tau_alpha * alpha_2_raw;
    vector<lower=0>[J_2] autofluo_2 = autofluo_1[auto_index_1] + tau_auto * autofluo_2_raw;
    vector<lower=0>[J_2] deltafluo_2 = deltafluo_1[auto_index_1] + tau_delta * deltafluo_2_raw;
    vector<lower=0>[J_2] log_autofluo_sigma_2 = log_autofluo_sigma[auto_index_1] + tau_auto_sigma *  log_autofluo_sigma_2_raw;
    vector<lower=0>[J_2] log_deltafluo_sigma_2 = log_deltafluo_sigma[delta_index_1] + tau_delta_sigma * log_deltafluo_sigma_2_raw;

    // Transform sigmas to linear space
    vector<lower=0>[J_2] autofluo_sigma = exp(log_autofluo_sigma);
    vector<lower=0>[J_2] autofluo_sigma_2 = exp(log_autofluo_sigma_2);
    vector<lower=0>[J_2] deltafluo_sigma_2 = exp(log_deltafluo_sigma_2);
  }

model {
    // Define the hyperpriors.
    alpha_1 ~ lognormal(2, 2);
    tau_alpha ~ normal(0, 1);
    alpha_2_tilde ~ normal(0, 10);

    //  Likelihood for calibration factor
    I_1 ~ GammaApproxBinom(I_2, alpha_2[index_2]);     

    // Likelihood for mean auto and mean delta fluorescence. 
    autofluo ~ normal(autofluo_2[auto_index_2], autofluo_sigma_2[auto_index_2]);
    deltafluo ~ normal(deltafluo_2[delta_index_2], deltafluo_sigma_2[delta_index_2]);
}

generated quantities {
    // Compute the number of repressors 
    vector[N] rep_per_cell = (mCherry_exp .* area) ./ alpha[fc_index_2];
    
    // Compute the fold-change in gene expression
    vector[N] foldchange = (yfp_exp - autofluo) ./ (deltafluo - autofluo);
}
