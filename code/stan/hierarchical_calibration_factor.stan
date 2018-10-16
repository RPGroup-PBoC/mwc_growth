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
    * @param alpha: Fluorescenc calibration factor in units of a.u. / molecule
    * @param N: Total number of measurements 
    **/
    real GammaApproxBinom_lpdf(real I1, real I2, real alpha, real p) { 
            return -log(alpha) + lgamma(((I1 + I2) / alpha) + 1) - lgamma((I1 / alpha) + 1)
                        - lgamma((I2 / alpha) + 1) + (I1 / alpha) * log(p) + (I2 / alpha) * log(1-p)
        }
    }

data {
    //Dimensional parameters
    int<lower=1> J_media; // Number of unique growth media
    int<lower=1> J_run; // Number of unique experimental across entire data set
    int<lower=1> N; // total number of measurements for fluctuations
    int<lower=1, upper=J_media> media_idx[N];
    int<lower=1, upper=J_run>  run_idx[N];
    
    // Experimental parameters
    real<lower=0> I_1[N]; // Observed mean pixel intensity of daughter cell 1
    real<lower=0> A_1[N]; // Area of cell 1 in square pixels
    real<lower=0> I_2[N]; // Observed mean pixel intensity of daughter cell 2 
    real<lower=0> A_2[N]; // Area of cell 2 in square pixels
    real<lower=0, upper=1> frac_area[N];
    
}
   
parameters {
    // Hyper parameters
    real<lower=0, upper=2^12>  alpha_mu[J_media]; // Hyperparameter for alpha

    // Low-level parameters
    real<lower=0, upper=2^12>  alpha_run[J_run]; // Low-level parameter for experimental alpha 
    real<lower=0> sigma[J_media]; // Hyperparameter for variance

    // Single-cell parameters
    real<lower=0> I1_tot[N];
    real<lower=0> I1_sigma[N];
    real<lower=0> I2_tot[N];
    real<lower=0> I2_sigma[N];
    real<lower=0, upper=1> prob[N];
    real<lower=0, upper=1> sigma_prob[N];
}

model {
    // Define the hyperpriors. 
    alpha_mu ~ lognormal(0, 5);
    sigma ~ normal(0, 100);

    // Define low-level priors
    I1_sigma ~ lognormal(0, 3);
    I2_sigma ~ lognormal(0, 3);
    
    for (i in 1:J_media) {
        for (j in 1:J_run) {
            alpha_run[j] ~ normal(alpha_mu[i], sigma[i]);
        }
    }
    
    // Iterate through each measurement and compute the likelihood
    for (i in 1:N) {
        // Define prior for Intensity. 
        I1_tot[i] ~ normal(I_1[i] * A_1[i], I1_sigma[i]);
        I2_tot[i] ~ normal(I_2[i] * A_2[i], I2_sigma[i]);

        // Define prior for partitioning probability
        prob[i] ~ normal(frac_area[i], sigma_prob[i]);

        // Evaluate likelihood.
        I1_tot[i] ~ GammaApproxBinom(I2_tot[i],alpha_run[run_idx[i]], prob[i]);
    } 
}
