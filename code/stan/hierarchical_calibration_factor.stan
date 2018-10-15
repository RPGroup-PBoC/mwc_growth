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
    real GammaApproxBinom_lpdf(real I1, real I2, real alpha) { 
            return -log(alpha) + lgamma(((I1 + I2) / alpha) + 1) - lgamma((I1 / alpha) + 1)
                        - lgamma((I2 / alpha) + 1) - ((I1 + I2) / alpha) * log(2);
        }
    }

data {
    //Dimensional parameters
    int<lower=1> J_media; // Number of unique growth media
    int<lower=1> J_run; // Number of unique experimental across entire data set
    int<lower=1> N; // total number of measurements
    int<lower=1, upper=J_media> media_idx[N];
    int<lower=1, upper=J_run>  run_idx[N];
    
    // Experimental parameters
    real<lower=0> I_1[N]; // Observed intensity of daughter cell 1
    real<lower=0> I_2[N]; // Observed intensity of daughter cell 2 
}
   
parameters {
    real<lower=0, upper=2^12>  alpha_mu[J_media]; // Hyperparameter for alpha
    real<lower=0, upper=2^12>  alpha_run[J_run]; // Low-level parameter for experimental alpha 
    real<lower=0> sigma[J_media]; // Hyperparameter for variance
}

model {
    // Define the priors. 
    alpha_mu ~ lognormal(0, 5);
    sigma ~ normal(0, 100);
    
    
    for (i in 1:J_media) {
        for (j in 1:J_run) {
        alpha_run[j] ~ normal(alpha_mu[i], sigma[i]);
    }
    }
    
    // Iterate through each measurement and compute the likelihood
    for (i in 1:N) {
        I_1[i] ~ GammaApproxBinom(I_2[i],alpha_run[run_idx[i]]);
    } 
}