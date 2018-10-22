/* 
* Calibration Factor Estimation
* ---------------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description 
* ---------------------------------------------------
* This model samples the posterior probability distribution
* for the calibration factor between observed fluorescence 
* and fluorophore copy number. This model assumes that 
* the noise in measurement is negligible compared to the 
* observed fluorescence, making individual fluorescence 
* measurements behave as delta functions. 
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
    real GammaApproxBinom_lpdf(vector I1, vector I2, real alpha) {
        return sum(-log(alpha) + lgamma(((I1 + I2) ./ alpha) + 1) - lgamma((I1 ./ alpha) + 1) - lgamma((I2 ./ alpha) + 1) - ((I1 + I2) ./ alpha) * log(2));
    } 
}
     
data {
    int<lower=0> N; // Number of data points
    vector<lower=0>[N] I1; // Observed fluorescence of daughter cell 1
    vector<lower=0>[N] I2; // Observed fluorescence of daughter cell 2
}

parameters {
    // Generate non-centered modifiers
    real<lower=0> alpha_mu;
    real alpha_raw; 
    real<lower=0> tau;
}

transformed parameters{
    real alpha = alpha_mu + tau * alpha_raw;
}


model {    
    alpha_raw ~ normal(0, 2);
    alpha_mu ~ lognormal(2, 2);
    tau ~ normal(0, 1);
    I1 ~ GammaApproxBinom(I2, alpha);  
}
