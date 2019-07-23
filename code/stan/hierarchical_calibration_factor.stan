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
    real GammaApproxBinom_lpdf(vector I1, vector I2, vector alpha, int N) { 
            return  sum(N * log(alpha)  + lgamma(((I1 + I2) ./ alpha) + 1) -
                     lgamma((I1 ./ alpha) + 1) - lgamma((I2 ./ alpha) + 1) - 
                     ((I1 + I2) ./ alpha) * log(2));
        }
    
}
data {
    // Dimensional parameters for calibration factor determination
    int<lower=1> J; // Number of unique days of experiments
    int<lower=1> N; // Total number of measurements for fluctuations
    int<lower=1, upper=J> idx[N];
 
    // Experimental data for calibration factor determination
    vector<lower=0>[N] I_1; // Observed mean pixel intensity of daughter cell 1
    vector<lower=0>[N] I_2; // Observed mean pixel intensity of daughter cell 2 
}
   
parameters {
    // Top-level parameters
    real<lower=1, upper=2^16> alpha_1; // Calibration factor for particular growth medium 
    real<lower=0> sigma;
    vector<lower=1>[J] alpha_2; // Non-centered parameterization for day cal factor

}
  
model {
    // Define the hyperpriors.
    alpha_1 ~ normal(0, 500); 
    sigma ~ normal(0, 10);
    
    // Define priors on non-centering
    alpha_2 ~ normal(alpha_1, sigma);    

    //  Likelihood for calibration factor
    I_1 ~ GammaApproxBinom(I_2, alpha_2[idx], N);        
}

 