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
            return sum(-N * log(alpha)  + sum(lgamma(((I1 + I2) ./ alpha) + 1) -
                     lgamma((I1 ./ alpha) + 1) - lgamma((I2 ./ alpha) + 1) - 
                     ((I1 + I2) ./ alpha) * log(2)));
        }
    
}
data {
    // Dimensional parameters for calibration factor determination
    int<lower=1> J_exp; // Number of unique experiments across entire data set
    int<lower=1> N_fluct; // Total number of measurements for fluctuations
    int<lower=1, upper=J_exp> index_1[N_fluct]; // Indices for fluctuation measurements
 
    // Experimental data for calibration factor determination
    vector<lower=0>[N_fluct] I_1; // Observed mean pixel intensity of daughter cell 1
    vector<lower=0>[N_fluct] I_2; // Observed mean pixel intensity of daughter cell 2 
}
   
parameters {
    // Define how hyperparameters vary
    real tau_alpha; 
    
    // Top-level parameters
    real<lower=0> log_alpha_1; // Calibration factor for particular growth medium 
    
    // Level-1 Parameters
    vector[J_exp] alpha_2_raw; // Non-centered parameterization for cal factor
}

transformed parameters {
    // Non-centered parameterization for means
    real alpha_1 = exp(log_alpha_1); 
    vector[J_exp] alpha_2 = alpha_1 + tau_alpha * alpha_2_raw; 
    }
  
model {
    // Define the hyperpriors.
    log_alpha_1 ~ gamma(2, 0.6); 
   
    // Priors on hyperparameter variation
    tau_alpha ~ normal(0, 1); 
    
    // Define priors on non-centering
    alpha_2_raw ~ normal(0, 10);    

    //  Likelihood for calibration factor
    I_1[index_1] ~ GammaApproxBinom(I_2[index_1], alpha_2[index_1], N_fluct);        
}

 