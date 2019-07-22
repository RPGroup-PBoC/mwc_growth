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
    int<lower=1> J_day; // Number of unique days of experiments
    int<lower=1> K_rep; // Number of unique replicates for that day
    int<lower=1> N_fluct; // Total number of measurements for fluctuations
    int<lower=1, upper=J_day> day_idx[K_rep]; // Indices for the days 
    int<lower=1, upper=K_rep> rep_idx[N_fluct];
 
    // Experimental data for calibration factor determination
    vector<lower=0>[N_fluct] I_1; // Observed mean pixel intensity of daughter cell 1
    vector<lower=0>[N_fluct] I_2; // Observed mean pixel intensity of daughter cell 2 
}
   
parameters {
    // Define how hyperparameters vary
    real tau_alpha; 
    
    // Top-level parameters
    real<lower=1, upper=2^16> alpha_1; // Calibration factor for particular growth medium 
    
    // Level-1 Parameters
    vector[J_day] alpha_2_raw; // Non-centered parameterization for day cal factor
    vector[K_rep] alpha_3_raw; // Non-centered parameterization for replicate factor
}

transformed parameters {
    // Non-centered parameterization for means
    vector[J_day] alpha_2 = alpha_1 + tau_alpha * alpha_2_raw; 
    vector[K_rep] alpha_3 = alpha_2[day_idx] + tau_alpha * alpha_3_raw; 
    }
  
model {
    // Define the hyperpriors.
    alpha_1 ~ normal(0, 1000); 
   
    // Priors on hyperparameter variation
    tau_alpha ~ normal(0, 1); 
    
    // Define priors on non-centering
    alpha_2_raw ~ normal(0, 10);    
    alpha_3_raw ~ normal(0, 10);

    //  Likelihood for calibration factor
    I_1[rep_idx] ~ GammaApproxBinom(I_2[rep_idx], alpha_3[rep_idx], N_fluct);        
}

 