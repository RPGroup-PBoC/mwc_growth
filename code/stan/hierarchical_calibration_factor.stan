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
    int<lower=1> J_exp; // Number of unique experiments across entire data set
    int<lower=1> N_fluct; // Total number of measurements for fluctuations
    int<lower=1, upper=J_exp> index_1[N_fluct]; // Indices for fluctuation measurements

    // Dimensional parameters for repressor copy number calculation
    int<lower=1> J_conc; // Number of unique ATC concentrations    
    int<lower=1> J_conc_exp;
    int<lower=1> N_mch; // Total number of repressor copy number measurements
    int<lower=1, upper=J_conc> index_2[J_conc_exp]; // Indices for unique concentrations
    int<lower=1, upper=J_conc_exp> index_3[N_mch]; // Indices for individual measurements
        
    // Experimental data for calibration factor determination
    vector<lower=0>[N_fluct] I_1; // Observed mean pixel intensity of daughter cell 1
    vector<lower=0>[N_fluct] I_2; // Observed mean pixel intensity of daughter cell 2 

    // Experimental data for repressor copy number calculation
    vector<lower=0>[N_mch] mcherry; // Integrated mCherry expression for each cell
}
   
parameters {
    // Define how hyperparameters vary
    real<lower=0> tau_alpha; 
    real<lower=0> tau_mCherry;
    real<lower=0> tau_sigma;
    
    // Top-level parameters
    real<lower=0> log_alpha_1; // Calibration factor for particular growth medium
    vector[J_conc] log_mCherry_1; // mCherry expression for each concentration
    vector[J_conc] log_sigma_1; //Homoscedastic error in mCherry expression for each concentration
    
    // Level-1 Parameters
    vector[J_exp] log_alpha_2_raw; // Non-centered parameterization for cal factor
    
    // Level-2 parameters 
    vector[J_conc_exp] log_mCherry_2_raw;
    vector[J_conc_exp] log_sigma_2_raw;    
    
}

transformed parameters {
    // Non-centered parameterization for means
    vector<lower=0>[J_exp] log_alpha_2 = log_alpha_1 + tau_alpha * log_alpha_2_raw; 
    vector[J_conc_exp] log_mCherry_2 = log_mCherry_1[index_2] + tau_mCherry * log_mCherry_2_raw;
    vector[J_conc_exp] log_sigma_2 = log_mCherry_2[index_2] + tau_sigma * log_sigma_2_raw;
    
    // Bring parameters back to linear scale. 
    vector[J_conc_exp] sigma_2 = exp(log_sigma_2);
    vector[J_conc_exp] mCherry_2 = exp(log_mCherry_2);
    vector[J_conc] mCherry_1 = exp(log_mCherry_1);
    real alpha_1 = exp(log_alpha_1);
    vector[J_exp] alpha_2 = exp(log_alpha_2);
  }
  
model {
    // Define the hyperpriors.
    log_alpha_1 ~ normal(2, 2);
    log_mCherry_1 ~ normal(4, 2);
    log_sigma_1 ~ normal(1, 2);
    
    // Priors on hyperparameter variation
    tau_alpha ~ normal(2, 2);
    tau_mCherry ~ normal(0, 10);
    tau_sigma ~ normal(0, 10);
    
    // Define priors on non-centering
    log_alpha_2_raw ~ normal(0, 1);
    log_mCherry_2_raw ~ normal(0, 1);
    log_sigma_2_raw ~ normal(0, 1);
    

    //  Likelihood for calibration factor
    I_1 ~ GammaApproxBinom(I_2, alpha_2[index_1]);     
    
    // Likelihood for mCherry measurements. 
    mcherry ~ normal(mCherry_2[index_3], sigma_2[index_3]);
}

//generated quantities {
 //   vector[J_conc] avg_rep = mCherry_1 ./ alpha_1;
// }
