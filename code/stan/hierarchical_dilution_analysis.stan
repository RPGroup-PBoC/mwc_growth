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
    int<lower=1> N_meas; // Total number of repressor copy number measurements
    int<lower=1, upper=J_conc> index_2[N_meas]; // Indices for unique concentrations
        
    // Experimental data for calibration factor determination
    vector<lower=0>[N_fluct] I_1; // Observed mean pixel intensity of daughter cell 1
    vector<lower=0>[N_fluct] I_2; // Observed mean pixel intensity of daughter cell 2 

    // Experimental data for repressor copy number calculation
    vector<lower=0>[N_meas] mcherry; // Integrated mCherry expression for each cell
    vector<lower=0>[N_meas] yfp; // Integrated mCherry expression for each cell
}
  
transformed data {
    vector[N_meas] log_mch_data = log(mcherry);
    vector[N_meas] log_yfp_data = log(yfp);
}
   
   
parameters {
    // Define how hyperparameter varies
    real<lower=0> tau_alpha; 
    
    // Top-level parameters
    real log_alpha_1; // Calibration factor for particular growth medium
    vector[J_conc] log_mch_1; // mCherry expression for each concentration
    vector[J_conc] log_yfp_1; // mCherry expression for each concentration
    vector[J_conc] log_mch_sigma_1; //Homoscedastic error in mCherry expression for each concentration
    vector[J_conc] log_yfp_sigma_1; //Homoscedastic error in mCherry expression for each concentration
    
    // Level-1 Parameters
    vector[J_exp] log_alpha_2_raw; // Non-centered parameterization for cal factor 
}

transformed parameters {
    // Non-centered parameterization for means
    vector[J_exp] log_alpha_2 = log_alpha_1 + tau_alpha * log_alpha_2_raw;  
    real alpha_1 = exp(log_alpha_1);
    vector[J_exp] alpha_2 = exp(log_alpha_2);
    vector[J_conc] mch_sigma_1 = exp(log_mch_sigma_1);
    vector[J_conc] yfp_sigma_1 = exp(log_yfp_sigma_1);   
}
  
model {
    // Define the hyperpriors.
    log_alpha_1 ~ normal(0, 3);
    log_mch_1 ~ normal(0, 3);
    log_yfp_1 ~ normal(0, 3);
    log_mch_sigma_1 ~ normal(0, 1);
    log_yfp_sigma_1 ~ normal(0, 1);
    
    // Priors on hyperparameter variation
    tau_alpha ~ normal(2, 2);
  
    // Define priors on non-centering
    log_alpha_2_raw ~ normal(0, 1);

    //  Likelihood for calibration factor
    I_1 ~ GammaApproxBinom(I_2, alpha_2[index_1]);     
    
    // Likelihood for mCherry measurements. 
    log_mch_data ~ normal(log_mch_1[index_2], mch_sigma_1[index_2]);
    log_yfp_data ~ normal(log_yfp_1[index_2], yfp_sigma_1[index_2]);
}

generated quantities {
   vector[J_conc] mch_1 = exp(log_mch_1);
   vector[J_conc] yfp_1 = exp(log_yfp_1);
   vector[J_conc - 1] avg_rep = mch_1[2:] ./ alpha_1;
   vector[J_conc - 1] fc = yfp_1[2:] ./ yfp_1[1];
}
