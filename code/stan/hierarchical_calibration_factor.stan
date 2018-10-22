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
    //Dimensional parameters
    int<lower=1> J_1; // Number of unique growth media
    int<lower=1> J_2; // Number of unique experiments across entire data set
    int<lower=1> N; // total number of measurements for fluctuations
    int<lower=1, upper=J_1> index_1[J_2];
    int<lower=1, upper=J_2> index_2[N];
    
    // Experimental parameters
    vector<lower=0>[N] I_1; // Observed mean pixel intensity of daughter cell 1
    vector<lower=0>[N] I_2; // Observed mean pixel intensity of daughter cell 2 
}
   
parameters {
    // Define how hyperparameters vary
    real<lower=0> tau_alpha;
    
    // Level-1 parameters 
    vector<lower=0>[J_1] alpha_1; 
    
    // Level-2 parameters
    vector[J_2] alpha_2_tilde;
}

transformed parameters {
    vector<lower=0>[J_2] alpha_2 = alpha_1[index_1] + tau_alpha * alpha_2_tilde;
  }

model {
    // Define the hyperpriors.
    alpha_1 ~ lognormal(2, 2);
    tau_alpha ~ normal(0, 1);
    alpha_2_tilde ~ normal(0, 10);

    // Iterate through each measurement and compute the likelihood 
    I_1 ~ GammaApproxBinom(I_2, alpha_2[index_2]);     
}
