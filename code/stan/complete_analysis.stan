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
    // Dimensional Parameters
    int<lower=1> N_dil; // Total number of dilution measurements
    int<lower=1> N_con; // Total number of control measurements
    int<lower=1> J; // Total number of replicates
    int<lower=1, upper=J> idx_dil[N_dil]; // Identification vector for dilution measurements
    int<lower=1, upper=J> idx_con[N_con]; // Identification vector for control measumrements
    
    // Observed parameters
    vector<lower=0>[N_dil] I1;
    vector<lower=0>[N_dil] I2;
    vector<lower=0>[N_con] auto_fluo;
}

parameters {
    // Hyper parameters
    real<lower=0> alpha_mu;
    real<lower=0> alpha_sigma;
    real<lower=0> auto_mu_hyper;
    real<lower=0> auto_sigma_hyper;
    
    // Low level parameters -- Dilution Measurements
    vector<lower=0>[J] alpha;
    
    // Low level parameters -- Control Measurements. 
    vector<lower=0>[J] auto_mu;
    vector<lower=0>[J] sigma;
}

model {
    // Define the hyper priors
    alpha_mu ~ lognormal(3, 3);
    alpha_sigma ~ normal(0, 100);
    auto_mu_hyper ~ lognormal(3, 3);
    auto_sigma_hyper ~ normal(0, 100);
    
    // Low Level priors -- Dilution Measurements
    alpha ~ normal(alpha_mu, alpha_sigma);
 
    // Low level priors -- Control Measurements
    auto_mu ~ normal(auto_mu_hyper, auto_sigma_hyper);
    
    // Evaluate the likelihoods
    auto_fluo ~ normal(auto_mu[idx_con], sigma[idx_con]);
    I1 - auto_mu_hyper ~ GammaApproxBinom(I2 - auto_mu_hyper, alpha[idx_dil]);
}
