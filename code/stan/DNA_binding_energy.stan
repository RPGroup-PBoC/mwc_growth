/* ****************************************************************************** 
*  Estimation of DNA binding Energy
* ------------------------------------------------------------------------------ 
* Author: Griffin Chure
* Last Modified: September 27, 2019
* Licewnse: MIT
* 
* Description
* ------------------------------------------------------------------------------ 
* This model infers the DNA binding energy from a collection of measurements of 
* the fold-change in gene expression.  Note that this model *does not* sample the
* allosteric energy difference and it is assumed that 100% of the repressors are 
* active at all points in time. 
* *****************************************************************************/ 
data { 
    // Dimensional parameters
    int<lower=1> J; // Number of unique conditions
    int<lower=1> N; // Total number of measurements
    int<lower=1, upper=J> idx[N]; // ID vector for condition

    // Input parameters
    vector<lower=1>[N] repressors; // Number of repressors per cell
    real<lower=1> Nns; // Number of nonspecific binding sites

    // Observed quantities
    vector[N] foldchange;
}

parameters {
    vector[J] epRA; // DNA binding energy
    vector<lower=0>[J] sigma; // Homoscedastic error
}

model {
    // Compute the mean fold-change in gene expression. 
    vector[N] mu = 1 ./ (1 + (repressors ./ Nns) .* exp(-epRA[idx]));

    // Define the priors. 
    epRA ~ normal(12, 6); // Same prior as used in Chure et al 2019 PNAS 
    sigma ~ normal(0, 0.1); // Sampe prior as used in Chure et al 2019 PNAS

    // Evaluate the likelihood. 
    foldchange ~ normal(mu, sigma[idx]);
    }