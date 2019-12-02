/* *****************************************************************************
* Pooled Estimation of Temperature Dependent Entropic Change
* ------------------------------------------------------------------------------- 
* Author: Griffin Chure
* Last Modified: October 19, 2019
* License: MIT
*
* Description
* ------------------------------------------------------------------------------- 
* This model infers a temperature dependent entropic term for the DNA binding
* energy and the allosteric energy, given measurements of the fold-change in 
* gene expression. It does so using all supplied data.
* *****************************************************************************/
data {
    // Define dimensional parameters
    int<lower=1> J; // Number of unique conditions
    int<lower=1> N; // Total number of fold-change measurements
    int<lower=1, upper=J> idx[N]; // Identification vector for each point. 

    // Define input constants. 
    real ref_temp; // Reference temperature in K
    real ref_epRA; // Reference state DNA binding energy
    real<lower=1> Nns; // Number of nonspecific binding sites. 
    real ref_epAI; // Reference allosteric energy difference. 
    vector[J] temp; // Temperature vector in K 
    vector<lower=0>[N] repressors; // Repressor count per replicate

    // Define the observed quantity
    vector[N] foldchange; // Observed fold-change in gene expression
}

transformed data { 
    vector[N] log_fc = log(foldchange);
}

parameters { 
    real delta_SR;
    real delta_SAI;
    real<lower=0> sigma; // Homoscedastic error
}

transformed parameters { 
    // Compute the modified binding energies;
    vector[J] epRA_star = delta_SR * (ref_temp - temp) + ref_epRA; 
    vector[J] epAI_star = delta_SAI * (ref_temp - temp) + ref_epAI;
}

model { 

    //  probability of a repressor being active
     vector[N] pact = 1 ./ (1 + exp(-epAI_star[idx]));

     // Compute the mean fold-change in gene expression
     vector[N] mu = -1 * log(1 + pact .* (repressors ./ Nns) .* exp(-epRA_star[idx])); 

     // Define the priors
     delta_SR ~ normal(0, 0.1);
     delta_SR ~ normal(0, 0.1);
     sigma ~ normal(0, 0.1);

     // Evaluate the likelihood
     log_fc ~ normal(mu, sigma);

 }
