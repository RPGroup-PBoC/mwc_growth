/* *****************************************************************************
* Estimation of Temperature Dependent Entropic Change
* ------------------------------------------------------------------------------- 
* Author: Griffin Chure
* Last Modified: September 27, 2019
* License: MIT
*
* Description
* ------------------------------------------------------------------------------- 
* This model infers a temperature dependent entropic term for the DNA binding
* energy and the allosteric energy, given measurements of the fold-change in 
* gene expression. This model can be run hierarchically if J > 1
* *****************************************************************************/

data {
    // Define dimensional parameters
    int<lower=1> J; // Number of unique conditions
    int<lower=1> N; // Total number of fold-change measurements
    int<lower=1, upper=J> idx[N]; // Identification vector for each point. 

    // Define input constants. 
    real ref_temp; // Reference temperature in C
    real ref_epRA; // Reference state DNA binding energy
    real<lower=1> Nns; // Number of nonspecific binding sites. 
    real ref_epAI; // Reference allosteric energy difference. 
    vector<lower=32, upper=42>[J] temp; // Temperature vector in C 
    vector<lower=1>[N] repressors; // Repressor count per replicate

    // Define the observed quantity
    vector[N] foldchange; // Observed fold-change in gene expression
}

transformed data {
    // Compute the relative temperature difference. 
    vector[J] rel_T = (temp + 273.15) ./ (ref_temp + 273.15); 
}

parameters { 
    real true_epRA;
    vector[J] delta_S_DNA; // Change in entropy for DNA binding
    // vector[J] delta_S_ALLO; // Change in entropy for allosteric energy difference
    vector<lower=0>[J]  sigma; // Homoscedastic error
}

transformed parameters { 
    // Compute the modified binding energies
    vector[J] epRA_star = rel_T .* (true_epRA - (temp + 273.15) .* delta_S_DNA);
    // vector[J] epAI_star = rel_T .* (true_epRA - (temp + 273.15) .* delta_S_ALLO);
}

model { 
     // probability of a repressor being active
    // vector[N] pact = 1 ./ (1 + exp(-epAI_star[idx]));

    // Compute the mean fold-change in gene expression
    vector[N] mu = 1 ./ (1 + (repressors ./ Nns) .* exp(-epRA_star[idx])); 

    // Define the priors
    true_epRA ~ normal(-12, 6);
    delta_S_DNA ~ normal(0, .1);
    // delta_S_ALLO ~ normal(0, .1);
    sigma ~ normal(0, 0.1);

    // Evaluate the likelihood
    foldchange ~ normal(mu, sigma[idx]);
}