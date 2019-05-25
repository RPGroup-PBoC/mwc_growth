/* ***************************************************************************
* DNA binding energy inference from single-cell data
* -----------------------------------------------------------------------------
*  Author: Griffin Chure
*  License: MIT
*   
* Description
* -----------------------------------------------------------------------------
* This stan model takes a collection og single-cell measurements of fold-change 
* values and estimates the DNA binding energy assuming the mean follows the
* prescribed theory. This model assumes that the repressor copy number is a
* delta function at the given value. This is an assumpting that should be relaxed
* in a more proper case. Thus, this model is more for exploratory purposes.
*/

data {
    int<lower=1> N; // Total number of data points
    vector[N] fc; // Set of observed fold-change measurements
    vector[N] R; // Median repressor copy number count 
}

parameters {
    real ep_RA;
    real<lower=0> sigma;
}

model {
    vector[N] mu;
    ep_RA ~ normal(-12, 6);
    sigma ~ normal(0, 0.1);
    mu = 1 ./ (1 + (1 / (1 + exp(-4.5))) * (R ./ 4600000) * exp(-ep_RA));
    fc ~ normal(mu, sigma);
}


