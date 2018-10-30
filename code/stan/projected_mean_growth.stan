/*
* Computed Mean Area Growth Curve Fitting
* ---------------------------------------
* Author: Griffin Chure
* License: MIT
* 
*/

data {
   int<lower=1> J; // Number of growth media
   int<lower=1> N; // Total number of measurements
   int<lower=1, upper=J> index_1[N];
   vector<lower=0>[N] area;
   vector<lower=0>[N] time;
}

transformed data {
    vector[N] log_area = log(area);
}

parameters {
    vector[J] log_A_0;
    vector[J] log_r;
    vector[J] log_sigma;
}

model {
    // Set the priors
    log_A_0 ~ normal(0, 1);
    log_r ~ normal(0, 1);
    log_sigma ~ normal(0, 2); 
    log_area ~ normal(log_A_0[index_1]  + exp(log_r[index_1]) .* time, exp(log_sigma[index_1]));
}