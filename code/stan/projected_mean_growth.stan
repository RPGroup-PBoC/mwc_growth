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
    vector[J] log_A0;
    vector[J] log_r;
    vector[J] log_sigma;
}

transformed parameters {  
    vector[J] sigma = exp(log_sigma);
    vector[J] r = exp(log_r);
    vector[J] A0 = exp(log_A0);    
}

model {
    // Set the priors
    log_A0 ~ normal(0, 1);
    log_r ~ normal(0, 2);
    log_sigma ~ normal(0, 2); 
    log_area ~ normal(log_A0[index_1] + time ./ r[index_1], sigma[index_1]);
}