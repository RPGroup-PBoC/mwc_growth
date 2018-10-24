/* 
* Hierarchical Model for Growth Rate Determination
* ------------------------------------------------
* Author: Griffin Chure
* License: MIT
*
* Description
* -------------------------------------------------
* Defines a hierarchical model for growth rate determined
* via microscopy for J different media, each taken on D 
* different days, with R individual colonies
*/

data {
    int<lower=1> J;
    int<lower=1> N;
    int<lower=1, upper=J> idx[N];
    vector<lower=0>[N] area;
    vector<lower=0>[N] time;
    }
   
transformed data {
vector[N] log_area = log(area);
}

parameters {
    real r;
    real log_sigma;
    real<lower=0> tau_r;
    real<lower=0> tau_sigma;
    vector[J] area_0;
    vector[J] r_2_raw;
    vector[J] log_sigma_2_raw;
}

transformed parameters {
    vector[J] r_2 = r + tau_r * r_2_raw;
    vector[J] log_sigma_2 = log_sigma + tau_sigma * log_sigma_2_raw;
    vector[J] sigma_2 = exp(log_sigma_2); 
}

model {
    vector[N] mu  = area_0[idx] .* exp(r_2[idx] .* time);
    r ~ lognormal(-1, 3);
    log_sigma ~ normal(0, 2);
    area_0 ~ normal(0, 100);
    log_sigma_2_raw ~ normal(0, 1);
    r_2_raw ~ normal(0, 1);
    tau_r ~ normal(0, 1);
    tau_sigma ~ normal(0, 10);
    area ~ normal(mu, sigma_2[idx]);    
}
