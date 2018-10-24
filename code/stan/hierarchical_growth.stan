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
    real log_r;
    real log_sigma;
    real<lower=0> tau_r;
    real<lower=0> tau_sigma;
    vector<lower=0>[J] area_0;
    vector[J] log_r_2_raw;
    vector[J] log_sigma_2_raw;
}

transformed parameters {
    vector[J] log_r_2 = log_r + tau_r * log_r_2_raw;
    vector[J] log_sigma_2 = log_sigma + tau_sigma * log_sigma_2_raw;
    vector[J] sigma_2 = exp(log_sigma_2); 
    vector[J] r_2 = exp(log_r_2);
    real r = exp(log_r);
}

model {
    vector[N] mu  = log(area_0[idx])  + r_2[idx] .* time;
    r ~ normal(0, 1);
    log_sigma ~ normal(0, 2);
    area_0 ~ normal(0, 100);
    log_sigma_2_raw ~ normal(0, 1);
    log_r_2_raw ~ normal(0, 1);
    tau_r ~ normal(0, 2);
    tau_sigma ~ normal(0, 2);
    log_area ~ normal(mu, log_sigma_2[idx]);    
}

generated quantities {
    real lambda = 1 / r;
}
