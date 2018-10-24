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
    vector[J] log_area_0;
    vector[J] log_r_2;
    vector[J] log_sigma_2;
}

transformed parameters {
    real r = exp(log_r);
    real sigma = exp(log_sigma);
    vector[J] r_2 = exp(log_r_2);
    vector[J] sigma_2 = exp(log_sigma_2); 
}

model {
    vector[N] mu  = log_area_0[idx] + r_2[idx] .* time;
    log_r ~ normal(0, 1);
    sigma ~ normal(0, 1);
    log_sigma_2 ~ normal(0, 1);
    log_r_2 ~ normal(log_r, log_sigma);
    log_area ~ normal(mu, log_sigma_2[idx]);    
}
