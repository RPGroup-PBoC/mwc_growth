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
    int<lower=1> J_1; // Number of biological replicates
    int<lower=1> J_2; // Number of microcolonies
    int<lower=1> N; // Number of area measurements
    int<lower=1, upper=J_1> index_1[J_2]; // ID vector of colonies to days
    int<lower=1, upper=J_2> index_2[N];  // ID vector of measurements to colonies
    vector<lower=0>[N] area; // Measured cell area in units of square microns
    vector<lower=0>[N] time; // Time in units of minutes
    }
   
transformed data {
vector[N] log_area = log(area);
}

parameters {
    real log_r;
    real log_sigma;
    real<lower=0> tau_r;
    real<lower=0> tau_sigma;
    vector<lower=0>[J_2] area_0;
    vector[J_1] log_r_2_raw;
    vector[J_2] log_r_3_raw;
    vector[J_1] log_sigma_2_raw;
    vector[J_2] log_sigma_3_raw;
}

transformed parameters {
    vector[J_1] log_r_2 = log_r + tau_r * log_r_2_raw;
    vector[J_1] log_sigma_2 = log_sigma + tau_sigma * log_sigma_2_raw;
    vector[J_2] log_r_3 = log_r_2[index_1] + tau_r * log_r_3_raw;
    vector[J_2] log_sigma_3 = log_sigma_2[index_1] + tau_sigma * log_sigma_3_raw;
    vector[J_1] sigma_2 = exp(log_sigma_2); 
    vector[J_2] sigma_3 = exp(log_sigma_3); 
    vector[J_1] r_2 = exp(log_r_2);
    vector[J_2] r_3 = exp(log_r_3);
    real r = exp(log_r);
}

model {
    vector[N] mu  = log(area_0[index_2])  + time ./ r_3[index_2];
    r ~ normal(0, 1);
    log_sigma ~ normal(0, 2);
    area_0 ~ normal(0, 10);
    log_sigma_2_raw ~ normal(0, 1);
    log_sigma_2_raw ~ normal(0, 1);
    log_r_2_raw ~ normal(0, 1);
    log_r_3_raw ~ normal(0, 1);
    tau_r ~ lognormal(0, 2);
    tau_sigma ~ lognormal(0, 2);
    log_area ~ normal(mu, sigma_3[index_2]);    
}


