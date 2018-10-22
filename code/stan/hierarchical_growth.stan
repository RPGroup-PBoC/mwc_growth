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

data  {
    //Dimensional parameters
    int<lower=1> J_1; // Number of unique media
    int<lower=1> J_2; // Number of unique days
    int<lower=1> J_3; // Number of unique colonies
    int<lower=1> N; // Number of measurements in data set
    
    // Indices for hierarchy
    int<lower=1, upper=J_1> index_1[J_2];
    int<lower=1, upper=J_2> index_2[J_3];
    int<lower=1, upper=J_3> index_3[N];

    // Experimental measurements
    vector<lower=0>[N] time; // in minutes
    vector<lower=0>[N] area; // in µm^2
    }

parameters {
    // Centered parameters
    vector<lower=0>[J_1] lambda; 
    vector<lower=0>[J_1] log_sigma;
    vector<lower=0>[J_3] area_mu;
    vector[J_3] area_raw;

    // How the hyperparameters vary
    real<lower=0> tau_lambda;
    real<lower=0> tau_sigma;
    real<lower=0> tau_area;

    // Non-centered parameters.
    vector[J_2] lambda_2_raw;
    vector[J_2] log_sigma_2_raw;
    vector[J_3] lambda_3_raw;
    vector[J_3] log_sigma_3_raw;
}

transformed parameters {
    // Level-1  parameters
    vector[J_2] lambda_2 = lambda[index_1] + tau_lambda * lambda_2_raw;
    vector[J_2] log_sigma_2 = log_sigma[index_1] + tau_sigma * log_sigma_2_raw;

    // Level 2 parameters
    vector[J_3] lambda_3 = lambda_2[index_2] + tau_lambda * lambda_3_raw;
    vector[J_3] log_sigma_3 = log_sigma_2[index_2] + tau_sigma * log_sigma_3_raw;
    vector[J_3] area0 = area_mu[index_2] + area_raw * tau_area;
   
    vector[J_1] sigma = exp(log_sigma);
    vector[J_2] sigma_2 = exp(log_sigma_2);
    vector[J_3] sigma_3 = exp(log_sigma_3);
    
}

model {
    vector[N] mu;
    
    // Define priors for uncentering offsets
    tau_lambda ~ normal(0, 1);
    tau_sigma ~ normal(0, 1);
    tau_area ~ normal(0, 1);

    // Define priors for uncentered parameters
    lambda_2_raw ~ normal(0, 1);
    log_sigma_2_raw ~ normal(0, 2);
    lambda_3_raw ~ normal(0, 1);
    log_sigma_3_raw ~ normal(0, 2);
    area_raw ~ normal(0, 2);

    // Define cenetered parameters
    area_mu ~ normal(0, 2);
    lambda ~ normal(0, 10);
    
    mu = area0[index_3] .* exp(time ./ lambda_3[index_3]);
    area ~ normal(mu, sigma_3[index_3]);
}