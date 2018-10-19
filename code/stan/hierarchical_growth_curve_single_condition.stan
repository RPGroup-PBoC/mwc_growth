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
    //Dimensional parameters
    int<lower=1> J_1; // Number of unique growth media
    int<lower=1> J_2; // Total number of days 
    int<lower=1> J_3; // Total number of microcolonies
    int<lower=1> N; // Number of individual measurements

    // Identifier arrays for measurements
    int<lower=1, upper=J_1> idx_1[J_2];
    int<lower=1, upper=J_2> idx_2[J_3];
    int<lower=1, upper=J_3> idx_3[N]; 

    // Observed parameters
    vector<lower=0>[N] time; // in minutes
    vector<lower=0>[N] area; // in square microns
}


parameters {
    // Level-0 parameters (per medium)
    vector<lower=0>[J_1] lambda_1; 
    real<lower=0> sigma_1;

    // Level-1 parameters (per day)
    vector[J_2] theta_2;
    real<lower=0> tau_2;

    // Level-2 parameters (per cell)
    vector[J_3]  theta_3;
    real<lower=0> tau_3;
    real<lower=0> sigma[J_3]; // Homoscedastic error
    vector[J_3] area_0;
}

transformed parameters {
    // Generate uncentered versions of the parameters
    vector<lower=0>[J_2] lambda_2 = lambda_1[idx_1] + tau_2 * theta_2[idx_1];
    vector<lower=0>[J_3] lambda_3 = lambda_2[idx_2] + tau_3 * theta_3[idx_2];
}


model {
    real theo;
    // Assign level-0 prior
    lambda_1 ~ lognormal(2, 3);
    sigma_1 ~ normal(0, 10);

    // Aissign static level-1 and 2 priors.
    tau_2 ~ normal(0, 1);
    tau_3 ~ normal(0, 1);
    sigma ~ normal(0, 10);

    // Evaluate the likelihood
    for (i in 1:N) {
    area[i] ~ normal(area_0[idx_3[i]] * exp(time[i] / lambda_3[idx_3[i]]), sigma[idx_3[i]]);
    }
    
}


