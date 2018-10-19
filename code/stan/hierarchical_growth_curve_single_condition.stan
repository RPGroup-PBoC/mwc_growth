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
    int<lower=1> J; // Number of unique growth media
    int<lower=1> D; // Total number of days 
    int<lower=1> R; // Total number of microcolonies
    int<lower=1> N; // Number of individual measurements
    int<lower=1> M; // Number of unique parameters

    // Identifier arrays for measurements
    int<lower=1, upper=J> media_idx[N];
    int<lower=1, upper=D> date_idx[N];
    int<lower=1, upper=R> rep_idx[N]; 

    // Identifier arrays for parameters
    int<lower=1, upper=J> unique_media_idx[M];
    int<lower=1, upper=D> unique_date_idx[M];
    int<lower=1, upper=R> unique_rep_idx[M];

    // Observed parameters
    int<lower=0> time[N]; // in minutes
    int<lower=0> area[N]; // in square pixels
}


parameters {
    // Level-0 parameters (per medium)
    real<lower=0> lambda_0[J]; 
    real<lower=0> sigma_0[J];

    // Level-1 parameters (per day)
    real<lower=0> lambda_1[D];
    real<lower=0> sigma_1[D];

    // Level-2 parameters (per cell)
    real<lower=0> lambda_2[R];
    real<lower=0> sigma[R]; // Homoscedastic error
    real<lower=0> area_0[R];
}

model {
    real theo;
    // Assign level-0 prior
    lambda_0 ~ normal(0, 10);
    sigma_0 ~ normal(0, 10);

    // Aissign static level-1 and 2 priors.
    sigma_1 ~ normal(0, 10);
    area_0 ~ lognormal(0, 2);
    sigma ~ normal(0, 1);

    // Assign level-1 and level-2 priors 
    for (i in 1:M) { 
        lambda_1[unique_date_idx[i]] ~ normal(lambda_0[unique_media_idx[i]], sigma_0[unique_media_idx[i]]);
        lambda_2[unique_rep_idx[i]] ~ normal(lambda_1[unique_date_idx[i]], sigma_1[unique_date_idx[i]]);
    }

    // Evaluate the likelihood
    for (i in 1:N) {
        theo = time[i] / lambda_2[rep_idx[i]];
        log(area[i] / area_0[rep_idx[i]]) ~ normal(theo, sigma[rep_idx[i]]);
    }
}


