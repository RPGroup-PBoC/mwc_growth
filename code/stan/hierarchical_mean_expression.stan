/*
* Hierarchical Estimation of Mean Expression
* ------------------------------------------
* Author: Griffin Chure
* License: MIT
*/ 

data {
    int<lower=1> J; // Number of unique experiments
    int<lower=1> N; // Total number of measurements
    int<lower=1> idx[N]; // Identification vector for experiment number
    real mCherry[N];
    real YFP[N];
}

parameters {
    // Top level parameters
    real<lower=0> log_mCherry_1;
    real<lower=0> log_YFP_1;
    real log_mCherry_sigma_1;
    real log_YFP_sigma_1;

    // Low-level uncentered parameters
    vector[J] log_mCherry_2_raw;
    vector[J] log_mCherry_sigma_2_raw;
    vector[J] log_YFP_2_raw;
    vector[J] log_YFP_sigma_2_raw;

    // How the hyperparameters vary
    real<lower=0> tau_log_mCherry;
    real<lower=0> tau_log_YFP;
    real<lower=0> tau_log_mCherry_sigma;
    real<lower=0> tau_log_YFP_sigma;
}

transformed parameters {
    // Un-center the parameter distributions
    vector[J] log_mCherry_2 = log_mCherry_1 + tau_log_mCherry * log_mCherry_2_raw;
    vector[J] log_mCherry_sigma_2 = log_mCherry_sigma_1 + tau_log_mCherry_sigma * log_mCherry_sigma_2_raw;
    vector[J] log_YFP_2 = log_YFP_1 + tau_log_YFP * log_YFP_2_raw;
    vector[J] log_YFP_sigma_2 = log_YFP_sigma_1 + tau_log_YFP_sigma * log_YFP_sigma_2_raw;

    // Convert parameters back to linear space
    vector[J] mCherry_2 = exp(log_mCherry_2);
    vector[J] mCherry_sigma_2 = exp(log_mCherry_sigma_2);
    vector[J] YFP_2 = exp(log_YFP_2);
    vector[J] YFP_sigma_2 = exp(log_YFP_sigma_2);
    real mCherry_1 = exp(log_mCherry_1);
    real YFP_1 = exp(log_YFP_1);
}

model {
    // Define hyperpriors
    log_mCherry_1 ~ normal(0, 3);
    log_YFP_1 ~ normal(0, 1);
    log_mCherry_sigma_1 ~ normal(0, 3);
    log_YFP_sigma_1 ~ normal(0, 1);

    // Define priors on uncentering parameters. 
    tau_log_mCherry ~ normal(0, 2);
    tau_log_YFP ~ normal(0, 1);
    tau_log_mCherry_sigma ~ normal(0, 2);
    tau_log_YFP_sigma ~ normal(0, 1);
    log_mCherry_2_raw ~ normal(0, 1);
    log_YFP_2_raw ~ normal(0, 1);
    log_mCherry_sigma_2_raw ~ normal(0, 1);
    log_YFP_sigma_2_raw ~ normal(0, 1);

    // Evaluate the likellihood
    mCherry ~ normal(mCherry_2[idx], mCherry_sigma_2[idx]);
    YFP ~ normal(YFP_2[idx], YFP_sigma_2[idx]);
}
