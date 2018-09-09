data  {
    // Dimensional parameters
    int<lower=1> J; // Number of unique samples
    int<lower=1> N; // Total number of measurements
    int<lower=1, upper=J> idx[N]; // Identification vector for samples

    // Physical parameters
    vector<lower=0>[N] time; // total growth time in units of minutes

    // Observed parameters
    vector<lower=0>[N] rel_abs; // Relative absorbance measurement at 600nm relative to t=0. 
}

parameters {
    // Hyper parameters
    real<lower=0> lambda_mu; // Hyperparameter for growth rate
    real<lower=0> lambda_sigma; // Hyperparameter for homoscedastic error

    // Low level parameters
    real<lower=0> lambda[J]; // Growth rate for individual samples
    real<lower=0> sigma[J]; // Homoscedastic error for samples
}

transformed parameters {
    vector<lower=0>[N] log_rel_abs; // Log absorbance ratio at 600nm relative to t=0
    log_rel_abs = log(rel_abs);
}

model {
    // Define a vector to compute the theoretical value
    vector[N] mu; 

    // Define the priors. 
    lambda_mu ~ normal(0, 1);
    lambda_sigma ~ normal(0, 1);
    lambda ~ normal(lambda_mu, lambda_sigma);
    sigma ~ normal(0, 1);

    // Evaluate the likelihood
    for (i in 1:N) {
        mu[i] = time[i] * lambda[idx[i]];
        log_rel_abs[i] ~ normal(mu[i], sigma[idx[i]]);
    }
}