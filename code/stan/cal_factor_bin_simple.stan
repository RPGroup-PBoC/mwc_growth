data {
    int<lower=1> J; // Number of bins
    int<lower=1> N; // Number of measurements
    int<lower=1, upper=J> idx[N]; // ID vector for bins
    vector<lower=0>[N] summed;
    vector<lower=0>[N] fluct;
}

transformed data {
    vector[N] log_summed = log(summed);
    vector[N] log_fluct = log(fluct);

}


parameters {
    vector[J] summed_mu;
    vector<lower=0>[J] summed_sigma;
    vector[J] fluct_mu;
    vector<lower=0>[J] fluct_sigma;
    real log_alpha;
    real<lower=0> sigma;
}

transformed parameters {
    real alpha = exp(log_alpha);
}

model {
    summed_mu ~ normal(0, 100);
    summed_sigma ~ normal(0, 10);
    fluct_mu ~ normal(0, 100);
    fluct_sigma ~ normal(0, 10);
    log_alpha ~ normal(0, 4);
    sigma ~ normal(0, 1);
    // Likelihood
    log_summed ~ normal(summed_mu[idx], summed_sigma[idx]);
    log_fluct ~ normal(fluct_mu[idx], fluct_sigma[idx]);
    log(fluct_mu) ~ normal(log_alpha + log(summed_mu, sigma);
}