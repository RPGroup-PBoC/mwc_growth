data {
    // Dimensional parameters
    int<lower=1> J; // Number of unique colonies
    int<lower=1> N; // Total number of data points
    int<lower=1, upper=J> colony[N]; // Colony identifier
    
    // measured parameters 
    vector<lower=0>[N] fractional_area; // A(t) / A(t=0)
    vector[N] time; // In meaningful units
    }

parameters { 
    // Hyperparameters for lambda
    real<lower=0> lambda_mu;
    real<lower=0> lambda_sigma;

    // Hyperparameters for sigma
    real<lower=0> sigma_mu;
    real<lower=0> sigma_sigma;

    // Low-level parameters
    vector<lower=0>[J] lambda;
    vector<lower=0>[J] sigma;
    }

transformed parameters {
    // Convert to log for more effective sampling
    vector[N] log_area;
    log_area = log(fractional_area);
}

model {
    // Define a vector to compute the mean of the likelihood
    vector[N] mu;

    // Define hyperpriors
    lambda_mu ~ normal(0, 100);
    lambda_sigma ~ normal(0, 1);
    sigma_mu ~ normal(0, 1);
    sigma_sigma ~ normal(0, 1);

    // Define low level priors.
    lambda ~ cauchy(lambda_mu, lambda_sigma);
    sigma ~ normal(sigma_mu, sigma_sigma);

    // Evaluate the likelihood
    for (i in 1:N) {
        mu[i] = time[i] * lambda[colony[i]];
        log_area[i] ~ normal(mu[i], sigma[colony[i]]);
    }
}
