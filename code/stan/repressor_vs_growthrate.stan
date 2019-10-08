

data {
    // Dimensional info
    int<lower=1> J; // Number of unique ATC concentrations
    int<lower=1> N; // Number of measurements
    int<lower=1, upper=J> idx[N]; // Maps the induction condition

    // Define the observables
    vector<lower=0, upper=1>[N] growth_rate; // growth rate in hr^-1
    vector<lower=0>[N] repressors; // Number of repressors per cell 
}

parameters { 
    real slope; 
    real<lower=0> intercept[J]; 
    real<lower=0> sigma; // Homoscedstic error about the trend
}

model { 
    // Define a vector to compute the mean of thenormal distribution
    vector[N] mu;

    // Define the priors
    slope ~ normal(0, 1); 
    intercept ~ normal(0, 100);

    
    for (i in 1:N) { 
        mu[i] =  intercept[idx[i]] + slope * growth_rate[i];
    }

    repressors ~ normal(mu, sigma);
}
