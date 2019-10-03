data {
    // Dimensional information
    int<lower=1> J; // total number of ATC bins      
    int<lower=1> N; // total number of measurements

    // Identification vectors 
    int<lower=1, upper=J> idx[N]; // ID vector mapping point to atc_bin

    // Define the observables. 
    vector[N] foldchange; 
    vector<lower=1>[N] repressors;
}

parameters {
    // Define the parameters
    vector<lower=0, upper=1>[J] fc_mu;
    vector<lower=0>[J] fc_sigma;
    vector<lower=0, upper=1>[J] rep_mu;
    vector<lower=0>[J] rep_sigma;

}

model { 
    // Define the priors
    fc_mu ~ beta(1.1, 3);
    fc_sigma ~ normal(0, 0.1);
    rep_mu ~ lognormal(0, 5);
    rep_sigma ~ normal(0, 100);

    // Define the likelihood
    foldchange ~ normal(fc_mu[idx], fc_sigma[idx]);
    repressors ~ normal(rep_mu[idx], rep_sigma[idx]);
}

generated quantities {
    // Compute the empirical bohr parameter
    vector[J] empirical_F = - log(-1 + 1 ./ fc_mu);
}