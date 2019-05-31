


data { 
    // Dimensional parameters
    int<lower=1> J; // Number of repressor copy numbers (grouped to integers)
    int<lower=1> N_dilution; // Number of measurements of dilution strain
    int<lower=1> N_auto; // Number of autofluorescence measurements
    int<lower=1> N_delta; // Number of constitutive expression measurements
    int<lower=1, upper=J> idx[N_dilution]; // Repressor idx

    // Fluorescence measurements.
    vector<lower=0>[N_auto] auto_yfp;
    vector<lower=0>[N_delta] delta_yfp;
    vector<lower=0>[N_dilution] dilution_yfp;
}

parameters { 
    real<lower=0> auto_mu;
    real<lower=0> auto_sig;
    real<lower=0> delta_mu;
    real<lower=0> delta_sig;
    vector<lower=0>[J] dilution_mu;
    vector<lower=0>[J] dilution_sig;
}

model { 
    auto_mu ~ uniform(0, 4095);
    auto_sig ~ normal(0, 10);
    delta_mu ~ uniform(0, 4095);
    delta_sig ~ normal(0, 10);
    dilution_mu ~ uniform(0, 4095);
    dilution_sig ~ normal(0, 10);

    auto_yfp ~ normal(auto_mu, auto_sig);
    delta_yfp ~ normal(delta_mu + auto_mu, delta_sig);
    dilution_yfp[idx] ~ normal(dilution_mu[idx] + auto_mu, dilution_sig[idx]);
}

generated quantities {
    vector[J] fold_change = dilution_mu ./ delta_mu;
}