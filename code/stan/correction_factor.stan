data { 
    int<lower=1> N; // Number of measurements
    vector<lower=0>[N] R; // Inferred number of repressors
    real epRA; 
    real epAI;
    real Nns;
    vector[N] foldchange;
}

transformed data {
    vector[N] logfc = log(foldchange);
}

parameters { 
    real corr_factor;
    real<lower=0> sigma;
}

transformed parameters { 
    vector[N] corr_R= corr_factor * R;
}

model { 
    vector[N] pact = corr_R / (1 + exp(-epAI));

    vector[N] mu = -log( 1 + pact * exp(-epRA) / Nns);
    corr_factor ~ normal(0, 10);
    sigma ~ normal(0, .1);
    logfc ~ normal(mu, sigma);
}