data {
    // Define dimensional parameters
    int<lower=1> J; // Number of unique conditions
    int<lower=1, upper=J> J_REF; // Index of the reference state
    int<lower=1> N; // Total number of fold-change measurements
    int<lower=1, upper=J> idx[N]; // Identification vector for each point. 

    // Define input constants. 
    real<lower=1> Nns; // Number of nonspecific binding sites. 
    vector[J] temp; // Temperature vector in K 
    vector<lower=0>[N] repressors; // Repressor count per replicate

    // Define the observed quantity
    vector[N] foldchange; // Observed fold-change in gene expression
}

transformed data { 
    vector[N] log_fc = log(foldchange);
}

parameters { 
    real delta_SR;
    real ref_epRA;
    real<lower=0> sigma; // Homoscedastic error
}

transformed parameters { 
    // Compute the modified binding energies;
    vector[J] epRA_star;
    for (i in 1:J) {
        if (i == J)
            epRA_star[i] = ref_epRA;
        else 
            epRA_star[i] = delta_SR * (temp[J_REF] - temp[i]) + ref_epRA; 
       
    }
}

model { 

     // Compute the mean fold-change in gene expression
     vector[N] mu = -1 * log(1 + (repressors ./ Nns) .* exp(-epRA_star[idx])); 

     // Define the priors
     ref_epRA ~ normal(-12, 10);
     delta_SR ~ normal(0, 10);
     sigma ~ normal(0, 0.1);

     // Evaluate the likelihood
     log_fc ~ cauchy(mu, sigma);

 }
