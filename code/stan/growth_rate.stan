data  {
    int<lower=1> J;
    int<lower=1> N;
    int<lower=1, upper=J> idx[N];
    vector<lower=0>[N] A;
    vector<lower=0>[N] time;
    }
    
transformed data {
    vector[N] log_a = log(A);
    }

parameters {
      real<lower=0> r;
      real log_sigma;
      vector[J] r_raw;
      vector[J] log_sigma_raw;
      real<lower=0> tau_r;
      real<lower=0> tau_sigma;
}

transformed parameters {
    vector[J] r_1 = r + r_raw * tau_r;
    vector[J] log_sigma_1 = log_sigma + log_sigma_raw * tau_sigma; 
    vector[J] sigma_1 = exp(log_sigma_1);
}

model {
    r ~ normal(0, 1);
    log_sigma ~ normal(0, 1);
    log_sigma_raw ~ normal(0, 1);
    r_raw ~ normal(0, 1);
    tau_r ~ normal(0, 1);
    tau_sigma ~ normal(0, 1);
    log_a ~ normal(r_1[idx] .* time, sigma_1[idx]);
}