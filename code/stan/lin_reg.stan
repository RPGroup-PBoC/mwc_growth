data {
    int<lower=0> N; 
    vector[N] mean_vals;
    vector[N] temps;
}

parameters {
    real<lower=0> slope;
    real<upper=0> intercept;
    real<lower=0> sigma;
}

model {
    // slope ~ normal(10, 10);
    // intercept ~ normal(10, 10);
    sigma ~ normal(0, 0.1);
    mean_vals ~ normal((slope ./ temps) - intercept, sigma);
}