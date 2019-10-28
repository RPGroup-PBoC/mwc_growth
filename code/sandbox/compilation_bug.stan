data {
    int N;
    vector[N] theta;
    vector[N] y;
}

parameters {
    real alpha[N];
    real<lower=0> sigma;
}

transformed parameters { 
    vector[N] mu = alpha .* theta;
}

model { 
    alpha ~ normal(0, 1);
    sigma ~ normal(0, 1);
    y ~ normal(mu, sigma);
}

