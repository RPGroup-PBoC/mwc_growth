data {
    int<lower=1> N;
    int<lower=1> J;
    vector<lower=1>[N] noise;
    vector<lower=1>[J] signal;
}

parameters {
    real<lower=0> mu1;
    real<lower=0> sigma1;
    real<lower=0> mu2;
    real<lower=0> sigma2;
}

model {
    mu1 ~ normal(0, 500);
    sigma1 ~ normal(0, 100);
    mu2 ~ normal(0, 500);
    sigma2 ~ normal(0, 100);
    noise ~ normal(mu1, sigma1);
    signal ~ normal(mu1 + mu2, sigma2);
}