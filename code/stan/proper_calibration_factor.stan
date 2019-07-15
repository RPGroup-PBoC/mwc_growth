functions {
    real approx_binom_lpdf(vector n1, vector n_tot) {
        return sum(lgamma(n_tot + 1) - lgamma(n1 + 1) 
                - lgamma(n_tot - n1 +1) - n_tot * log(2));
}
}
data {
    int <lower=1> N; // Number of measurements
    vector<lower=0>[N] I1;
    vector<lower=1>[N] I2;
}

parameters {
   real<lower=0> alpha;     
   positive_ordered[2*N] n_prot;
   real<lower=0> sigma;
}

model {
    n_prot[N+1:] ~ normal(0, 5000); 
    n_prot[:N] ~ approx_binom_lpdf(n_prot[N+1:]) ;
    alpha ~ normal(0, 1000);
    sigma ~ normal(0, 10);
    I1 ~ normal(n_prot[:N] * alpha, sigma);
    I2 ~ normal((n_prot[N+1:] - n_prot[:N]) * alpha, sigma);
}