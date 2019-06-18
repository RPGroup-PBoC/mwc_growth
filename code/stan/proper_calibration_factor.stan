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
   vector<lower=0>[N] n_tot;
   vector<lower=0>[N] n1;
   real<lower=0> sigma;
}

// transformed parameters {
//     vector<lower=0>[N] n_tot_round = floor(n_tot);
//     vector<lower=0>[N] n1_round = floor(n1);
//     vector<lower=0>[N] n2 = n_tot_round - n1_round;
// }

model {
    n_tot ~ normal(0, 1000); 
    // n1_round ~ approx_binom(n_tot_round);
    n1 ~ normal(n_tot * 0.5, n_tot * (0.5 * 0.5));
    sigma ~ normal(0, 10);
    I1 ~ normal(n1 * alpha, sigma);
    I2 ~ normal((n_tot - n1) * alpha, sigma);
}