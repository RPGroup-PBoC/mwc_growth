functions {
    // real marginalized_normal_lpdf(vector I1, vector I2, vector n1, vector n_tot,
    // real alpha, int N)
    // {
    //     vector[N] aI1 = alpha * n1;
    //     vector[N] aI2 = alpha * (n_tot - n1);
    //     real out = 0;
    //     real logpost = 0;
    //     real binom = 0;
    //     for (i in 1:N) {
    //     if (n1[i] > n_tot[i])
    //         return  -1E10;
    //     else
    //         out += (I1[i] - aI1[i])^2 + (I2[i] - aI2[i])^2;
    //         binom += lgamma(n_tot[i] + 1) - lgamma(n1[i] + 1) - lgamma((n_tot[i]-n1[i]) + 1) - n_tot[i] * log(2);
    //     }
    //     return binom - N * 0.5 * log(out);
    // }
    real gamma_appx_binom_lpdf(vector n_1, vector n_tot) { 
        return sum(lgamma(n_tot + 1) - lgamma(n_1 + 1) - lgamma(n_tot - n_1 +1)
        - n_tot * log(2));
    }
}
data {
    int <lower=1> N; // Number of measurements
    vector<lower=0>[N] I1;
    vector<lower=1>[N] I2;
}


parameters {
   real<lower=0> log_alpha;
   real<lower=0> sigma;
   vector<lower=0>[2*N] n_prot;
}

transformed parameters {
    real alpha = exp(log_alpha);

}


model {
    n_prot[:N] ~ normal(0, 2000);
    log_alpha ~ normal(0, 5);
    sigma ~ normal(0, 1);
    n_prot[N+1:] ~ gamma_appx_binom(n_prot[:N]);
    I1 ~ normal(alpha *n_prot[N+1:], sigma);
    I2 ~ normal(alpha * (n_prot[:N] - n_prot[N+1:]), sigma);

}