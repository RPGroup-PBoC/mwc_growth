functions {
    // real marginalized_normal_lpdf(vector I1, vector I2, vector n1, vector n_tot,
    //  real alpha, int N)
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
    real gamma_appx_binom_lpdf(vector n_1, vector n_tot, int N) { 
        real lp = 0;
        for (i in 1:N) {
            if (n_1[i] > N) 
                return -1E10;  
            else
                lp += -n_tot[i] * log(2) + lgamma(n_tot[i] + 1) - lgamma(n_1[i]
                + 1) - lgamma(n_tot[i] - n_1[i] + 1);
        }
        return lp;
    }
}

}
data {
    int <lower=1> N; // Number of measurements
    vector<lower=0>[N] I1;
    vector<lower=1>[N] I2;
}


parameters {
   real<lower=0> alpha;
   real<lower=0> sigma;
//    vector<lower=0>[2*N] n_prot;
   vector<lower=0>[N] n_tot;
   vector<lower=0>[N] n_1;
}


model {
    alpha ~ normal(0, 500);
    n_tot ~ normal(0, 2000);
    n_1 ~ gamma_appx_binom(n_tot, N);
    sigma ~ normal(0, 100);
    I1 ~ normal(alpha * n_1, sigma);
    I2 ~ normal(alpha * (n_tot - n_1), sigma);

}