functions {
    real marginalized_normal_lpdf(vector I1, vector I2, vector n1, vector n_tot,
    real alpha, int N)
    {
        vector[N] aI1 = alpha * n1;
        vector[N] aI2 = alpha * (n_tot - n1);
        real out = 0;
        real logpost = 0;
        real binom = 0;
        for (i in 1:N) {
        if (n1[i] > n_tot[i])
            return  -1E10;
        else
            out += (I1[i] - aI1[i])^2 + (I2[i] - aI2[i])^2;
            binom += lgamma(n_tot[i] + 1) - lgamma(n1[i] + 1) - lgamma((n_tot[i]-n1[i]) + 1) - n_tot[i] * log(2);
        }
        return binom - N * 0.5 * log(out);
    }
}

data {
    int<lower=1> J; // Number of replicates
    int <lower=1> N; // Number of measurements

    vector<lower=0>[N] I1;
    vector<lower=1>[N] I2;
}

parameters {
   // Define the hyper parameter
   real log_alpha_hyper_mu;
   vector[J] log_alpha_mu_tilde;     
   real tau;
   vector<lower=0>[2*N] n_prot;
}

transformed parameters {
    vector[J] log_alpha_1 = log_alpha_hyper_mu + log_alpha_mu_tilde * tau;
    vector[J] alpha_1 = exp(log_alpha_1);
    real alpha_hyper = exp(log_alpha_hyper_mu);
}

model {
    log_alpha ~ normal(0, 4);
    n_prot ~ normal(0, 1500);
    I1 ~ marginalized_normal(I2, n_prot[:N], n_prot[N+1:], alpha, N);

}