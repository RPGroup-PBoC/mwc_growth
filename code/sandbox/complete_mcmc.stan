

data {
    int<lower=1> N; // Number of data points.
    vector<lower=0>[N] I1; // Intensity of daughter 1
    vector<lower=0>[N] I2; // Intensity of daughter 2
    real<lower=0, upper=1> p;
}

parameters {
    real<lower=0> alpha; // Calibration factor. 
    real<lower=1E-9> sigma; // Homoscedastic error
    vector<lower=1>[N] N_tot; // Total number of proteins
    vector<lower=0, upper=N_tot[N]>[N] N1; // proteins in daugter 1
}

model {
   alpha ~ normal(0, 1000);
   N_tot ~ normal(0, 100);
   sigma ~ normal(0, 1);
   N1 ~ normal(N_tot * p, sqrt(N_tot * p * (1 - p)));
   I1 ~ normal(alpha * N1, sigma);
   I2 ~ normal(alpha * (N_tot - N1), sigma);
}












