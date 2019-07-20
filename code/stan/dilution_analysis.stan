functions{
    /** 
    * Approximate the Binomial distirubution for continuous variables 
    * as a ratio of Gamma functions 
    * 
    * @param I1: Observed fluorescence of daughter cell 1. 
    * @param I2: Observed fluorescence of daughter cell 2.
    * @param alpha: Fluorescenc calibration factor in units of a.u. / molecule
    * @param N: Total number of measurements 
    **/
    real GammaApproxBinom_lpdf(vector I1, vector I2, real alpha, int N)
    {
        vector[N] n1 = I1 ./ alpha;
        vector[N] n2 = I2 ./ alpha;
        vector[N] ntot = n1 + n2;
        return -N * log(alpha) + sum(lgamma(ntot + 1) - 
                    lgamma(n1 + 1) - lgamma(n2 + 1) - 
                    ntot * log(2));
    } 
}
     
data {
    int<lower=0> N_fluct; // Number of data points
    int<lower=0> J; 
    int<lower=0> N_fc;
    int<lower=0, upper=J> fc_idx[N_fc];
    vector<lower=0> yfp[N_fc]
    vector<lower=0> mch[N_fc]
    vector<lower=0>[N_fluct] I1; // Observed fluorescence of daughter cell 1
    vector<lower=0>[N_fluct] I2; // Observed fluorescence of daughter cell 2
}
parameters {
    real alpha;
    vector<lower=0>[J] rep_mu;
    vector<lower=0>[J] rep_sigma;
    vector<lower=0>[J-2] fc_mu;
    vector<lower=0>[J-2] fc_sigma;
}

model {   
    alpha ~ normal(0, 500);
    I1 ~ GammaApproxBinom(I2, alpha, N);  

}
