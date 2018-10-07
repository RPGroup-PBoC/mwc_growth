/* 
* Calibration Factor Estimation
* ---------------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description 
* ---------------------------------------------------
* This model samples the posterior probability distribution
* for the calibration factor between observed fluorescence 
* and fluorophore copy number. This model assumes that 
* the noise in measurement is negligible compared to the 
* observed fluorescence, making individual fluorescence 
* measurements behave as delta functions. 
*/
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
    real GammaApproxBinom_lpdf(real[] I1, real[] I2, real alpha, int N) {
        vector[N] lprob;
        for (i in 1:N){
            lprob[i] = lgamma(((I1[i] + I2[i]) / alpha) + 1) - lgamma((I1[i] / alpha) + 1)
                        - lgamma((I2[i] / alpha) + 1) - ((I1[i] + I2[i]) / alpha) * log(2);
        }
        return -N * log(alpha) + sum(lprob);
    }
     
}
data {
    int<lower=0> N; // Number of data points
    real I1[N]; // Observed fluorescence of daughter cell 1
    real I2[N]; // Observed fluorescence of daughter cell 2
}

parameters {
    real<lower=0> alpha; // Calibraiton factor in units of a.u. /molecule
}

model {    
    alpha ~ lognormal(0, 10);
    I1 ~ GammaApproxBinom(I2, alpha, N);  
}
