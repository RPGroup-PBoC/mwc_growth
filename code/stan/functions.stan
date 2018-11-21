/*
* Functions
* -------------------------------------
* Author: Griffin Chure
* License: MIT
*
* Description 
* ------------------------------------
* This file contains an array of functions
* used for the inference used in the 
* MWC Growth Project, github.com/gchure/mwc_growth
*/

functions{
    /** 
    * Approximate the Binomial distirubution for continuous variables 
    * as a ratio of Gamma functions 
    * 
    * @param I1: Observed fluorescence of daughter cell 1. 
    * @param I2: Observed fluorescence of daughter cell 2.
    * @param alpha: Fluorescence calibration factor in units of a.u. / molecule
    * @param N: Total number of measurements 
    **/
    real GammaApproxBinom_lpdf(vector I1, vector I2, vector alpha ) { 
            return sum(-log(alpha))  + sum(lgamma(((I1 + I2) ./ alpha) + 1) - lgamma((I1 ./ alpha) + 1)
                        - lgamma((I2 ./ alpha) + 1) - ((I1 + I2) ./ alpha) * log(2));
        }
        
    /**
    * Calculate the fold-change in gene expression.
    *
    * @param R The number of repressors per cell
    * @param Nns The number of nonspecific repressor binding sites.
    * @param ep_r The  binding energy of the repressor to the DNA in kBT.
    * @param ep_ai The energetic difference between the active and inactive
    *        states of the repressor in kBT.
    **/
    vector fold_change(vector R, real Nns, real ep_r, real ep_ai) {
        // Compute the various componenets piecewise for simplicity.
        real pact = 1 / (1 + exp(-ep_ai));
        return 1 ./ (1 + pact * (R / Nns) * exp(-ep_r));
        } 
}