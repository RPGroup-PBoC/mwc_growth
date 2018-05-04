/**
* Ka/Ki Fit Model
* --------------------------------------------------------------------------
* This Stan model fits the log transform of the dissociation constants Ka and
* Ki for the active and inactive state of the repressor, respectively.
*
* Author: Griffin Chure
* Creation Date: March 11, 2018
* License: MIT
* Contact: gchure@caltech.edu
**/

functions {
    /**
    * Compute the probability of a repressor being active given an inducer
    * concentration c.
    *
    * @param c Concentration of allosteric effector.
    * @param ep_a Log transform of effector dissociation constant from active
    *        repressor, Ka, in kBT.
    * @param ep_i Log transform of effector dissociation constant from inactive
    *        repressor, Ki, in kBT.
    * @param ep_ai Energy difference between the active and inactive state of
    *        the repressor in kBT.
    * @param n_sites The number of allosterically dependent sites.
    * @return prob_act The probability of a repressor being active with the
    *         given parameters.
    **/
    real prob_act(real c, real ep_a, real ep_i, real ep_ai, int n_sites) {
        // Calculate the relevant components piecewise for simplicity.
        real numerator;
        real denominator;
        numerator = (1 + c * exp(-ep_a))^n_sites;
        denominator = numerator + exp(-ep_ai) * (1 + c * exp(-ep_i))^n_sites;
        return numerator / denominator;}

    /**
    * Compute the level of repression in a simple repression architecture.
    *
    * @param pact The probability of an active repressor.
    * @param R The number of repressors per cell.
    * @param Nns The number of nonspecific binding sites.
    * @param ep_r The binding energy of the repressor to the DNA in kBT.
    * @return repression The level of repression given these parameters.
    **/
    real repression(real pact, real R, real Nns, real ep_r) {
        return 1 + pact * (R / Nns) * exp(-ep_r);
      }

    /**
    * Calculate the fold-change in gene expression.
    *
    * @param R The number of repressors per cell
    * @param Nns The number of nonspecific repressor binding sites.
    * @param ep_r The  binding energy of the repressor to the DNA in kBT.
    * @param c The concentration of allosteric effector.
    * @param ep_a The log transform of the effector dissociation constant from
    *        the active repressor, Ka, in kBT.
    * @param ep_i The log tranform of the effector dissociation constant from
    *        the active repressor, Ki, in kBT.
    * @param ep_ai The energetic difference between the active and inactive
    *        states of the repressor in kBT.
    * @param n_sites The number of allostericaly dependent effector binding
    *        sites.
    **/
    real fold_change(real R, real Nns, real ep_r, real c, real ep_a, real ep_i,
                    real ep_ai, int n_sites) {
        // Compute the various componenets piecewise for simplicity.
        real pact;
        real rep;
        pact = prob_act(c, ep_a, ep_i, ep_ai, n_sites);
        rep = repression(pact, R, Nns, ep_r);
        return rep^-1;
        }
      }

data {
    int<lower=0> N; // Number of data points
    vector[N] fc; // Experimentally measured fold-change
    vector[N] c; // Inducer concentration
    real<lower=0> R; //Number of repressors per cell.
    real<lower=0> Nns;  // Number of nonspecific binding sites for repressor.
    real ep_ai; // Allosteric energy difference.
    real ep_r;  // DNA binding energy of repressor.
    int<lower=1> n_sites; // Number of allosterically dependent sites.
    }

parameters {
    real ep_a; // Log transform of Ka.
    real ep_i; // Log transform of Ki.
    real<lower=0> sigma; // Constant measurement error.
    }

model {
    // Compute the expected fold-change. Stan does not allow for elementwise
    // operations on row vectors, so this is split into a simple loop.
    vector[N] fc_theo;
    for (k in 1:N) {
        fc_theo[k] = fold_change(R, Nns, ep_r, c[k], ep_a, ep_i, ep_ai, n_sites);
        }

    // Define the priors for the binding constants.
    ep_a ~ normal(0, 10);
    ep_i ~ normal(0, 10);
    sigma ~ beta(0.5, 0.5);

    // Define the likelihood.
    fc ~ normal(fc_theo, sigma);
  }
