
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
    // Dimensional quantities
    int<lower=1> N; // Total Number of measurements.
    int<lower=1> NR; // Number of unique RBS sequences.
    int<lower=1> Nep;// Number of unique operator sequences.
    int<lower=1, upper=Nep> ep_idx[N];
    int<lower=1, upper=NR> R_idx[N];

    // Parametric quantities
    vector<lower=1>[NR] R_mu; // Experimentally determined mean repressors
    vector<lower=0>[NR] R_sig; // Experimental error on repressors

    // Constants
    real ep_ai; // Energy difference between active/inactive states
    int<lower=1> n_sites; // Number of allosterically coupled sites
    real<lower=0> n_ns; // Number of nonspecific binding sites

    // Observed quantities
    vector[N] fc_obs; // Observed fold-change 
    vector<lower=0>[N] c; // IPTG concentration
}

transformed data {
    vector[NR] logR_mu = log(R_mu);
    vector[NR] logR_sig = log(R_sig);
}

parameters {
    vector<lower=0>[NR] log_R; // Number of repressors
    vector[Nep] ep_r; // DNA binding affinity
    real ep_a; // Active repressor inducer dissociation constant
    real ep_i; // Inactive repressor inducer dissociation constant
    real<lower=0> sigma; //Homoscedastic error term
}

transformed parameters {
    vector[NR] R = exp(log_R);
}

model {
    vector[N] mu;

    // Priors
    log_R ~ normal(logR_mu, logR_sig);
    ep_r ~ normal(-12, 6);
    ep_a ~ normal(0, 10);
    ep_i ~ normal(0, 10);
    sigma ~ normal(0, 0.1);

    // Likelihood
    for (i in 1:N) {
    mu[i] = fold_change(R[R_idx[i]], n_ns, ep_r[ep_idx[i]], c[i], 
                    ep_a, ep_i, ep_ai, n_sites);

    }
    fc_obs ~ normal(mu, sigma); 
}

generated quantities {
    real ka = exp(ep_a);
    real ki = exp(ep_i);
}