

data { 
    // Dimensional parameters
    int<lower=1> J; // Number of unique replicates
    int<lower=1> N_fluct; // Number of fluctuation measurements
    int<lower=1> N_foldchange; // Number of dilution strain fold-chang measurements
    int<lower=1> N_delta; // Number of constitituve expression measurements

    // Identification vectors
    int<lower=1, upper=J> fluct_idx[N_fluct] ; 
    int<lower=1, upper=J> foldchange_idx[N_foldchange]; 
    int<lower=1, upper=J> delta_idx[N_delta];

    // Observed data
    vector<lower=0>[N_fluct] I1;
    vector[N_delta] delta_yfp;
    vector[N_foldchange] foldchange_mcherry;
    vector[N_foldchange] foldchange_yfp;



}