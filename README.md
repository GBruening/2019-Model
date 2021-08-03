# Data analysis scripts for paper

License - CC0 Universal

Resample trajectory data used to simulate movements is stored in Resamp_data_Nov2020.
 - This is aggregated data from a previous study, code can be found [here](https://github.com/GBruening/Met_mass_r).
 
 The main file to compute the simulation is supercompute_arm_model.m.
 This will output multiple .mat files, which are used in Create_Plots.m to generate the .csv files used in the main analysis.
 met_model_stats.R will compute all stats and generate figures used in the paper.
