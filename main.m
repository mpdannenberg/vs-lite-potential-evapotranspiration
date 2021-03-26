% Main file

% Data/Methods figures
make_itrdb_sites_map;
make_vslite_conceptual_figure;

% Calibrate VS-Lite parameters (takes a long time... a few hours)
calibrate_itrdb_parameters;

% Calculate basic validation stats for VS-Lite model over independent
% validation period 
calculate_vslite_stats;

% Make validation maps
make_vslite_r2_maps;
make_vslite_r2_histograms;

% Get downscaled CMIP5 projections and simulate change in PET, soil
% moisture, moisture stress, and growth
get_bcsd_projections;
calculate_vslite_bcsd_projections;
plot_vslite_bcsd_projections; 
supplemental_plot_vslite_bcsd_projections_rcp45; 
plot_vslite_bcsd_gM_trajectory;

