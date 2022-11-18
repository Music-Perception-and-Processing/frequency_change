READ ME 

This is to explain the most important scripts and routines accompanying the data analysis and modeling part of the manuscript "A unitary model of auditory frequency change perception" by Kai Siedenburg, Jackson Graves, and Daniel Pressnitzer. Should run with Matlab 2020b or higher. 

Contact: kai.siedenburg [] uol.de

EXAMPLES

Contains a simple and isolated example to see how the basic synthesis and model is working. 

- minimal example for stimulus generation and modeling: shepard_example.m 


SCRIPTS

Contains main scripts to analyse data. 

- load participant data into structures: shepard_load_data.m 

- get info about participant bio data: shepard_participants_info.m

- evaluation of participant data: shepard_exps_eval.m

- run model simulation and save data: shepard_model_simulate.m

- evaluate model simulation from saved simulation data: shepard_model_simu_eval.m

- plot modeling results from from saved evaluation data: shepard_model_simu_plot.m



DATA_RAW

Contains raw data from experiments. 

- results files from experiments 1-3 with human participants 


FUNCTIONS

Contains most important synthesis and modeling functions. 

- compute shepard stimuli : shepard_spectrum_2d

- compute model output: shepard_model_2d

- plot weights in triangle space: triangle_plot3

- boostrap CIs: boot_CI

- fade signals appropriately: cos_ramp

- matrix helper: make_design_matrix


DATA 

Contains aggregated data from experiments with human participants or from simulations.

- index sets from simulation: simu_indmat_v8c.mat

- simulated response behavior of model: simu_resp_mat_v8c.mat

- bootstrap simulations: res_mat_revision.mat


RECREATION OF FIGURES 

Because Github does not allow to upload files of the required size, we cannot provide the folder DATA that is required to directly reproducing all figures and results [without running all simulations first]. If you want to get direct access to these data, please contact kai.siedenburg [] uol.de. 


- To recreate Figs 2-4 and Fig S1, run shepard_exps_eval.m

- To recreate Figs 5-8 and Figs. S2-S5, run shepard_model_simu_plot.m

