% get it all together 

%% LOAD DATA
cd(fullfile('..', 'data_raw'))

%% Exp 1
cd('Data_Exp1')
% cd('/Users/kaisiedenburg/Dropbox (Personal)/2dShepard/ExperimentScriptsData/Data_Exp1')
filens = {
        'shep2d_Exp1a_part_1_2019_2_19_14_21.mat';
        'shep2d_Exp1a_2_part_1_2019_2_19_15_18.mat';
        'shep2d_Exp1a_1_part_2_2019_2_20_12_47.mat';
        'shep2d_Exp1a_2_part_2_2019_2_20_13_44.mat';
        'shep2d_Exp1a_1_part_3_2019_2_20_15_57.mat';
        'shep2d_Exp1a_2_part_3_2019_2_20_17_6.mat';
        'shep2d_Exp1a_1_part_4_2019_2_21_13_14.mat';
        'shep2d_Exp1a_2_part_4_2019_2_21_14_6.mat'
        'shep2d_Exp1a_1_part_5_2019_3_5_8_16.mat'
        'shep2d_Exp1a_2_part_5_2019_3_5_9_19.mat'
        'shep2d_Exp1a_1_part_6_2019_3_6_16_10.mat'
        'shep2d_Exp1a_2_part_6_2019_3_6_17_7.mat'
        'shep2d_Exp1a_1_part_7_2019_3_8_8_9.mat'
        'shep2d_Exp1a_2_part_7_2019_3_8_9_5.mat'
        'shep2d_Exp1a_1_part_8_2019_3_8_11_39.mat'
        'shep2d_Exp1a_2_part_8_2019_3_8_12_39.mat'
        'shep2d_Exp1a_1_part_9_2019_3_8_16_17.mat'
        'shep2d_Exp1a_2_part_9_2019_3_8_17_28.mat'
        'shep2d_Exp1a_1_part_10_2019_3_11_16_49.mat'
        'shep2d_Exp1a_2_part_10_2019_3_11_17_41.mat'
        'shep2d_Exp1a_1_part_11_2019_3_13_13_23.mat'
        'shep2d_Exp1a_2_part_11_2019_3_13_14_21.mat'
        'shep2d_Exp1a_1_part_12_2019_3_13_15_39.mat'
        'shep2d_Exp1a_2_part_12_2019_3_13_16_45.mat'
        'shep2d_Exp1a_1_part_13_2019_3_14_13_45.mat'
        'shep2d_Exp1a_2_part_13_2019_3_14_14_43.mat'
        };
    
resp_mat.e1 = [];    
for k = 1:length(filens)
    xx = load(filens{k});
    resp_mat.e1 = [resp_mat.e1; xx.resp_mat];
end

%% Exp 2
cd(fullfile('..','Data_Exp2'))
cd('/Users/kaisiedenburg/Dropbox (Personal)/2dShepard/ExperimentScriptsData/Data_Exp2')
filens = {
        'shep2d_Exp2_part_1_2019_4_16_10_49.mat'
        'shep2d_Exp2_part_2_2019_4_16_14_22.mat'
        'shep2d_Exp2_part_3_2019_4_26_15_3.mat'
        'shep2d_Exp2_part_4_2019_4_30_8_22.mat'
        'shep2d_Exp2_part_5_2019_4_30_9_54.mat'
        'shep2d_Exp2_part_6_2019_5_2_12_12.mat'
        'shep2d_Exp2_part_7_2019_5_2_14_31.mat'
        'shep2d_Exp2_part_8_2019_5_2_16_11.mat'
        'shep2d_Exp2_part_9_2019_5_7_8_20.mat'
        'shep2d_Exp2_part_10_2019_5_8_10_25.mat'
        'shep2d_Exp2_part_11_2019_5_8_15_22.mat'
        'shep2d_Exp2_part_12_2019_5_9_14_20.mat'
        };
    
resp_mat.e2 = [];
for k = 1:length(filens)
    xx = load(filens{k});
    resp_mat.e2 = [resp_mat.e2; xx.resp_mat];
end

%% Exp 3
cd(fullfile('..','Data_Exp3'))
parts = [1:12];
  
filens = {
    'shep2d_Exp3A_SFS_part_1_2019_6_3_10_20.mat'
    'shep2d_Exp3A_SFS_part_2_2019_6_7_10_41.mat'
    'shep2d_Exp3A_SFS_part_3_2019_6_10_11_54.mat'
    'shep2d_Exp3A_SFS_part_4_2019_6_11_14_27.mat'
    'shep2d_Exp3A_SFS_part_5_2019_6_11_16_1.mat'
    'shep2d_Exp3A_SFS_part_6_2019_6_13_10_46.mat'
    'shep2d_Exp3A_SFS_part_7_2019_6_13_15_13.mat'
    'shep2d_Exp3A_SFS_part_8_2019_6_14_8_54.mat'
    'shep2d_Exp3A_SFS_part_9_2019_6_17_8_32.mat'
    'shep2d_Exp3A_SFS_part_10_2019_6_19_8_23.mat'
    'shep2d_Exp3A_SFS_part_11_2019_6_19_16_18.mat'
    'shep2d_Exp3A_SFS_part_12_2019_6_22_11_1.mat'
    'shep2d_Exp3A_ENV_part_1_2019_6_3_10_26.mat'
    'shep2d_Exp3A_ENV_part_2_2019_6_7_10_35.mat'
    'shep2d_Exp3A_ENV_part_3_2019_6_10_11_59.mat'
    'shep2d_Exp3A_ENV_part_4_2019_6_11_14_21.mat'
    'shep2d_Exp3A_ENV_part_5_2019_6_11_16_6.mat'
    'shep2d_Exp3A_ENV_part_6_2019_6_13_10_40.mat'
    'shep2d_Exp3A_ENV_part_7_2019_6_13_15_19.mat'
    'shep2d_Exp3A_ENV_part_8_2019_6_14_8_48.mat'
    'shep2d_Exp3A_ENV_part_9_2019_6_17_8_38.mat'
    'shep2d_Exp3A_ENV_part_10_2019_6_19_8_16.mat'
    'shep2d_Exp3A_ENV_part_11_2019_6_19_16_24.mat'
    'shep2d_Exp3A_ENV_part_12_2019_6_22_10_56.mat'
    };
       
resp_mat.e3a = [];
for k = 1:length(filens)
    xx = load(filens{k});
    resp_mat.e3a = [resp_mat.e3a; xx.resp_mat];
end


% SFS-ENV combined 
filens = {
  'shep2d_Exp3B_part_1_2019_6_3_10_32.mat'
  'shep2d_Exp3B_part_2_2019_6_7_10_47.mat'
  'shep2d_Exp3B_part_3_2019_6_10_12_15.mat'
  'shep2d_Exp3B_part_4_2019_6_11_14_34.mat'
  'shep2d_Exp3B_part_5_2019_6_11_16_12.mat'
  'shep2d_Exp3B_part_6_2019_6_13_10_53.mat'
  'shep2d_Exp3B_part_7_2019_6_13_15_29.mat'
  'shep2d_Exp3B_part_8_2019_6_14_9_1.mat'
  'shep2d_Exp3B_part_9_2019_6_17_8_48.mat'
  'shep2d_Exp3B_part_10_2019_6_19_8_32.mat'
  'shep2d_Exp3B_part_11_2019_6_19_16_32.mat'
  'shep2d_Exp3B_part_12_2019_6_22_11_9.mat'
    };
       
resp_mat.e3b = [];
for k = 1:length(filens)
    xx = load(filens{k});
    resp_mat.e3b = [resp_mat.e3b; xx.resp_mat];
end
