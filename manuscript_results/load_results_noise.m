clear 
%close all

%% Load data
%cd('Results');

K_vec = [1 2 5 10 20 50 100];
G_vec = [1 2 5 10 20 50 100];
K_str = {'001', '002', '005', '010', '020', '050', '100'};
G_str = {'001', '002', '005', '010', '020', '050', '100'};
VarList = getParamList(K_vec, G_vec, K_str, G_str);
NumSimuls = 10;

% Each data structure contains results for each K and G combination. 
% For each K and G combination, results are stored as:
%   RefSt_fwd: forward problem results for unstrained reference state
%   DefSt_fwd: forward problem results for strained drop
%   RefSt_inv: inverse problem results for unstrained reference state
%   DefSt_inv: inverse problem results for strained drop
%   Params_phys: physical parameters of the simulation
%   Params_num: numerical parameters of the simulation

% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90_Noiem3 = assembleResults_Noise(0.1, 1, 60, 'pepicelli', 0.9, K_vec, G_vec, 1e-3, 8, NumSimuls, VarList);
Wo01_Ar5_Sig60_Pep90_Noiem3 = assembleResults_Noise(0.1, 5, 60, 'pepicelli', 0.9, K_vec, G_vec, 1e-3, 8, NumSimuls, VarList);
Wo1_Ar1_Sig60_Pep90_Noiem3  = assembleResults_Noise(1, 1, 60, 'pepicelli', 0.9, K_vec, G_vec, 1e-3, 8, NumSimuls, VarList);
Wo1_Ar5_Sig60_Pep90_Noiem3  = assembleResults_Noise(1, 5, 60, 'pepicelli', 0.9, K_vec, G_vec, 1e-3, 8, NumSimuls, VarList);
return
%% Fix Sigma and compute average area reference conf.
% Because I saved the average and standard deviations for the reference
% state as the deformed state in example_elastic_inverse_noise.m, so I need
% to recompute the true average and standard deviation for sigma at the
% drop apex for the deformed configuration.

Wo01_Ar1_Sig60_Pep90_Noiem3 = ComputeRelevantQuantities_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, VarList, NumSimuls);
Wo01_Ar5_Sig60_Pep90_Noiem3 = ComputeRelevantQuantities_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, VarList, NumSimuls);
Wo1_Ar1_Sig60_Pep90_Noiem3 = ComputeRelevantQuantities_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, VarList, NumSimuls);
Wo1_Ar5_Sig60_Pep90_Noiem3 = ComputeRelevantQuantities_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, VarList, NumSimuls);

%% Fit circle for Capillary Pressure Microtensiometry

fraction = 0.75;

% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90_Noiem3 = calculateStressCircle_Ref_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, VarList, NumSimuls, fraction, 0);
Wo01_Ar5_Sig60_Pep90_Noiem3 = calculateStressCircle_Ref_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, VarList, NumSimuls, fraction, 0);
Wo1_Ar1_Sig60_Pep90_Noiem3 = calculateStressCircle_Ref_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, VarList, NumSimuls, fraction, 0);
Wo1_Ar5_Sig60_Pep90_Noiem3 = calculateStressCircle_Ref_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, VarList, NumSimuls, fraction, 0);

Wo01_Ar1_Sig60_Pep90_Noiem3 = calculateStressCircle_Def_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, VarList, NumSimuls, fraction, 0);
Wo01_Ar5_Sig60_Pep90_Noiem3 = calculateStressCircle_Def_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, VarList, NumSimuls, fraction, 0);
Wo1_Ar1_Sig60_Pep90_Noiem3 = calculateStressCircle_Def_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, VarList, NumSimuls, fraction, 0);
Wo1_Ar5_Sig60_Pep90_Noiem3 = calculateStressCircle_Def_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, VarList, NumSimuls, fraction, 0);

Wo01_Ar1_Sig60_Pep90_Noiem3 = calculateModulusCircle_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo01_Ar5_Sig60_Pep90_Noiem3 = calculateModulusCircle_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo1_Ar1_Sig60_Pep90_Noiem3 = calculateModulusCircle_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo1_Ar5_Sig60_Pep90_Noiem3 = calculateModulusCircle_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, NumSimuls, VarList);

%% Make tables with analyzed data

% Pepicelli, sigma = 60, strain = 0.9
makeTable_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, NumSimuls, VarList);
makeTable_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, NumSimuls, VarList);
makeTable_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, NumSimuls, VarList);
makeTable_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, NumSimuls, VarList);

%% Plot apical strain parameters for inverse problem
x = [1 2 5 10 20 50 100]; y = x;
[X,Y] = meshgrid(x,y);

% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90_Noiem3_ApexStrainTable = getStrainApexNonhomogen_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo01_Ar5_Sig60_Pep90_Noiem3_ApexStrainTable = getStrainApexNonhomogen_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo1_Ar1_Sig60_Pep90_Noiem3_ApexStrainTable = getStrainApexNonhomogen_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo1_Ar5_Sig60_Pep90_Noiem3_ApexStrainTable = getStrainApexNonhomogen_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, NumSimuls, VarList);

% plot different Ar and Wo combinations, for a given strain 
% Pepicelli, sigma = 60, strain = 0.9
figure
getStrainApexNonhomogenPlot_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, VarList, Wo01_Ar1_Sig60_Pep90_Noiem3_ApexStrainTable, X, Y, 1);
getStrainApexNonhomogenPlot_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, VarList, Wo01_Ar5_Sig60_Pep90_Noiem3_ApexStrainTable, X, Y, 2);
getStrainApexNonhomogenPlot_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, VarList, Wo1_Ar1_Sig60_Pep90_Noiem3_ApexStrainTable, X, Y, 3);
getStrainApexNonhomogenPlot_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, VarList, Wo1_Ar5_Sig60_Pep90_Noiem3_ApexStrainTable, X, Y, 4);
figure
getStrainApexNonhomogenPlotSTD_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, VarList, Wo01_Ar1_Sig60_Pep90_Noiem3_ApexStrainTable, X, Y, 1);
getStrainApexNonhomogenPlotSTD_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, VarList, Wo01_Ar5_Sig60_Pep90_Noiem3_ApexStrainTable, X, Y, 2);
getStrainApexNonhomogenPlotSTD_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, VarList, Wo1_Ar1_Sig60_Pep90_Noiem3_ApexStrainTable, X, Y, 3);
getStrainApexNonhomogenPlotSTD_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, VarList, Wo1_Ar5_Sig60_Pep90_Noiem3_ApexStrainTable, X, Y, 4);

%% Plot apical stress parameters for inverse problem

% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90_Noiem3_ApexStressTable = getStressApexNonhomogen_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo01_Ar5_Sig60_Pep90_Noiem3_ApexStressTable = getStressApexNonhomogen_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo1_Ar1_Sig60_Pep90_Noiem3_ApexStressTable = getStressApexNonhomogen_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo1_Ar5_Sig60_Pep90_Noiem3_ApexStressTable = getStressApexNonhomogen_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, NumSimuls, VarList);

% plot different Ar and Wo combinations, for a given strain 
% Pepicelli, sigma = 60, strain = 0.9
figure
getStressApexNonhomogenPlot_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, VarList, Wo01_Ar1_Sig60_Pep90_Noiem3_ApexStressTable, X, Y, 1);
getStressApexNonhomogenPlot_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, VarList, Wo01_Ar5_Sig60_Pep90_Noiem3_ApexStressTable, X, Y, 2);
getStressApexNonhomogenPlot_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, VarList, Wo1_Ar1_Sig60_Pep90_Noiem3_ApexStressTable, X, Y, 3);
getStressApexNonhomogenPlot_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, VarList, Wo1_Ar5_Sig60_Pep90_Noiem3_ApexStressTable, X, Y, 4);
figure
getStressApexNonhomogenPlotSTD_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, VarList, Wo01_Ar1_Sig60_Pep90_Noiem3_ApexStressTable, X, Y, 1);
getStressApexNonhomogenPlotSTD_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, VarList, Wo01_Ar5_Sig60_Pep90_Noiem3_ApexStressTable, X, Y, 2);
getStressApexNonhomogenPlotSTD_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, VarList, Wo1_Ar1_Sig60_Pep90_Noiem3_ApexStressTable, X, Y, 3);
getStressApexNonhomogenPlotSTD_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, VarList, Wo1_Ar5_Sig60_Pep90_Noiem3_ApexStressTable, X, Y, 4);

%% Plot modulus parameters

% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90_Noiem3_ModulusTable = getModulusNonhomogen_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo01_Ar5_Sig60_Pep90_Noiem3_ModulusTable = getModulusNonhomogen_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo1_Ar1_Sig60_Pep90_Noiem3_ModulusTable = getModulusNonhomogen_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, NumSimuls, VarList);
Wo1_Ar5_Sig60_Pep90_Noiem3_ModulusTable = getModulusNonhomogen_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, NumSimuls, VarList);

% plot different Ar and Wo combinations, for a given strain 
% Pepicelli, sigma = 60, strain = 0.9
figure
getModulusNonhomogenPlot_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, VarList, Wo01_Ar1_Sig60_Pep90_Noiem3_ModulusTable, X, Y, 1);
getModulusNonhomogenPlot_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, VarList, Wo01_Ar5_Sig60_Pep90_Noiem3_ModulusTable, X, Y, 2);
getModulusNonhomogenPlot_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, VarList, Wo1_Ar1_Sig60_Pep90_Noiem3_ModulusTable, X, Y, 3);
getModulusNonhomogenPlot_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, VarList, Wo1_Ar5_Sig60_Pep90_Noiem3_ModulusTable, X, Y, 4);
figure
getModulusNonhomogenPlotSTD_Noise(Wo01_Ar1_Sig60_Pep90_Noiem3, VarList, Wo01_Ar1_Sig60_Pep90_Noiem3_ModulusTable, X, Y, 1);
getModulusNonhomogenPlotSTD_Noise(Wo01_Ar5_Sig60_Pep90_Noiem3, VarList, Wo01_Ar5_Sig60_Pep90_Noiem3_ModulusTable, X, Y, 2);
getModulusNonhomogenPlotSTD_Noise(Wo1_Ar1_Sig60_Pep90_Noiem3, VarList, Wo1_Ar1_Sig60_Pep90_Noiem3_ModulusTable, X, Y, 3);
getModulusNonhomogenPlotSTD_Noise(Wo1_Ar5_Sig60_Pep90_Noiem3, VarList, Wo1_Ar5_Sig60_Pep90_Noiem3_ModulusTable, X, Y, 4);

