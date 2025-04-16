clear 
%close all

%% Load data
%cd('Results');


K_vec = [1 2 5 10 20 50 100];
G_vec = [1 2 5 10 20 50 100];
K_str = {'001', '002', '005', '010', '020', '050', '100'};
G_str = {'001', '002', '005', '010', '020', '050', '100'};
VarList = getParamList(K_vec, G_vec, K_str, G_str);


% Each data structure contains results for each K and G combination. 
% For each K and G combination, results are stored as:
%   RefSt_fwd: forward problem results for unstrained reference state
%   DefSt_fwd: forward problem results for strained drop
%   RefSt_inv: inverse problem results for unstrained reference state
%   DefSt_inv: inverse problem results for strained drop
%   Params_phys: physical parameters of the simulation
%   Params_num: numerical parameters of the simulation

% Pepicelli, sigma = 20, strain = 0.8
Wo01_Ar1_Sig20_Pep80 = assembleResults(0.1, 1, 20, 'pepicelli', 0.8, K_vec, G_vec, VarList);
Wo01_Ar5_Sig20_Pep80 = assembleResults(0.1, 5, 20, 'pepicelli', 0.8, K_vec, G_vec, VarList);
Wo1_Ar1_Sig20_Pep80  = assembleResults(1, 1, 20, 'pepicelli', 0.8, K_vec, G_vec, VarList);
Wo1_Ar5_Sig20_Pep80  = assembleResults(1, 5, 20, 'pepicelli', 0.8, K_vec, G_vec, VarList);
% Pepicelli, sigma = 20, strain = 0.9
Wo01_Ar1_Sig20_Pep90 = assembleResults(0.1, 1, 20, 'pepicelli', 0.9, K_vec, G_vec, VarList);
Wo01_Ar5_Sig20_Pep90 = assembleResults(0.1, 5, 20, 'pepicelli', 0.9, K_vec, G_vec, VarList);
Wo1_Ar1_Sig20_Pep90  = assembleResults(1, 1, 20, 'pepicelli', 0.9, K_vec, G_vec, VarList);
Wo1_Ar5_Sig20_Pep90  = assembleResults(1, 5, 20, 'pepicelli', 0.9, K_vec, G_vec, VarList);
% Pepicelli, sigma = 20, strain = 0.95
Wo01_Ar1_Sig20_Pep95 = assembleResults(0.1, 1, 20, 'pepicelli', 0.95, K_vec, G_vec, VarList);
Wo01_Ar5_Sig20_Pep95 = assembleResults(0.1, 5, 20, 'pepicelli', 0.95, K_vec, G_vec, VarList);
Wo1_Ar1_Sig20_Pep95  = assembleResults(1, 1, 20, 'pepicelli', 0.95, K_vec, G_vec, VarList);
Wo1_Ar5_Sig20_Pep95  = assembleResults(1, 5, 20, 'pepicelli', 0.95, K_vec, G_vec, VarList);
% Pepicelli, sigma = 60, strain = 0.8
Wo01_Ar1_Sig60_Pep80 = assembleResults(0.1, 1, 60, 'pepicelli', 0.8, K_vec, G_vec, VarList);
Wo01_Ar5_Sig60_Pep80 = assembleResults(0.1, 5, 60, 'pepicelli', 0.8, K_vec, G_vec, VarList);
Wo1_Ar1_Sig60_Pep80  = assembleResults(1, 1, 60, 'pepicelli', 0.8, K_vec, G_vec, VarList);
Wo1_Ar5_Sig60_Pep80  = assembleResults(1, 5, 60, 'pepicelli', 0.8, K_vec, G_vec, VarList);
% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90 = assembleResults(0.1, 1, 60, 'pepicelli', 0.9, K_vec, G_vec, VarList);
Wo01_Ar5_Sig60_Pep90 = assembleResults(0.1, 5, 60, 'pepicelli', 0.9, K_vec, G_vec, VarList);
Wo1_Ar1_Sig60_Pep90  = assembleResults(1, 1, 60, 'pepicelli', 0.9, K_vec, G_vec, VarList);
Wo1_Ar5_Sig60_Pep90  = assembleResults(1, 5, 60, 'pepicelli', 0.9, K_vec, G_vec, VarList);
% Pepicelli, sigma = 60, strain = 0.95
Wo01_Ar1_Sig60_Pep95 = assembleResults(0.1, 1, 60, 'pepicelli', 0.95, K_vec, G_vec, VarList);
Wo01_Ar5_Sig60_Pep95 = assembleResults(0.1, 5, 60, 'pepicelli', 0.95, K_vec, G_vec, VarList);
Wo1_Ar1_Sig60_Pep95  = assembleResults(1, 1, 60, 'pepicelli', 0.95, K_vec, G_vec, VarList);
Wo1_Ar5_Sig60_Pep95  = assembleResults(1, 5, 60, 'pepicelli', 0.95, K_vec, G_vec, VarList);
return

%% Plot Shapes example

figure
pl_Wo01_Ar1_Sig60_Ref = plotShapes(Wo01_Ar1_Sig60_Pep80, 1, 'Ref');
pl_Wo01_Ar1_Sig60_Pep95 = plotShapes(Wo01_Ar1_Sig60_Pep95, 1, 'Def');
pl_Wo01_Ar1_Sig60_Pep90 = plotShapes(Wo01_Ar1_Sig60_Pep90, 1, 'Def');
pl_Wo01_Ar1_Sig60_Pep80 = plotShapes(Wo01_Ar1_Sig60_Pep80, 1, 'Def');
legend([pl_Wo01_Ar1_Sig60_Ref, pl_Wo01_Ar1_Sig60_Pep95, pl_Wo01_Ar1_Sig60_Pep90, pl_Wo01_Ar1_Sig60_Pep80],...
    {'Ref.', '$\Delta A / A_0 = 5$ \%', '$\Delta A / A_0 = 10$ \%', '$\Delta A / A_0 = 20$ \%'}, ...
    'Location', 'SouthEast', 'Interpreter', 'latex');
figure
pl_Wo01_Ar5_Sig60_Ref = plotShapes(Wo01_Ar5_Sig60_Pep80, 2, 'Ref');
pl_Wo01_Ar5_Sig60_Pep95 = plotShapes(Wo01_Ar5_Sig60_Pep95, 2, 'Def');
pl_Wo01_Ar5_Sig60_Pep90 = plotShapes(Wo01_Ar5_Sig60_Pep90, 2, 'Def');
pl_Wo01_Ar5_Sig60_Pep80 = plotShapes(Wo01_Ar5_Sig60_Pep80, 2, 'Def');
figure
pl_Wo1_Ar1_Sig60_Ref = plotShapes(Wo1_Ar1_Sig60_Pep80, 3, 'Ref');
pl_Wo1_Ar1_Sig60_Pep95 = plotShapes(Wo1_Ar1_Sig60_Pep95, 3, 'Def');
pl_Wo1_Ar1_Sig60_Pep90 = plotShapes(Wo1_Ar1_Sig60_Pep90, 3, 'Def');
pl_Wo1_Ar1_Sig60_Pep80 = plotShapes(Wo1_Ar1_Sig60_Pep80, 3, 'Def');
legend([pl_Wo1_Ar1_Sig60_Ref, pl_Wo1_Ar1_Sig60_Pep95, pl_Wo1_Ar1_Sig60_Pep90, pl_Wo1_Ar1_Sig60_Pep80],...
    {'Ref.', '$\Delta A / A_0 = 5$ \%', '$\Delta A / A_0 = 10$ \%', '$\Delta A / A_0 = 20$ \%'}, ...
    'Location', 'SouthEast', 'Interpreter', 'latex');
figure
pl_Wo1_Ar5_Sig60_Ref = plotShapes(Wo1_Ar5_Sig60_Pep80, 4, 'Ref');
pl_Wo1_Ar5_Sig60_Pep95 = plotShapes(Wo1_Ar5_Sig60_Pep95, 4, 'Def');
pl_Wo1_Ar5_Sig60_Pep90 = plotShapes(Wo1_Ar5_Sig60_Pep90, 4, 'Def');
pl_Wo1_Ar5_Sig60_Pep80 = plotShapes(Wo1_Ar5_Sig60_Pep80, 4, 'Def');

%% Find dilatational nonhomogeneity
x = [1 2 5 10 20 50 100]; y = x;
[X,Y] = meshgrid(x,y);

% Pepicelli, sigma = 20, strain = 0.8
Wo01_Ar1_Sig20_Pep80_DilParamTable = getDilParam(Wo01_Ar1_Sig20_Pep80, VarList);
Wo01_Ar5_Sig20_Pep80_DilParamTable = getDilParam(Wo01_Ar5_Sig20_Pep80, VarList);
Wo1_Ar1_Sig20_Pep80_DilParamTable = getDilParam(Wo1_Ar1_Sig20_Pep80, VarList);
Wo1_Ar5_Sig20_Pep80_DilParamTable = getDilParam(Wo1_Ar5_Sig20_Pep80, VarList);
% Pepicelli, sigma = 20, strain = 0.9
Wo01_Ar1_Sig20_Pep90_DilParamTable = getDilParam(Wo01_Ar1_Sig20_Pep90, VarList);
Wo01_Ar5_Sig20_Pep90_DilParamTable = getDilParam(Wo01_Ar5_Sig20_Pep90, VarList);
Wo1_Ar1_Sig20_Pep90_DilParamTable = getDilParam(Wo1_Ar1_Sig20_Pep90, VarList);
Wo1_Ar5_Sig20_Pep90_DilParamTable = getDilParam(Wo1_Ar5_Sig20_Pep90, VarList);
% Pepicelli, sigma = 20, strain = 0.95
Wo01_Ar1_Sig20_Pep95_DilParamTable = getDilParam(Wo01_Ar1_Sig20_Pep95, VarList);
Wo01_Ar5_Sig20_Pep95_DilParamTable = getDilParam(Wo01_Ar5_Sig20_Pep95, VarList);
Wo1_Ar1_Sig20_Pep95_DilParamTable = getDilParam(Wo1_Ar1_Sig20_Pep95, VarList);
Wo1_Ar5_Sig20_Pep95_DilParamTable = getDilParam(Wo1_Ar5_Sig20_Pep95, VarList);
% Pepicelli, sigma = 60, strain = 0.8
Wo01_Ar1_Sig60_Pep80_DilParamTable = getDilParam(Wo01_Ar1_Sig60_Pep80, VarList);
Wo01_Ar5_Sig60_Pep80_DilParamTable = getDilParam(Wo01_Ar5_Sig60_Pep80, VarList);
Wo1_Ar1_Sig60_Pep80_DilParamTable = getDilParam(Wo1_Ar1_Sig60_Pep80, VarList);
Wo1_Ar5_Sig60_Pep80_DilParamTable = getDilParam(Wo1_Ar5_Sig60_Pep80, VarList);
% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90_DilParamTable = getDilParam(Wo01_Ar1_Sig60_Pep90, VarList);
Wo01_Ar5_Sig60_Pep90_DilParamTable = getDilParam(Wo01_Ar5_Sig60_Pep90, VarList);
Wo1_Ar1_Sig60_Pep90_DilParamTable = getDilParam(Wo1_Ar1_Sig60_Pep90, VarList);
Wo1_Ar5_Sig60_Pep90_DilParamTable = getDilParam(Wo1_Ar5_Sig60_Pep90, VarList);
% Pepicelli, sigma = 60, strain = 0.95
Wo01_Ar1_Sig60_Pep95_DilParamTable = getDilParam(Wo01_Ar1_Sig60_Pep95, VarList);
Wo01_Ar5_Sig60_Pep95_DilParamTable = getDilParam(Wo01_Ar5_Sig60_Pep95, VarList);
Wo1_Ar1_Sig60_Pep95_DilParamTable = getDilParam(Wo1_Ar1_Sig60_Pep95, VarList);
Wo1_Ar5_Sig60_Pep95_DilParamTable = getDilParam(Wo1_Ar5_Sig60_Pep95, VarList);

% Pepicelli, sigma = 20, strain = 0.8
figure
getDilParamPlot(Wo01_Ar1_Sig20_Pep80, VarList, Wo01_Ar1_Sig20_Pep80_DilParamTable, X, Y, 1);
getDilParamPlot(Wo01_Ar5_Sig20_Pep80, VarList, Wo01_Ar5_Sig20_Pep80_DilParamTable, X, Y, 2);
getDilParamPlot(Wo1_Ar1_Sig20_Pep80, VarList, Wo1_Ar1_Sig20_Pep80_DilParamTable, X, Y, 3);
getDilParamPlot(Wo1_Ar5_Sig20_Pep80, VarList, Wo1_Ar5_Sig20_Pep80_DilParamTable, X, Y, 4);
% Pepicelli, sigma = 20, strain = 0.9
figure
getDilParamPlot(Wo01_Ar1_Sig20_Pep90, VarList, Wo01_Ar1_Sig20_Pep90_DilParamTable, X, Y, 1);
getDilParamPlot(Wo01_Ar5_Sig20_Pep90, VarList, Wo01_Ar5_Sig20_Pep90_DilParamTable, X, Y, 2);
getDilParamPlot(Wo1_Ar1_Sig20_Pep90, VarList, Wo1_Ar1_Sig20_Pep90_DilParamTable, X, Y, 3);
getDilParamPlot(Wo1_Ar5_Sig20_Pep90, VarList, Wo1_Ar5_Sig20_Pep90_DilParamTable, X, Y, 4);
% Pepicelli, sigma = 20, strain = 0.95
figure
getDilParamPlot(Wo01_Ar1_Sig20_Pep95, VarList, Wo01_Ar1_Sig20_Pep95_DilParamTable, X, Y, 1);
getDilParamPlot(Wo01_Ar5_Sig20_Pep95, VarList, Wo01_Ar5_Sig20_Pep95_DilParamTable, X, Y, 2);
getDilParamPlot(Wo1_Ar1_Sig20_Pep95, VarList, Wo1_Ar1_Sig20_Pep95_DilParamTable, X, Y, 3);
getDilParamPlot(Wo1_Ar5_Sig20_Pep95, VarList, Wo1_Ar5_Sig20_Pep95_DilParamTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.8
figure
getDilParamPlot(Wo01_Ar1_Sig60_Pep80, VarList, Wo01_Ar1_Sig60_Pep80_DilParamTable, X, Y, 1);
getDilParamPlot(Wo01_Ar5_Sig60_Pep80, VarList, Wo01_Ar5_Sig60_Pep80_DilParamTable, X, Y, 2);
getDilParamPlot(Wo1_Ar1_Sig60_Pep80, VarList, Wo1_Ar1_Sig60_Pep80_DilParamTable, X, Y, 3);
getDilParamPlot(Wo1_Ar5_Sig60_Pep80, VarList, Wo1_Ar5_Sig60_Pep80_DilParamTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.9
figure
getDilParamPlot(Wo01_Ar1_Sig60_Pep90, VarList, Wo01_Ar1_Sig60_Pep90_DilParamTable, X, Y, 1);
getDilParamPlot(Wo01_Ar5_Sig60_Pep90, VarList, Wo01_Ar5_Sig60_Pep90_DilParamTable, X, Y, 2);
getDilParamPlot(Wo1_Ar1_Sig60_Pep90, VarList, Wo1_Ar1_Sig60_Pep90_DilParamTable, X, Y, 3);
getDilParamPlot(Wo1_Ar5_Sig60_Pep90, VarList, Wo1_Ar5_Sig60_Pep90_DilParamTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.95
figure
getDilParamPlot(Wo01_Ar1_Sig60_Pep95, VarList, Wo01_Ar1_Sig60_Pep95_DilParamTable, X, Y, 1);
getDilParamPlot(Wo01_Ar5_Sig60_Pep95, VarList, Wo01_Ar5_Sig60_Pep95_DilParamTable, X, Y, 2);
getDilParamPlot(Wo1_Ar1_Sig60_Pep95, VarList, Wo1_Ar1_Sig60_Pep95_DilParamTable, X, Y, 3);
getDilParamPlot(Wo1_Ar5_Sig60_Pep95, VarList, Wo1_Ar5_Sig60_Pep95_DilParamTable, X, Y, 4);

%% Find stress nonhomogeneity

% Pepicelli, sigma = 20, strain = 0.8
Wo01_Ar1_Sig20_Pep80_StressParamTable = getStressParam(Wo01_Ar1_Sig20_Pep80, VarList);
Wo01_Ar5_Sig20_Pep80_StressParamTable = getStressParam(Wo01_Ar5_Sig20_Pep80, VarList);
Wo1_Ar1_Sig20_Pep80_StressParamTable = getStressParam(Wo1_Ar1_Sig20_Pep80, VarList);
Wo1_Ar5_Sig20_Pep80_StressParamTable = getStressParam(Wo1_Ar5_Sig20_Pep80, VarList);
% Pepicelli, sigma = 20, strain = 0.9
Wo01_Ar1_Sig20_Pep90_StressParamTable = getStressParam(Wo01_Ar1_Sig20_Pep90, VarList);
Wo01_Ar5_Sig20_Pep90_StressParamTable = getStressParam(Wo01_Ar5_Sig20_Pep90, VarList);
Wo1_Ar1_Sig20_Pep90_StressParamTable = getStressParam(Wo1_Ar1_Sig20_Pep90, VarList);
Wo1_Ar5_Sig20_Pep90_StressParamTable = getStressParam(Wo1_Ar5_Sig20_Pep90, VarList);
% Pepicelli, sigma = 20, strain = 0.95
Wo01_Ar1_Sig20_Pep95_StressParamTable = getStressParam(Wo01_Ar1_Sig20_Pep95, VarList);
Wo01_Ar5_Sig20_Pep95_StressParamTable = getStressParam(Wo01_Ar5_Sig20_Pep95, VarList);
Wo1_Ar1_Sig20_Pep95_StressParamTable = getStressParam(Wo1_Ar1_Sig20_Pep95, VarList);
Wo1_Ar5_Sig20_Pep95_StressParamTable = getStressParam(Wo1_Ar5_Sig20_Pep95, VarList);
% Pepicelli, sigma = 60, strain = 0.8
Wo01_Ar1_Sig60_Pep80_StressParamTable = getStressParam(Wo01_Ar1_Sig60_Pep80, VarList);
Wo01_Ar5_Sig60_Pep80_StressParamTable = getStressParam(Wo01_Ar5_Sig60_Pep80, VarList);
Wo1_Ar1_Sig60_Pep80_StressParamTable = getStressParam(Wo1_Ar1_Sig60_Pep80, VarList);
Wo1_Ar5_Sig60_Pep80_StressParamTable = getStressParam(Wo1_Ar5_Sig60_Pep80, VarList);
% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90_StressParamTable = getStressParam(Wo01_Ar1_Sig60_Pep90, VarList);
Wo01_Ar5_Sig60_Pep90_StressParamTable = getStressParam(Wo01_Ar5_Sig60_Pep90, VarList);
Wo1_Ar1_Sig60_Pep90_StressParamTable = getStressParam(Wo1_Ar1_Sig60_Pep90, VarList);
Wo1_Ar5_Sig60_Pep90_StressParamTable = getStressParam(Wo1_Ar5_Sig60_Pep90, VarList);
% Pepicelli, sigma = 60, strain = 0.95
Wo01_Ar1_Sig60_Pep95_StressParamTable = getStressParam(Wo01_Ar1_Sig60_Pep95, VarList);
Wo01_Ar5_Sig60_Pep95_StressParamTable = getStressParam(Wo01_Ar5_Sig60_Pep95, VarList);
Wo1_Ar1_Sig60_Pep95_StressParamTable = getStressParam(Wo1_Ar1_Sig60_Pep95, VarList);
Wo1_Ar5_Sig60_Pep95_StressParamTable = getStressParam(Wo1_Ar5_Sig60_Pep95, VarList);

% Pepicelli, sigma = 20, strain = 0.8
figure
getStressParamPlot(Wo01_Ar1_Sig20_Pep80, VarList, Wo01_Ar1_Sig20_Pep80_StressParamTable, X, Y, 1);
getStressParamPlot(Wo01_Ar5_Sig20_Pep80, VarList, Wo01_Ar5_Sig20_Pep80_StressParamTable, X, Y, 2);
getStressParamPlot(Wo1_Ar1_Sig20_Pep80, VarList, Wo1_Ar1_Sig20_Pep80_StressParamTable, X, Y, 3);
getStressParamPlot(Wo1_Ar5_Sig20_Pep80, VarList, Wo1_Ar5_Sig20_Pep80_StressParamTable, X, Y, 4);
% Pepicelli, sigma = 20, strain = 0.9
figure
getStressParamPlot(Wo01_Ar1_Sig20_Pep90, VarList, Wo01_Ar1_Sig20_Pep90_StressParamTable, X, Y, 1);
getStressParamPlot(Wo01_Ar5_Sig20_Pep90, VarList, Wo01_Ar5_Sig20_Pep90_StressParamTable, X, Y, 2);
getStressParamPlot(Wo1_Ar1_Sig20_Pep90, VarList, Wo1_Ar1_Sig20_Pep90_StressParamTable, X, Y, 3);
getStressParamPlot(Wo1_Ar5_Sig20_Pep90, VarList, Wo1_Ar5_Sig20_Pep90_StressParamTable, X, Y, 4);
% Pepicelli, sigma = 20, strain = 0.95
figure
getStressParamPlot(Wo01_Ar1_Sig20_Pep95, VarList, Wo01_Ar1_Sig20_Pep95_StressParamTable, X, Y, 1);
getStressParamPlot(Wo01_Ar5_Sig20_Pep95, VarList, Wo01_Ar5_Sig20_Pep95_StressParamTable, X, Y, 2);
getStressParamPlot(Wo1_Ar1_Sig20_Pep95, VarList, Wo1_Ar1_Sig20_Pep95_StressParamTable, X, Y, 3);
getStressParamPlot(Wo1_Ar5_Sig20_Pep95, VarList, Wo1_Ar5_Sig20_Pep95_StressParamTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.8
figure
getStressParamPlot(Wo01_Ar1_Sig60_Pep80, VarList, Wo01_Ar1_Sig60_Pep80_StressParamTable, X, Y, 1);
getStressParamPlot(Wo01_Ar5_Sig60_Pep80, VarList, Wo01_Ar5_Sig60_Pep80_StressParamTable, X, Y, 2);
getStressParamPlot(Wo1_Ar1_Sig60_Pep80, VarList, Wo1_Ar1_Sig60_Pep80_StressParamTable, X, Y, 3);
getStressParamPlot(Wo1_Ar5_Sig60_Pep80, VarList, Wo1_Ar5_Sig60_Pep80_StressParamTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.9
figure
getStressParamPlot(Wo01_Ar1_Sig60_Pep90, VarList, Wo01_Ar1_Sig60_Pep90_StressParamTable, X, Y, 1);
getStressParamPlot(Wo01_Ar5_Sig60_Pep90, VarList, Wo01_Ar5_Sig60_Pep90_StressParamTable, X, Y, 2);
getStressParamPlot(Wo1_Ar1_Sig60_Pep90, VarList, Wo1_Ar1_Sig60_Pep90_StressParamTable, X, Y, 3);
getStressParamPlot(Wo1_Ar5_Sig60_Pep90, VarList, Wo1_Ar5_Sig60_Pep90_StressParamTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.95
figure
getStressParamPlot(Wo01_Ar1_Sig60_Pep95, VarList, Wo01_Ar1_Sig60_Pep95_StressParamTable, X, Y, 1);
getStressParamPlot(Wo01_Ar5_Sig60_Pep95, VarList, Wo01_Ar5_Sig60_Pep95_StressParamTable, X, Y, 2);
getStressParamPlot(Wo1_Ar1_Sig60_Pep95, VarList, Wo1_Ar1_Sig60_Pep95_StressParamTable, X, Y, 3);
getStressParamPlot(Wo1_Ar5_Sig60_Pep95, VarList, Wo1_Ar5_Sig60_Pep95_StressParamTable, X, Y, 4);

%% Evaluate 3 Test cases
%  1. G = 50, K = 1
%  2. G = 1,  K = 50
%  3. G = 50, K = 50

% Make input data structures for all Wo and Ar combinations
TestCases = {'K001_G050', 'K050_G001', 'K050_G050'};
InputVarList = cell(1,3);
for i = 1:length(TestCases)
    InputVarList{1,i} = matlab.lang.makeValidName(TestCases{i});
end


% Pepicelli, sigma = 20, strain = 0.8
Wo01_Ar1_Sig20_Pep80_Comparison = getInputDataStruc(Wo01_Ar1_Sig20_Pep80, InputVarList);
Wo01_Ar5_Sig20_Pep80_Comparison = getInputDataStruc(Wo01_Ar5_Sig20_Pep80, InputVarList);
Wo1_Ar1_Sig20_Pep80_Comparison = getInputDataStruc(Wo1_Ar1_Sig20_Pep80, InputVarList);
Wo1_Ar5_Sig20_Pep80_Comparison = getInputDataStruc(Wo1_Ar5_Sig20_Pep80, InputVarList);
% Pepicelli, sigma = 20, strain = 0.9
Wo01_Ar1_Sig20_Pep90_Comparison = getInputDataStruc(Wo01_Ar1_Sig20_Pep90, InputVarList);
Wo01_Ar5_Sig20_Pep90_Comparison = getInputDataStruc(Wo01_Ar5_Sig20_Pep90, InputVarList);
Wo1_Ar1_Sig20_Pep90_Comparison = getInputDataStruc(Wo1_Ar1_Sig20_Pep90, InputVarList);
Wo1_Ar5_Sig20_Pep90_Comparison = getInputDataStruc(Wo1_Ar5_Sig20_Pep90, InputVarList);
% Pepicelli, sigma = 20, strain = 0.95
Wo01_Ar1_Sig20_Pep95_Comparison = getInputDataStruc(Wo01_Ar1_Sig20_Pep95, InputVarList);
Wo01_Ar5_Sig20_Pep95_Comparison = getInputDataStruc(Wo01_Ar5_Sig20_Pep95, InputVarList);
Wo1_Ar1_Sig20_Pep95_Comparison = getInputDataStruc(Wo1_Ar1_Sig20_Pep95, InputVarList);
Wo1_Ar5_Sig20_Pep95_Comparison = getInputDataStruc(Wo1_Ar5_Sig20_Pep95, InputVarList);
% Pepicelli, sigma = 60, strain = 0.8
Wo01_Ar1_Sig60_Pep80_Comparison = getInputDataStruc(Wo01_Ar1_Sig60_Pep80, InputVarList);
Wo01_Ar5_Sig60_Pep80_Comparison = getInputDataStruc(Wo01_Ar5_Sig60_Pep80, InputVarList);
Wo1_Ar1_Sig60_Pep80_Comparison = getInputDataStruc(Wo1_Ar1_Sig60_Pep80, InputVarList);
Wo1_Ar5_Sig60_Pep80_Comparison = getInputDataStruc(Wo1_Ar5_Sig60_Pep80, InputVarList);
% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90_Comparison = getInputDataStruc(Wo01_Ar1_Sig60_Pep90, InputVarList);
Wo01_Ar5_Sig60_Pep90_Comparison = getInputDataStruc(Wo01_Ar5_Sig60_Pep90, InputVarList);
Wo1_Ar1_Sig60_Pep90_Comparison = getInputDataStruc(Wo1_Ar1_Sig60_Pep90, InputVarList);
Wo1_Ar5_Sig60_Pep90_Comparison = getInputDataStruc(Wo1_Ar5_Sig60_Pep90, InputVarList);
% Pepicelli, sigma = 60, strain = 0.95
Wo01_Ar1_Sig60_Pep95_Comparison = getInputDataStruc(Wo01_Ar1_Sig60_Pep95, InputVarList);
Wo01_Ar5_Sig60_Pep95_Comparison = getInputDataStruc(Wo01_Ar5_Sig60_Pep95, InputVarList);
Wo1_Ar1_Sig60_Pep95_Comparison = getInputDataStruc(Wo1_Ar1_Sig60_Pep95, InputVarList);
Wo1_Ar5_Sig60_Pep95_Comparison = getInputDataStruc(Wo1_Ar5_Sig60_Pep95, InputVarList);


% Pepicelli, sigma = 20, strain = 0.8
figure
compareStrainStressProfilesPlot(Wo01_Ar1_Sig20_Pep80_Comparison, InputVarList, 1);
compareStrainStressProfilesPlot(Wo01_Ar5_Sig20_Pep80_Comparison, InputVarList, 2);
compareStrainStressProfilesPlot(Wo1_Ar1_Sig20_Pep80_Comparison,  InputVarList, 3);
compareStrainStressProfilesPlot(Wo1_Ar5_Sig20_Pep80_Comparison,  InputVarList, 4);
% Pepicelli, sigma = 20, strain = 0.9
figure
compareStrainStressProfilesPlot(Wo01_Ar1_Sig20_Pep90_Comparison, InputVarList, 1);
compareStrainStressProfilesPlot(Wo01_Ar5_Sig20_Pep90_Comparison, InputVarList, 2);
compareStrainStressProfilesPlot(Wo1_Ar1_Sig20_Pep90_Comparison,  InputVarList, 3);
compareStrainStressProfilesPlot(Wo1_Ar5_Sig20_Pep90_Comparison,  InputVarList, 4);
% Pepicelli, sigma = 20, strain = 0.95
figure
compareStrainStressProfilesPlot(Wo01_Ar1_Sig20_Pep95_Comparison, InputVarList, 1);
compareStrainStressProfilesPlot(Wo01_Ar5_Sig20_Pep95_Comparison, InputVarList, 2);
compareStrainStressProfilesPlot(Wo1_Ar1_Sig20_Pep95_Comparison,  InputVarList, 3);
compareStrainStressProfilesPlot(Wo1_Ar5_Sig20_Pep95_Comparison,  InputVarList, 4);
% Pepicelli, sigma = 60, strain = 0.8
figure
compareStrainStressProfilesPlot(Wo01_Ar1_Sig60_Pep80_Comparison, InputVarList, 1);
compareStrainStressProfilesPlot(Wo01_Ar5_Sig60_Pep80_Comparison, InputVarList, 2);
compareStrainStressProfilesPlot(Wo1_Ar1_Sig60_Pep80_Comparison,  InputVarList, 3);
compareStrainStressProfilesPlot(Wo1_Ar5_Sig60_Pep80_Comparison,  InputVarList, 4);
% Pepicelli, sigma = 60, strain = 0.9
figure
compareStrainStressProfilesPlot(Wo01_Ar1_Sig60_Pep90_Comparison, InputVarList, 1);
compareStrainStressProfilesPlot(Wo01_Ar5_Sig60_Pep90_Comparison, InputVarList, 2);
compareStrainStressProfilesPlot(Wo1_Ar1_Sig60_Pep90_Comparison,  InputVarList, 3);
compareStrainStressProfilesPlot(Wo1_Ar5_Sig60_Pep90_Comparison,  InputVarList, 4);
% Pepicelli, sigma = 60, strain = 0.95
figure
compareStrainStressProfilesPlot(Wo01_Ar1_Sig60_Pep95_Comparison, InputVarList, 1);
compareStrainStressProfilesPlot(Wo01_Ar5_Sig60_Pep95_Comparison, InputVarList, 2);
compareStrainStressProfilesPlot(Wo1_Ar1_Sig60_Pep95_Comparison,  InputVarList, 3);
compareStrainStressProfilesPlot(Wo1_Ar5_Sig60_Pep95_Comparison,  InputVarList, 4);


%% Fit circle for Capillary Pressure Microtensiometry

fraction = 0.75;

% Pepicelli, sigma = 20, strain = 0.8
Wo01_Ar1_Sig20_Pep80 = calculateStressCircle_Ref(Wo01_Ar1_Sig20_Pep80, VarList, fraction, 0);
Wo01_Ar5_Sig20_Pep80 = calculateStressCircle_Ref(Wo01_Ar5_Sig20_Pep80, VarList, fraction, 0);
Wo1_Ar1_Sig20_Pep80 = calculateStressCircle_Ref(Wo1_Ar1_Sig20_Pep80, VarList, fraction, 0);
Wo1_Ar5_Sig20_Pep80 = calculateStressCircle_Ref(Wo1_Ar5_Sig20_Pep80, VarList, fraction, 0);
Wo01_Ar1_Sig20_Pep80 = calculateStressCircle_Def(Wo01_Ar1_Sig20_Pep80, VarList, fraction, 0);
Wo01_Ar5_Sig20_Pep80 = calculateStressCircle_Def(Wo01_Ar5_Sig20_Pep80, VarList, fraction, 0);
Wo1_Ar1_Sig20_Pep80 = calculateStressCircle_Def(Wo1_Ar1_Sig20_Pep80, VarList, fraction, 0);
Wo1_Ar5_Sig20_Pep80 = calculateStressCircle_Def(Wo1_Ar5_Sig20_Pep80, VarList, fraction, 0);
Wo01_Ar1_Sig20_Pep80 = calculateModulusCircle(Wo01_Ar1_Sig20_Pep80, VarList);
Wo01_Ar5_Sig20_Pep80 = calculateModulusCircle(Wo01_Ar5_Sig20_Pep80, VarList);
Wo1_Ar1_Sig20_Pep80 = calculateModulusCircle(Wo1_Ar1_Sig20_Pep80, VarList);
Wo1_Ar5_Sig20_Pep80 = calculateModulusCircle(Wo1_Ar5_Sig20_Pep80, VarList);
% Pepicelli, sigma = 20, strain = 0.9
Wo01_Ar1_Sig20_Pep90 = calculateStressCircle_Ref(Wo01_Ar1_Sig20_Pep90, VarList, fraction, 0);
Wo01_Ar5_Sig20_Pep90 = calculateStressCircle_Ref(Wo01_Ar5_Sig20_Pep90, VarList, fraction, 0);
Wo1_Ar1_Sig20_Pep90 = calculateStressCircle_Ref(Wo1_Ar1_Sig20_Pep90, VarList, fraction, 0);
Wo1_Ar5_Sig20_Pep90 = calculateStressCircle_Ref(Wo1_Ar5_Sig20_Pep90, VarList, fraction, 0);
Wo01_Ar1_Sig20_Pep90 = calculateStressCircle_Def(Wo01_Ar1_Sig20_Pep90, VarList, fraction, 0);
Wo01_Ar5_Sig20_Pep90 = calculateStressCircle_Def(Wo01_Ar5_Sig20_Pep90, VarList, fraction, 0);
Wo1_Ar1_Sig20_Pep90 = calculateStressCircle_Def(Wo1_Ar1_Sig20_Pep90, VarList, fraction, 0);
Wo1_Ar5_Sig20_Pep90 = calculateStressCircle_Def(Wo1_Ar5_Sig20_Pep90, VarList, fraction, 0);
Wo01_Ar1_Sig20_Pep90 = calculateModulusCircle(Wo01_Ar1_Sig20_Pep90, VarList);
Wo01_Ar5_Sig20_Pep90 = calculateModulusCircle(Wo01_Ar5_Sig20_Pep90, VarList);
Wo1_Ar1_Sig20_Pep90 = calculateModulusCircle(Wo1_Ar1_Sig20_Pep90, VarList);
Wo1_Ar5_Sig20_Pep90 = calculateModulusCircle(Wo1_Ar5_Sig20_Pep90, VarList);
% Pepicelli, sigma = 20, strain = 0.95
Wo01_Ar1_Sig20_Pep95 = calculateStressCircle_Ref(Wo01_Ar1_Sig20_Pep95, VarList, fraction, 0);
Wo01_Ar5_Sig20_Pep95 = calculateStressCircle_Ref(Wo01_Ar5_Sig20_Pep95, VarList, fraction, 0);
Wo1_Ar1_Sig20_Pep95 = calculateStressCircle_Ref(Wo1_Ar1_Sig20_Pep95, VarList, fraction, 0);
Wo1_Ar5_Sig20_Pep95 = calculateStressCircle_Ref(Wo1_Ar5_Sig20_Pep95, VarList, fraction, 0);
Wo01_Ar1_Sig20_Pep95 = calculateStressCircle_Def(Wo01_Ar1_Sig20_Pep95, VarList, fraction, 0);
Wo01_Ar5_Sig20_Pep95 = calculateStressCircle_Def(Wo01_Ar5_Sig20_Pep95, VarList, fraction, 0);
Wo1_Ar1_Sig20_Pep95 = calculateStressCircle_Def(Wo1_Ar1_Sig20_Pep95, VarList, fraction, 0);
Wo1_Ar5_Sig20_Pep95 = calculateStressCircle_Def(Wo1_Ar5_Sig20_Pep95, VarList, fraction, 0);
Wo01_Ar1_Sig20_Pep95 = calculateModulusCircle(Wo01_Ar1_Sig20_Pep95, VarList);
Wo01_Ar5_Sig20_Pep95 = calculateModulusCircle(Wo01_Ar5_Sig20_Pep95, VarList);
Wo1_Ar1_Sig20_Pep95 = calculateModulusCircle(Wo1_Ar1_Sig20_Pep95, VarList);
Wo1_Ar5_Sig20_Pep95 = calculateModulusCircle(Wo1_Ar5_Sig20_Pep95, VarList);
% Pepicelli, sigma = 60, strain = 0.8
Wo01_Ar1_Sig60_Pep80 = calculateStressCircle_Ref(Wo01_Ar1_Sig60_Pep80, VarList, fraction, 0);
Wo01_Ar5_Sig60_Pep80 = calculateStressCircle_Ref(Wo01_Ar5_Sig60_Pep80, VarList, fraction, 0);
Wo1_Ar1_Sig60_Pep80 = calculateStressCircle_Ref(Wo1_Ar1_Sig60_Pep80, VarList, fraction, 0);
Wo1_Ar5_Sig60_Pep80 = calculateStressCircle_Ref(Wo1_Ar5_Sig60_Pep80, VarList, fraction, 0);
Wo01_Ar1_Sig60_Pep80 = calculateStressCircle_Def(Wo01_Ar1_Sig60_Pep80, VarList, fraction, 0);
Wo01_Ar5_Sig60_Pep80 = calculateStressCircle_Def(Wo01_Ar5_Sig60_Pep80, VarList, fraction, 0);
Wo1_Ar1_Sig60_Pep80 = calculateStressCircle_Def(Wo1_Ar1_Sig60_Pep80, VarList, fraction, 0);
Wo1_Ar5_Sig60_Pep80 = calculateStressCircle_Def(Wo1_Ar5_Sig60_Pep80, VarList, fraction, 0);
Wo01_Ar1_Sig60_Pep80 = calculateModulusCircle(Wo01_Ar1_Sig60_Pep80, VarList);
Wo01_Ar5_Sig60_Pep80 = calculateModulusCircle(Wo01_Ar5_Sig60_Pep80, VarList);
Wo1_Ar1_Sig60_Pep80 = calculateModulusCircle(Wo1_Ar1_Sig60_Pep80, VarList);
Wo1_Ar5_Sig60_Pep80 = calculateModulusCircle(Wo1_Ar5_Sig60_Pep80, VarList);
% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90 = calculateStressCircle_Ref(Wo01_Ar1_Sig60_Pep90, VarList, fraction, 0);
Wo01_Ar5_Sig60_Pep90 = calculateStressCircle_Ref(Wo01_Ar5_Sig60_Pep90, VarList, fraction, 0);
Wo1_Ar1_Sig60_Pep90 = calculateStressCircle_Ref(Wo1_Ar1_Sig60_Pep90, VarList, fraction, 0);
Wo1_Ar5_Sig60_Pep90 = calculateStressCircle_Ref(Wo1_Ar5_Sig60_Pep90, VarList, fraction, 0);
Wo01_Ar1_Sig60_Pep90 = calculateStressCircle_Def(Wo01_Ar1_Sig60_Pep90, VarList, fraction, 0);
Wo01_Ar5_Sig60_Pep90 = calculateStressCircle_Def(Wo01_Ar5_Sig60_Pep90, VarList, fraction, 0);
Wo1_Ar1_Sig60_Pep90 = calculateStressCircle_Def(Wo1_Ar1_Sig60_Pep90, VarList, fraction, 0);
Wo1_Ar5_Sig60_Pep90 = calculateStressCircle_Def(Wo1_Ar5_Sig60_Pep90, VarList, fraction, 0);
Wo01_Ar1_Sig60_Pep90 = calculateModulusCircle(Wo01_Ar1_Sig60_Pep90, VarList);
Wo01_Ar5_Sig60_Pep90 = calculateModulusCircle(Wo01_Ar5_Sig60_Pep90, VarList);
Wo1_Ar1_Sig60_Pep90 = calculateModulusCircle(Wo1_Ar1_Sig60_Pep90, VarList);
Wo1_Ar5_Sig60_Pep90 = calculateModulusCircle(Wo1_Ar5_Sig60_Pep90, VarList);
% Pepicelli, sigma = 60, strain = 0.95
Wo01_Ar1_Sig60_Pep95 = calculateStressCircle_Ref(Wo01_Ar1_Sig60_Pep95, VarList, fraction, 0);
Wo01_Ar5_Sig60_Pep95 = calculateStressCircle_Ref(Wo01_Ar5_Sig60_Pep95, VarList, fraction, 0);
Wo1_Ar1_Sig60_Pep95 = calculateStressCircle_Ref(Wo1_Ar1_Sig60_Pep95, VarList, fraction, 0);
Wo1_Ar5_Sig60_Pep95 = calculateStressCircle_Ref(Wo1_Ar5_Sig60_Pep95, VarList, fraction, 0);
Wo01_Ar1_Sig60_Pep95 = calculateStressCircle_Def(Wo01_Ar1_Sig60_Pep95, VarList, fraction, 0);
Wo01_Ar5_Sig60_Pep95 = calculateStressCircle_Def(Wo01_Ar5_Sig60_Pep95, VarList, fraction, 0);
Wo1_Ar1_Sig60_Pep95 = calculateStressCircle_Def(Wo1_Ar1_Sig60_Pep95, VarList, fraction, 0);
Wo1_Ar5_Sig60_Pep95 = calculateStressCircle_Def(Wo1_Ar5_Sig60_Pep95, VarList, fraction, 0);
Wo01_Ar1_Sig60_Pep95 = calculateModulusCircle(Wo01_Ar1_Sig60_Pep95, VarList);
Wo01_Ar5_Sig60_Pep95 = calculateModulusCircle(Wo01_Ar5_Sig60_Pep95, VarList);
Wo1_Ar1_Sig60_Pep95 = calculateModulusCircle(Wo1_Ar1_Sig60_Pep95, VarList);
Wo1_Ar5_Sig60_Pep95 = calculateModulusCircle(Wo1_Ar5_Sig60_Pep95, VarList);


%% Make tables with analyzed data

% Pepicelli, sigma = 20, strain = 0.8
makeTable(Wo01_Ar1_Sig20_Pep80, VarList);
makeTable(Wo01_Ar5_Sig20_Pep80, VarList);
makeTable(Wo1_Ar1_Sig20_Pep80, VarList);
makeTable(Wo1_Ar5_Sig20_Pep80, VarList);
% Pepicelli, sigma = 20, strain = 0.9
makeTable(Wo01_Ar1_Sig20_Pep90, VarList);
makeTable(Wo01_Ar5_Sig20_Pep90, VarList);
makeTable(Wo1_Ar1_Sig20_Pep90, VarList);
makeTable(Wo1_Ar5_Sig20_Pep90, VarList);
% Pepicelli, sigma = 20, strain = 0.95
makeTable(Wo01_Ar1_Sig20_Pep95, VarList);
makeTable(Wo01_Ar5_Sig20_Pep95, VarList);
makeTable(Wo1_Ar1_Sig20_Pep95, VarList);
makeTable(Wo1_Ar5_Sig20_Pep95, VarList);
% Pepicelli, sigma = 60, strain = 0.8
makeTable(Wo01_Ar1_Sig60_Pep80, VarList);
makeTable(Wo01_Ar5_Sig60_Pep80, VarList);
makeTable(Wo1_Ar1_Sig60_Pep80, VarList);
makeTable(Wo1_Ar5_Sig60_Pep80, VarList);
% Pepicelli, sigma = 60, strain = 0.9
makeTable(Wo01_Ar1_Sig60_Pep90, VarList);
makeTable(Wo01_Ar5_Sig60_Pep90, VarList);
makeTable(Wo1_Ar1_Sig60_Pep90, VarList);
makeTable(Wo1_Ar5_Sig60_Pep90, VarList);
% Pepicelli, sigma = 60, strain = 0.95
makeTable(Wo01_Ar1_Sig60_Pep95, VarList);
makeTable(Wo01_Ar5_Sig60_Pep95, VarList);
makeTable(Wo1_Ar1_Sig60_Pep95, VarList);
makeTable(Wo1_Ar5_Sig60_Pep95, VarList);


%% Plot apical strain parameters for inverse problem
x = [1 2 5 10 20 50 100]; y = x;
[X,Y] = meshgrid(x,y);

% Pepicelli, sigma = 20, strain = 0.8
Wo01_Ar1_Sig20_Pep80_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar1_Sig20_Pep80, VarList);
Wo01_Ar5_Sig20_Pep80_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar5_Sig20_Pep80, VarList);
Wo1_Ar1_Sig20_Pep80_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar1_Sig20_Pep80, VarList);
Wo1_Ar5_Sig20_Pep80_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar5_Sig20_Pep80, VarList);
% Pepicelli, sigma = 20, strain = 0.9
Wo01_Ar1_Sig20_Pep90_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar1_Sig20_Pep90, VarList);
Wo01_Ar5_Sig20_Pep90_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar5_Sig20_Pep90, VarList);
Wo1_Ar1_Sig20_Pep90_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar1_Sig20_Pep90, VarList);
Wo1_Ar5_Sig20_Pep90_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar5_Sig20_Pep90, VarList);
% Pepicelli, sigma = 20, strain = 0.95
Wo01_Ar1_Sig20_Pep95_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar1_Sig20_Pep95, VarList);
Wo01_Ar5_Sig20_Pep95_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar5_Sig20_Pep95, VarList);
Wo1_Ar1_Sig20_Pep95_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar1_Sig20_Pep95, VarList);
Wo1_Ar5_Sig20_Pep95_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar5_Sig20_Pep95, VarList);
% Pepicelli, sigma = 60, strain = 0.8
Wo01_Ar1_Sig60_Pep80_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar1_Sig60_Pep80, VarList);
Wo01_Ar5_Sig60_Pep80_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar5_Sig60_Pep80, VarList);
Wo1_Ar1_Sig60_Pep80_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar1_Sig60_Pep80, VarList);
Wo1_Ar5_Sig60_Pep80_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar5_Sig60_Pep80, VarList);
% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar1_Sig60_Pep90, VarList);
Wo01_Ar5_Sig60_Pep90_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar5_Sig60_Pep90, VarList);
Wo1_Ar1_Sig60_Pep90_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar1_Sig60_Pep90, VarList);
Wo1_Ar5_Sig60_Pep90_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar5_Sig60_Pep90, VarList);
% Pepicelli, sigma = 60, strain = 0.95
Wo01_Ar1_Sig60_Pep95_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar1_Sig60_Pep95, VarList);
Wo01_Ar5_Sig60_Pep95_ApexStrainTable = getStrainApexNonhomogen(Wo01_Ar5_Sig60_Pep95, VarList);
Wo1_Ar1_Sig60_Pep95_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar1_Sig60_Pep95, VarList);
Wo1_Ar5_Sig60_Pep95_ApexStrainTable = getStrainApexNonhomogen(Wo1_Ar5_Sig60_Pep95, VarList);

% Pepicelli, sigma = 20, strain = 0.8
figure
getStrainApexNonhomogenPlot_v2(Wo01_Ar1_Sig20_Pep80, VarList, Wo01_Ar1_Sig20_Pep80_ApexStrainTable, X, Y, 1);
getStrainApexNonhomogenPlot_v2(Wo01_Ar5_Sig20_Pep80, VarList, Wo01_Ar5_Sig20_Pep80_ApexStrainTable, X, Y, 2);
getStrainApexNonhomogenPlot_v2(Wo1_Ar1_Sig20_Pep80, VarList, Wo1_Ar1_Sig20_Pep80_ApexStrainTable, X, Y, 3);
getStrainApexNonhomogenPlot_v2(Wo1_Ar5_Sig20_Pep80, VarList, Wo1_Ar5_Sig20_Pep80_ApexStrainTable, X, Y, 4);
% Pepicelli, sigma = 20, strain = 0.9
figure
getStrainApexNonhomogenPlot_v2(Wo01_Ar1_Sig20_Pep90, VarList, Wo01_Ar1_Sig20_Pep90_ApexStrainTable, X, Y, 1);
getStrainApexNonhomogenPlot_v2(Wo01_Ar5_Sig20_Pep90, VarList, Wo01_Ar5_Sig20_Pep90_ApexStrainTable, X, Y, 2);
getStrainApexNonhomogenPlot_v2(Wo1_Ar1_Sig20_Pep90, VarList, Wo1_Ar1_Sig20_Pep90_ApexStrainTable, X, Y, 3);
getStrainApexNonhomogenPlot_v2(Wo1_Ar5_Sig20_Pep90, VarList, Wo1_Ar5_Sig20_Pep90_ApexStrainTable, X, Y, 4);
% Pepicelli, sigma = 20, strain = 0.95
figure
getStrainApexNonhomogenPlot_v2(Wo01_Ar1_Sig20_Pep95, VarList, Wo01_Ar1_Sig20_Pep95_ApexStrainTable, X, Y, 1);
getStrainApexNonhomogenPlot_v2(Wo01_Ar5_Sig20_Pep95, VarList, Wo01_Ar5_Sig20_Pep95_ApexStrainTable, X, Y, 2);
getStrainApexNonhomogenPlot_v2(Wo1_Ar1_Sig20_Pep95, VarList, Wo1_Ar1_Sig20_Pep95_ApexStrainTable, X, Y, 3);
getStrainApexNonhomogenPlot_v2(Wo1_Ar5_Sig20_Pep95, VarList, Wo1_Ar5_Sig20_Pep95_ApexStrainTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.8
figure
getStrainApexNonhomogenPlot_v2(Wo01_Ar1_Sig60_Pep80, VarList, Wo01_Ar1_Sig60_Pep80_ApexStrainTable, X, Y, 1);
getStrainApexNonhomogenPlot_v2(Wo01_Ar5_Sig60_Pep80, VarList, Wo01_Ar5_Sig60_Pep80_ApexStrainTable, X, Y, 2);
getStrainApexNonhomogenPlot_v2(Wo1_Ar1_Sig60_Pep80, VarList, Wo1_Ar1_Sig60_Pep80_ApexStrainTable, X, Y, 3);
getStrainApexNonhomogenPlot_v2(Wo1_Ar5_Sig60_Pep80, VarList, Wo1_Ar5_Sig60_Pep80_ApexStrainTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.9
figure
getStrainApexNonhomogenPlot_v2(Wo01_Ar1_Sig60_Pep90, VarList, Wo01_Ar1_Sig60_Pep90_ApexStrainTable, X, Y, 1);
getStrainApexNonhomogenPlot_v2(Wo01_Ar5_Sig60_Pep90, VarList, Wo01_Ar5_Sig60_Pep90_ApexStrainTable, X, Y, 2);
getStrainApexNonhomogenPlot_v2(Wo1_Ar1_Sig60_Pep90, VarList, Wo1_Ar1_Sig60_Pep90_ApexStrainTable, X, Y, 3);
getStrainApexNonhomogenPlot_v2(Wo1_Ar5_Sig60_Pep90, VarList, Wo1_Ar5_Sig60_Pep90_ApexStrainTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.95
figure
getStrainApexNonhomogenPlot_v2(Wo01_Ar1_Sig60_Pep95, VarList, Wo01_Ar1_Sig60_Pep95_ApexStrainTable, X, Y, 1);
getStrainApexNonhomogenPlot_v2(Wo01_Ar5_Sig60_Pep95, VarList, Wo01_Ar5_Sig60_Pep95_ApexStrainTable, X, Y, 2);
getStrainApexNonhomogenPlot_v2(Wo1_Ar1_Sig60_Pep95, VarList, Wo1_Ar1_Sig60_Pep95_ApexStrainTable, X, Y, 3);
getStrainApexNonhomogenPlot_v2(Wo1_Ar5_Sig60_Pep95, VarList, Wo1_Ar5_Sig60_Pep95_ApexStrainTable, X, Y, 4);


%% Plot apical stress parameters for inverse problem
x = [1 2 5 10 20 50 100]; y = x;
[X,Y] = meshgrid(x,y);

% Pepicelli, sigma = 20, strain = 0.8
Wo01_Ar1_Sig20_Pep80_ApexStressTable = getStressApexNonhomogen(Wo01_Ar1_Sig20_Pep80, VarList);
Wo01_Ar5_Sig20_Pep80_ApexStressTable = getStressApexNonhomogen(Wo01_Ar5_Sig20_Pep80, VarList);
Wo1_Ar1_Sig20_Pep80_ApexStressTable = getStressApexNonhomogen(Wo1_Ar1_Sig20_Pep80, VarList);
Wo1_Ar5_Sig20_Pep80_ApexStressTable = getStressApexNonhomogen(Wo1_Ar5_Sig20_Pep80, VarList);
% Pepicelli, sigma = 20, strain = 0.9
Wo01_Ar1_Sig20_Pep90_ApexStressTable = getStressApexNonhomogen(Wo01_Ar1_Sig20_Pep90, VarList);
Wo01_Ar5_Sig20_Pep90_ApexStressTable = getStressApexNonhomogen(Wo01_Ar5_Sig20_Pep90, VarList);
Wo1_Ar1_Sig20_Pep90_ApexStressTable = getStressApexNonhomogen(Wo1_Ar1_Sig20_Pep90, VarList);
Wo1_Ar5_Sig20_Pep90_ApexStressTable = getStressApexNonhomogen(Wo1_Ar5_Sig20_Pep90, VarList);
% Pepicelli, sigma = 20, strain = 0.95
Wo01_Ar1_Sig20_Pep95_ApexStressTable = getStressApexNonhomogen(Wo01_Ar1_Sig20_Pep95, VarList);
Wo01_Ar5_Sig20_Pep95_ApexStressTable = getStressApexNonhomogen(Wo01_Ar5_Sig20_Pep95, VarList);
Wo1_Ar1_Sig20_Pep95_ApexStressTable = getStressApexNonhomogen(Wo1_Ar1_Sig20_Pep95, VarList);
Wo1_Ar5_Sig20_Pep95_ApexStressTable = getStressApexNonhomogen(Wo1_Ar5_Sig20_Pep95, VarList);
% Pepicelli, sigma = 60, strain = 0.8
Wo01_Ar1_Sig60_Pep80_ApexStressTable = getStressApexNonhomogen(Wo01_Ar1_Sig60_Pep80, VarList);
Wo01_Ar5_Sig60_Pep80_ApexStressTable = getStressApexNonhomogen(Wo01_Ar5_Sig60_Pep80, VarList);
Wo1_Ar1_Sig60_Pep80_ApexStressTable = getStressApexNonhomogen(Wo1_Ar1_Sig60_Pep80, VarList);
Wo1_Ar5_Sig60_Pep80_ApexStressTable = getStressApexNonhomogen(Wo1_Ar5_Sig60_Pep80, VarList);
% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90_ApexStressTable = getStressApexNonhomogen(Wo01_Ar1_Sig60_Pep90, VarList);
Wo01_Ar5_Sig60_Pep90_ApexStressTable = getStressApexNonhomogen(Wo01_Ar5_Sig60_Pep90, VarList);
Wo1_Ar1_Sig60_Pep90_ApexStressTable = getStressApexNonhomogen(Wo1_Ar1_Sig60_Pep90, VarList);
Wo1_Ar5_Sig60_Pep90_ApexStressTable = getStressApexNonhomogen(Wo1_Ar5_Sig60_Pep90, VarList);
% Pepicelli, sigma = 60, strain = 0.95
Wo01_Ar1_Sig60_Pep95_ApexStressTable = getStressApexNonhomogen(Wo01_Ar1_Sig60_Pep95, VarList);
Wo01_Ar5_Sig60_Pep95_ApexStressTable = getStressApexNonhomogen(Wo01_Ar5_Sig60_Pep95, VarList);
Wo1_Ar1_Sig60_Pep95_ApexStressTable = getStressApexNonhomogen(Wo1_Ar1_Sig60_Pep95, VarList);
Wo1_Ar5_Sig60_Pep95_ApexStressTable = getStressApexNonhomogen(Wo1_Ar5_Sig60_Pep95, VarList);

% Pepicelli, sigma = 20, strain = 0.8
figure
getStressApexNonhomogenPlot_v2(Wo01_Ar1_Sig20_Pep80, VarList, Wo01_Ar1_Sig20_Pep80_ApexStressTable, X, Y, 1);
getStressApexNonhomogenPlot_v2(Wo01_Ar5_Sig20_Pep80, VarList, Wo01_Ar5_Sig20_Pep80_ApexStressTable, X, Y, 2);
getStressApexNonhomogenPlot_v2(Wo1_Ar1_Sig20_Pep80, VarList, Wo1_Ar1_Sig20_Pep80_ApexStressTable, X, Y, 3);
getStressApexNonhomogenPlot_v2(Wo1_Ar5_Sig20_Pep80, VarList, Wo1_Ar5_Sig20_Pep80_ApexStressTable, X, Y, 4);
% Pepicelli, sigma = 20, strain = 0.9
figure
getStressApexNonhomogenPlot_v2(Wo01_Ar1_Sig20_Pep90, VarList, Wo01_Ar1_Sig20_Pep90_ApexStressTable, X, Y, 1);
getStressApexNonhomogenPlot_v2(Wo01_Ar5_Sig20_Pep90, VarList, Wo01_Ar5_Sig20_Pep90_ApexStressTable, X, Y, 2);
getStressApexNonhomogenPlot_v2(Wo1_Ar1_Sig20_Pep90, VarList, Wo1_Ar1_Sig20_Pep90_ApexStressTable, X, Y, 3);
getStressApexNonhomogenPlot_v2(Wo1_Ar5_Sig20_Pep90, VarList, Wo1_Ar5_Sig20_Pep90_ApexStressTable, X, Y, 4);
% Pepicelli, sigma = 20, strain = 0.95
figure
getStressApexNonhomogenPlot_v2(Wo01_Ar1_Sig20_Pep95, VarList, Wo01_Ar1_Sig20_Pep95_ApexStressTable, X, Y, 1);
getStressApexNonhomogenPlot_v2(Wo01_Ar5_Sig20_Pep95, VarList, Wo01_Ar5_Sig20_Pep95_ApexStressTable, X, Y, 2);
getStressApexNonhomogenPlot_v2(Wo1_Ar1_Sig20_Pep95, VarList, Wo1_Ar1_Sig20_Pep95_ApexStressTable, X, Y, 3);
getStressApexNonhomogenPlot_v2(Wo1_Ar5_Sig20_Pep95, VarList, Wo1_Ar5_Sig20_Pep95_ApexStressTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.8
figure
getStressApexNonhomogenPlot_v2(Wo01_Ar1_Sig60_Pep80, VarList, Wo01_Ar1_Sig60_Pep80_ApexStressTable, X, Y, 1);
getStressApexNonhomogenPlot_v2(Wo01_Ar5_Sig60_Pep80, VarList, Wo01_Ar5_Sig60_Pep80_ApexStressTable, X, Y, 2);
getStressApexNonhomogenPlot_v2(Wo1_Ar1_Sig60_Pep80, VarList, Wo1_Ar1_Sig60_Pep80_ApexStressTable, X, Y, 3);
getStressApexNonhomogenPlot_v2(Wo1_Ar5_Sig60_Pep80, VarList, Wo1_Ar5_Sig60_Pep80_ApexStressTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.9
figure
getStressApexNonhomogenPlot_v2(Wo01_Ar1_Sig60_Pep90, VarList, Wo01_Ar1_Sig60_Pep90_ApexStressTable, X, Y, 1);
getStressApexNonhomogenPlot_v2(Wo01_Ar5_Sig60_Pep90, VarList, Wo01_Ar5_Sig60_Pep90_ApexStressTable, X, Y, 2);
getStressApexNonhomogenPlot_v2(Wo1_Ar1_Sig60_Pep90, VarList, Wo1_Ar1_Sig60_Pep90_ApexStressTable, X, Y, 3);
getStressApexNonhomogenPlot_v2(Wo1_Ar5_Sig60_Pep90, VarList, Wo1_Ar5_Sig60_Pep90_ApexStressTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.95
figure
getStressApexNonhomogenPlot_v2(Wo01_Ar1_Sig60_Pep95, VarList, Wo01_Ar1_Sig60_Pep95_ApexStressTable, X, Y, 1);
getStressApexNonhomogenPlot_v2(Wo01_Ar5_Sig60_Pep95, VarList, Wo01_Ar5_Sig60_Pep95_ApexStressTable, X, Y, 2);
getStressApexNonhomogenPlot_v2(Wo1_Ar1_Sig60_Pep95, VarList, Wo1_Ar1_Sig60_Pep95_ApexStressTable, X, Y, 3);
getStressApexNonhomogenPlot_v2(Wo1_Ar5_Sig60_Pep95, VarList, Wo1_Ar5_Sig60_Pep95_ApexStressTable, X, Y, 4);


%% Plot modulus parameters
x = [1 2 5 10 20 50 100]; y = x;
[X,Y] = meshgrid(x,y);

% Pepicelli, sigma = 20, strain = 0.8
Wo01_Ar1_Sig20_Pep80_ModulusTable = getModulusNonhomogen(Wo01_Ar1_Sig20_Pep80, VarList);
Wo01_Ar5_Sig20_Pep80_ModulusTable = getModulusNonhomogen(Wo01_Ar5_Sig20_Pep80, VarList);
Wo1_Ar1_Sig20_Pep80_ModulusTable = getModulusNonhomogen(Wo1_Ar1_Sig20_Pep80, VarList);
Wo1_Ar5_Sig20_Pep80_ModulusTable = getModulusNonhomogen(Wo1_Ar5_Sig20_Pep80, VarList);
% Pepicelli, sigma = 20, strain = 0.9
Wo01_Ar1_Sig20_Pep90_ModulusTable = getModulusNonhomogen(Wo01_Ar1_Sig20_Pep90, VarList);
Wo01_Ar5_Sig20_Pep90_ModulusTable = getModulusNonhomogen(Wo01_Ar5_Sig20_Pep90, VarList);
Wo1_Ar1_Sig20_Pep90_ModulusTable = getModulusNonhomogen(Wo1_Ar1_Sig20_Pep90, VarList);
Wo1_Ar5_Sig20_Pep90_ModulusTable = getModulusNonhomogen(Wo1_Ar5_Sig20_Pep90, VarList);
% Pepicelli, sigma = 20, strain = 0.95
Wo01_Ar1_Sig20_Pep95_ModulusTable = getModulusNonhomogen(Wo01_Ar1_Sig20_Pep95, VarList);
Wo01_Ar5_Sig20_Pep95_ModulusTable = getModulusNonhomogen(Wo01_Ar5_Sig20_Pep95, VarList);
Wo1_Ar1_Sig20_Pep95_ModulusTable = getModulusNonhomogen(Wo1_Ar1_Sig20_Pep95, VarList);
Wo1_Ar5_Sig20_Pep95_ModulusTable = getModulusNonhomogen(Wo1_Ar5_Sig20_Pep95, VarList);
% Pepicelli, sigma = 60, strain = 0.8
Wo01_Ar1_Sig60_Pep80_ModulusTable = getModulusNonhomogen(Wo01_Ar1_Sig60_Pep80, VarList);
Wo01_Ar5_Sig60_Pep80_ModulusTable = getModulusNonhomogen(Wo01_Ar5_Sig60_Pep80, VarList);
Wo1_Ar1_Sig60_Pep80_ModulusTable = getModulusNonhomogen(Wo1_Ar1_Sig60_Pep80, VarList);
Wo1_Ar5_Sig60_Pep80_ModulusTable = getModulusNonhomogen(Wo1_Ar5_Sig60_Pep80, VarList);
% Pepicelli, sigma = 60, strain = 0.9
Wo01_Ar1_Sig60_Pep90_ModulusTable = getModulusNonhomogen(Wo01_Ar1_Sig60_Pep90, VarList);
Wo01_Ar5_Sig60_Pep90_ModulusTable = getModulusNonhomogen(Wo01_Ar5_Sig60_Pep90, VarList);
Wo1_Ar1_Sig60_Pep90_ModulusTable = getModulusNonhomogen(Wo1_Ar1_Sig60_Pep90, VarList);
Wo1_Ar5_Sig60_Pep90_ModulusTable = getModulusNonhomogen(Wo1_Ar5_Sig60_Pep90, VarList);
% Pepicelli, sigma = 60, strain = 0.95
Wo01_Ar1_Sig60_Pep95_ModulusTable = getModulusNonhomogen(Wo01_Ar1_Sig60_Pep95, VarList);
Wo01_Ar5_Sig60_Pep95_ModulusTable = getModulusNonhomogen(Wo01_Ar5_Sig60_Pep95, VarList);
Wo1_Ar1_Sig60_Pep95_ModulusTable = getModulusNonhomogen(Wo1_Ar1_Sig60_Pep95, VarList);
Wo1_Ar5_Sig60_Pep95_ModulusTable = getModulusNonhomogen(Wo1_Ar5_Sig60_Pep95, VarList);

% Pepicelli, sigma = 20, strain = 0.8
figure
getModulusNonhomogenPlot_v2(Wo01_Ar1_Sig20_Pep80, VarList, Wo01_Ar1_Sig20_Pep80_ModulusTable, X, Y, 1);
getModulusNonhomogenPlot_v2(Wo01_Ar5_Sig20_Pep80, VarList, Wo01_Ar5_Sig20_Pep80_ModulusTable, X, Y, 2);
getModulusNonhomogenPlot_v2(Wo1_Ar1_Sig20_Pep80, VarList, Wo1_Ar1_Sig20_Pep80_ModulusTable, X, Y, 3);
getModulusNonhomogenPlot_v2(Wo1_Ar5_Sig20_Pep80, VarList, Wo1_Ar5_Sig20_Pep80_ModulusTable, X, Y, 4);
% Pepicelli, sigma = 20, strain = 0.9
figure
getModulusNonhomogenPlot_v2(Wo01_Ar1_Sig20_Pep90, VarList, Wo01_Ar1_Sig20_Pep90_ModulusTable, X, Y, 1);
getModulusNonhomogenPlot_v2(Wo01_Ar5_Sig20_Pep90, VarList, Wo01_Ar5_Sig20_Pep90_ModulusTable, X, Y, 2);
getModulusNonhomogenPlot_v2(Wo1_Ar1_Sig20_Pep90, VarList, Wo1_Ar1_Sig20_Pep90_ModulusTable, X, Y, 3);
getModulusNonhomogenPlot_v2(Wo1_Ar5_Sig20_Pep90, VarList, Wo1_Ar5_Sig20_Pep90_ModulusTable, X, Y, 4);
% Pepicelli, sigma = 20, strain = 0.95
figure
getModulusNonhomogenPlot_v2(Wo01_Ar1_Sig20_Pep95, VarList, Wo01_Ar1_Sig20_Pep95_ModulusTable, X, Y, 1);
getModulusNonhomogenPlot_v2(Wo01_Ar5_Sig20_Pep95, VarList, Wo01_Ar5_Sig20_Pep95_ModulusTable, X, Y, 2);
getModulusNonhomogenPlot_v2(Wo1_Ar1_Sig20_Pep95, VarList, Wo1_Ar1_Sig20_Pep95_ModulusTable, X, Y, 3);
getModulusNonhomogenPlot_v2(Wo1_Ar5_Sig20_Pep95, VarList, Wo1_Ar5_Sig20_Pep95_ModulusTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.8
figure
getModulusNonhomogenPlot_v2(Wo01_Ar1_Sig60_Pep80, VarList, Wo01_Ar1_Sig60_Pep80_ModulusTable, X, Y, 1);
getModulusNonhomogenPlot_v2(Wo01_Ar5_Sig60_Pep80, VarList, Wo01_Ar5_Sig60_Pep80_ModulusTable, X, Y, 2);
getModulusNonhomogenPlot_v2(Wo1_Ar1_Sig60_Pep80, VarList, Wo1_Ar1_Sig60_Pep80_ModulusTable, X, Y, 3);
getModulusNonhomogenPlot_v2(Wo1_Ar5_Sig60_Pep80, VarList, Wo1_Ar5_Sig60_Pep80_ModulusTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.9
figure
getModulusNonhomogenPlot_v2(Wo01_Ar1_Sig60_Pep90, VarList, Wo01_Ar1_Sig60_Pep90_ModulusTable, X, Y, 1);
getModulusNonhomogenPlot_v2(Wo01_Ar5_Sig60_Pep90, VarList, Wo01_Ar5_Sig60_Pep90_ModulusTable, X, Y, 2);
getModulusNonhomogenPlot_v2(Wo1_Ar1_Sig60_Pep90, VarList, Wo1_Ar1_Sig60_Pep90_ModulusTable, X, Y, 3);
getModulusNonhomogenPlot_v2(Wo1_Ar5_Sig60_Pep90, VarList, Wo1_Ar5_Sig60_Pep90_ModulusTable, X, Y, 4);
% Pepicelli, sigma = 60, strain = 0.95
figure
getModulusNonhomogenPlot_v2(Wo01_Ar1_Sig60_Pep95, VarList, Wo01_Ar1_Sig60_Pep95_ModulusTable, X, Y, 1);
getModulusNonhomogenPlot_v2(Wo01_Ar5_Sig60_Pep95, VarList, Wo01_Ar5_Sig60_Pep95_ModulusTable, X, Y, 2);
getModulusNonhomogenPlot_v2(Wo1_Ar1_Sig60_Pep95, VarList, Wo1_Ar1_Sig60_Pep95_ModulusTable, X, Y, 3);
getModulusNonhomogenPlot_v2(Wo1_Ar5_Sig60_Pep95, VarList, Wo1_Ar5_Sig60_Pep95_ModulusTable, X, Y, 4);

