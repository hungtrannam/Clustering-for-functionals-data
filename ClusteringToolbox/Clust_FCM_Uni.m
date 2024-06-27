clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading

addpath Data Figure_Output Package/Clust Package/Vis/ Package/Val/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting

% Data
[Data, param.x, param.truelabels] = SimPDFAbnormal_Uniform( ...
    { ...
    linspace(2, 6, 10*20), ...
    linspace(7, 7.5, 10), ...
    }, ...
    sqrt([2, .2]));
% FCM
param.maxIter   = 200;                            % Maximum number of iterations
param.mFuzzy    = 2;                              % Fuzziness parameter
param.epsilon   = 1e-10;                          % Convergence criterion
param.kClust    = 2;                              % Number of clusters
param.val       = 2;                              % Validation parameter
param.FvIni     = 3;                              % Initial value for fuzziness parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering via alg

results = IFCM_(Data, param, 'Visualize', 'True'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
% figure;
% heatmap(results.Cluster.U);


% h = figure;
% subplot(1,2,1);
% PlotPDFeachIteration(Data, results.Cluster.IDX, param.x); hold on;
% plot(param.x,results.Data.fv, "LineWidth", 3, "DisplayName", "Representative PDF"); hold off;
% subplot(1,2,2);
% PlotCDEsContour(Data, results.Data.fv, results.Cluster.IDX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation
results = validityClustering(results, param);
















% options = fcmOptions(...
% NumClusters=2,...
% Exponent=2.0,...
% Verbose=true,
% DistanceMetric = 'fmle');
% [centers,U,objFcn,info] = fcm(Data',options);