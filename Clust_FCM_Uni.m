clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading

addpath Data Figure_Output Package/Clust Package/Vis/ Package/Val/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting

param.h         = 0.001;
param.x         = 0: param.h : 15;


% Data
[Data, param.truelabels] = SimPDFAbnormal_Uniform( ...
    { ...
    linspace(2, 5, 5*50), ...
    linspace(8, 8.5, 5), ...
    }, ...
    sqrt([1, .8]), ...
    param.x);


% FCM
param.maxIter   = 200;                            % Maximum number of iterations
param.mFuzzy    = 2;                              % Fuzziness parameter
param.epsilon   = 1e-10;                          % Convergence criterion
param.kClust    = 2;                              % Number of clusters
param.val       = 2;                              % Validation parameter
param.FvIni     = 3;                              % Initial value for fuzziness parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering via alg

results = IFCM_(Data, param); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

figure(1);
PlotPDFeachIteration(Data, results.Cluster.IDX, param.x); hold on;
plot(param.x,results.Data.fv, "LineWidth", 3, "DisplayName", "Representative PDF");
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation
results = validityClustering(results, param);
















% options = fcmOptions(...
% NumClusters=2,...
% Exponent=2.0,...
% Verbose=true,
% DistanceMetric = 'fmle');
% [centers,U,objFcn,info] = fcm(Data',options);