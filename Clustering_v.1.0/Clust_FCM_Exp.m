clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading

addpath Data Figure_Output Package/Clust Package/Vis/ Package/Val/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting

% Data
[Data, param.x, param.truelabels] = SimPDFAbnormal_Expo( ...
    { ...
    linspace(5, 6, 10*8), ...
    linspace(7, 7.5, 3), ...
    });
% FCM
param.maxIter = 1000;
param.mFuzzy = 2;
param.epsilon = 1e-10;
param.kClust = 2;
param.val = 2;
param.FvIni = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering via alg

results = IFCM_(Data, param); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
% figure;
% heatmap(results.Cluster.U);


h = PlotPDFeachIteration(Data, results.Cluster.IDX, param.x); hold on;
plot(param.x,results.Data.fv, "LineWidth", 3, "DisplayName", "Representative PDF"); hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation
results = validityClustering(results, param);
















% options = fcmOptions(...
% NumClusters=2,...
% Exponent=2.0,...
% Verbose=true,
% DistanceMetric = 'fmle');
% [centers,U,objFcn,info] = fcm(Data',options);