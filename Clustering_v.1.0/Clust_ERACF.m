clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading

addpath Data Figure_Output Package/Clust Package/Vis/ Package/Val/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting

% Data
% [Data, param.x, param.truelabels] = SimPDFAbnormal_Uniform( ...
%     { ...
%     linspace(2, 6, 10*9), ...
%     linspace(7, 7.5, 10), ...
%     }, ...
%     sqrt([2, .2]));


[Data, param.x, param.truelabels] = SimPDFAbnormal( ...
    { ...
    linspace(0, 3, 5*10), ...
    linspace(5, 6, 5)}, ...
    sqrt([.8, .5]));
% FCM
param.L       = 0.1;
param.epsilon = 1e-10;
param.maxIter = 500;
param.val     = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering via alg

results       = ERACF_(Data, param, 'Visualize', 'True'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

figure;
PlotPDFeachIteration(Data, results.Cluster.IDX, param.x); legend("Location","eastoutside"); hold on;
plot(param.x,results.Data.fv, "LineWidth", 3, "DisplayName", "Representative PDF"); hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation

results       = validityClustering(results, param);
















% options = fcmOptions(...
% NumClusters=2,...
% Exponent=2.0,...
% Verbose=true,
% DistanceMetric = 'fmle');
% [centers,U,objFcn,info] = fcm(Data',options);
