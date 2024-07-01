clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading

addpath Data Figure_Output Package/Clust Package/Vis/ Package/Val/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting

% Data
param.h         = .03;
param.x         = -5: param.h : 10;
param.val       = 2;

abnormal_params = {
    {[0.15, 2.5], sqrt([.4, .4])}
    {[0.3, 2.2], sqrt([.4, .4])}
    };

% Simulate data and true labels using the SimPDFAbnormal function
[Data, param.truelabels] = SimPDFAbnormal( ...
    { ...
    linspace(0, 1, 10*50), ...
    linspace(4, 5, 10), ...
    linspace(6, 6.5, 10*10)}, ...
    sqrt([.5, .5 .5]), ...
    param.x);


% FCM
param.L       = 0.2;
param.epsilon = 1e-10;
param.maxIter = 500;
param.val     = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering via alg

results       = ERACF_(Data, param); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

figure;
PlotPDFeachIteration(Data, results.Cluster.IDX, param.x); legend("Location","eastoutside"); hold on;
plot(param.x,results.Data.fv, "LineWidth", 3, "DisplayName", "Representative PDF"); hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation

results       = validityClustering(results, param);

