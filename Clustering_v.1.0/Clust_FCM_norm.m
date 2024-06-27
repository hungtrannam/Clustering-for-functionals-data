clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading

addpath Data Figure_Output Package/Clust Package/Vis/ Package/Val/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting

% Data
% abnormal_params = {
%     {[0.15, 2.5], sqrt([.4, .4])} 
%     {[0.3, 2.2], sqrt([.4, .4])}
%     };
% 
% [Data, param.x, param.truelabels] = SimPDFAbnormal( ...
%     { ...
%     linspace(-1, 1, 10*90), ...
%     linspace(3, 4, 10)}, ...
%     sqrt([.5, .5]));

[Data, param.x, param.truelabels] = SimPDFAbnormal( ...
    { ...
    linspace(0, 3, 5*15), ...
    linspace(5, 6, 5)}, ...
    sqrt([.8, .5]));

% FCM
param.maxIter   = 200;                            % Maximum number of iterations
param.mFuzzy    = 2;                              % Fuzziness parameter
param.epsilon   = 1e-10;                          % Convergence criterion
param.kClust    = 2;                              % Number of clusters
param.val       = 2;                              % Validation parameter
param.FvIni     = 2;                              % Initial value for fuzziness parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering via alg

results = IFCM_(Data, param, 'Visualize', 'True'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
% figure;
% heatmap(results.Cluster.U);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation
results = validityClustering(results, param);
















% options = fcmOptions(...
% NumClusters=2,...
% Exponent=2.0,...
% Verbose=true,
% DistanceMetric = 'fmle');
% [centers,U,objFcn,info] = fcm(Data',options);