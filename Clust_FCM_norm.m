%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hung Tran-Nam
% Date: Jun 23, 2024
% Description: This script performs data clustering using the IFCM_ algorithm,
%              plots the results, and validates the clustering outcome.
%
% Copyright (c) 2024, Van Lang University
% All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Clear the workspace, close all figures, and clear command window

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the statistics package
% pkg load statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading

% Add specified directories to the MATLAB path for access to custom functions and scripts
addpath Data Figure_Output Package/Clust Package/Vis/ Package/Val/ Package/ToolBox/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting

param.h         = .03;
param.x         = -5: param.h : 10;

abnormal_params = {
    {[0.15, 2.5], sqrt([.4, .4])}
    {[0.3, 2.2], sqrt([.4, .4])}
    };

% Simulate data and true labels using the SimPDFAbnormal function
[Data, param.truelabels] = SimPDFAbnormal( ...
    { ...
    linspace(0, 1, 10*100), ...
    linspace(4, 5, 10)}, ...
    sqrt([.5, .5]), ...
    param.x);




% FCM (Fuzzy C-Means) parameters
param.maxIter   = 200;                            % Maximum number of iterations
param.mFuzzy    = 2;                              % Fuzziness parameter
param.epsilon   = 1e-10;                          % Convergence criterion
param.kClust    = 2;                              % Number of clusters
param.val       = 2;                              % Validation parameter
param.FvIni     = 3;                              % Initial value for fuzziness parameter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering via algorithm

% Perform clustering using the IFCM_ function
results         = IFCM_(Data, param);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

% figure;
% heatmap(results.Cluster.U);

% Plot the PDF for each iteration of the clustering process
figure(1);
PlotPDFeachIteration(Data, results.Cluster.IDX, param.x); hold on;
plot(param.x,results.Data.fv, "LineWidth", 3, "DisplayName", "Representative PDF");
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation

% Validate the clustering results using the validityClustering function
results         = validityClustering(results, param);

% Determine the best number of clusters
% param.Klist     = 2:5;
% bestK           = evalclusters_(Data,param, param.Klist);
% param.kClust    = bestK;
% 



% Example FCM options (uncomment and modify if needed)
% options = fcmOptions(...
% NumClusters=2,...
% Exponent=2.0,...
% Verbose=true,
% DistanceMetric = 'fmle');
% [centers,U,objFcn,info] = fcm(Data',options);
