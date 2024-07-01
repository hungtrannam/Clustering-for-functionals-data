%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML110
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%
clc;
clear;
close all;

%% Load Data

addpath Data Figure_Output Package/Clust Package/Vis/ Package/Val

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting

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
    linspace(0, 1, 10*100), ...
    linspace(4, 5, 10)}, ...
    sqrt([.5, .5]), ...
    param.x);

%% Run DBSCAN Clustering Algorithm
param.epsilon = 1;
param.MinPts  = 10;
param.val = 2;

results = DBSCAN_(Data', param);

%% Plot Results
PlotPDFeachIteration(Data, results.Cluster.IDX+1, param.x);
title(['DBSCAN Clustering (\epsilon = ' num2str(param.epsilon) ', MinPts = ' num2str(param.MinPts) ')']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation
results = validityClustering(results, param);
