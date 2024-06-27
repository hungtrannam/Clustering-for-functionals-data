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

abnormal_params = {
    {[0.15, 2.7], sqrt([.4, .4])} 
    {[0.3, 2.9], sqrt([.4, .4])}
    };

[Data, param.x, param.truelabels] = SimPDFAbnormal( ...
    { ...
    linspace(0, 0.5, 5*25), ...
    linspace(2.5, 3, 5),
    }, ...
    sqrt([.4, .4]), abnormal_params);
    % abnormal_params ...


%% Run DBSCAN Clustering Algorithm
param.epsilon = 1;
param.MinPts  = 5;
param.val = 2;

results = DBSCAN_(Data', param);

%% Plot Results
PlotPDFeachIteration(Data, results.Cluster.IDX+1, param.x);
title(['DBSCAN Clustering (\epsilon = ' num2str(param.epsilon) ', MinPts = ' num2str(param.MinPts) ')']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation
results = validityClustering(results, param);
