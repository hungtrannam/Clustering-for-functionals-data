clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading

addpath Data Figure_Output Package/Clust Package/Vis/ Package/Val/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting

% Data

abnormal_params = {
    {[0.15, 2.7], sqrt([.2, .2])} 
    {[0.3, 2.9], sqrt([.2, .2])}
    };

[Data, param.x, param.truelabels] = SimPDFAbnormal( ...
    { ...
    linspace(0, 0.5, 50), ...
    linspace(2.5, 3, 50),
    }, ...
    sqrt([.2, .2]), abnormal_params);
    % abnormal_params ...
 

% PCM
param.maxIter    = 1000;
param.mFuzzy     = 2;
param.epsilon    = 1e-10;
param.kClust     = 2;
param.val        = 3;
param.FvIni      = 3;
param.K          = 1;
param.alphaCut   = .2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering via alg

results          = PCM_(Data, param); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
figure;
heatmap(results.Cluster.U);

h = PlotPDFeachIteration(Data, results.Cluster.IDX+1, param.x); hold on;
plot(param.x,results.Data.fv, "LineWidth", 3, "DisplayName", "Representative PDF"); hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation
results = validityClustering(results, param);



% options = fcmOptions(...
% NumClusters=2,...
% Exponent=2.0,...
% Verbose=true);
% [centers,U,objFcn,info] = fcm(Data',options);