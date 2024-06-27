clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading

addpath Data Figure_Output Package/Clust Package/Vis/ Package/Val/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting

% Data
% [Data, param.x, param.truelabels] = SimPDFAnormal( ...
%     { ...
%     0.35:0.003:0.45, ...
%     0.75:0.003:0.8, ...
%     0.10, 0.55, 0.56, 0.57, 0.90 ...
%     }, ...
%     sqrt([0.008, 0.008, 0.008,0.008,0.008,0.008,0.008,0.008]));

abnormal_params = {
    {[0.15, 2.5], sqrt([.1, .2])} 
    {[0.3, 2.2], sqrt([.1, .2])}
    };

[Data, param.x, param.truelabels] = SimPDFAbnormal( ...
    { ...
    0.15:.005:0.45, ...
    [2, 2.1, 2.2, 2.3, 2.5],
    }, ...
    sqrt([.1, .2]), ...
    abnormal_params ...
    );


% FCM
param.maxIter = 500;
param.mFuzzy = 2;
param.epsilon = 1e-10;
param.kClust = 2;
param.val = 2;
param.K = 1;
param.cF = 1;
param.cP = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering via alg

results = PFCM_(Data, param); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
figure;
heatmap(results.Cluster.U);


h = PlotPDFeachIteration(Data, results.Cluster.IDX, param.x); hold on;
plot(param.x,results.Data.fv, "LineWidth", 3, "DisplayName", "Representative PDF"); hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validation
results = validityClustering(results, param);
















% options = fcmOptions(...
% NumClusters=2,...
% Exponent=2.0,...
% Verbose=true);
% [centers,U,objFcn,info] = fcm(Data',options);