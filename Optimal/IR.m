%calling function to search the optimal number of clusters
close all
clear all


addpath /home/hung-tran-nam/Hung_matlab/Clustering_v.1.0/Data
addpath /home/hung-tran-nam/Hung_matlab/Clustering_v.1.0/Package/Clust
addpath /home/hung-tran-nam/Hung_matlab/Clustering_v.1.0/Package/Val

%%%%%%%%%%%%%%%%%%%%%%%
% Setting

% FCM
param.maxIter    = 1000;
param.mFuzzy     = 2;
param.epsilon    = 1e-10;
param.kClust     = 2;
param.val        = 3;
param.FvIni      = 3;


ment             = [];
nkmin            = 10;
nkmax            = 100;
index            = 1;

L                = nkmin:10:nkmax;

%%%%%%%%%%%%%%%%%%%%%%%
% Loot

for cln = L
    param.IR    = cln;
    [Data, param.x, param.truelabels] = SimPDFAbnormal( ...
    { ...
    linspace(6, 7, 5*param.IR), ...
    linspace(9, 10, 5)}, ...
    sqrt([0.25, 0.25]));
    % abnormal_params ...
  
    results     = IFCM_(Data, param);
    results     = validityClustering(results,param);

    ment{index} = results.validity;
    index       = index + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%
% Plotting


if param.val == 1 || param.val == 3
    PC=[]; CE=[]; SC=[]; XB=[];
    for i = 1:index-1
        PC = [PC ment{i}.PC];
        CE = [CE ment{i}.CE];
        SC = [SC ment{i}.SC];
        XB = [XB ment{i}.XB];
    end
    figure(1)
    clf
    subplot(4,1,1), plot(L,PC, 'k-o')
    title('Partition Coefficient (PC)')
    subplot(4,1,2), plot(L,CE,'r-o')
    title('Classification Entropy (CE)')
    subplot(4,1,3), plot(L,SC,'g-o')
    title('Separation Index (S)')
    subplot(4,1,4), plot(L,XB, 'b-o')
    title('Xie and Beni Index (XB)')

end

if param.val == 2 || param.val == 3
    Gmean=[]; ARI=[]; RI=[];

    for i = 1:index-1
        Gmean = [Gmean ment{i}.Gmean];
        ARI   = [ARI ment{i}.ARI];
        RI    = [RI ment{i}.RI];

    end
    figure(2)
    clf
    subplot(3,1,1), plot(L, Gmean, 'r-o')
    title('Gmean')
    subplot(3,1,2), plot(L, ARI,'g-o')
    title('ARI')
    subplot(3,1,3), plot(L, RI,'b-o')
    title('RI')
end

