%calling function to search the optimal number of clusters
close all
clear all

addpath /home/hung-tran-nam/Hung_matlab/Clustering_v.1.0/Data
addpath /home/hung-tran-nam/Hung_matlab/Clustering_v.1.0/Package/Clust
addpath /home/hung-tran-nam/Hung_matlab/Clustering_v.1.0/Package/Val

%%%%%%%%%%%%%%%%%%%%%%%
% Setting

% Data
[Data, param.x, param.truelabels] = SimPDFAbnormal( ...
    { ...
    linspace(0, 1, 10*100), ...
    linspace(3, 4, 10)}, ...
    sqrt([.5, .5]));


% FCM
param.maxIter    = 1000;
param.mFuzzy     = 2;
param.epsilon    = 1e-10;
param.kClust     = 2;
param.val        = 2;
param.FvIni      = 3;

ment             = [];
nkmax            = 5;


%%%%%%%%%%%%%%%%%%%%%%%
% Loot

index = 1;

for cln = 2:.5:nkmax
    param.mFuzzy  = cln;
    results       = IFCM_(Data,param);
    results       = validityClustering(results,param);
    % validation
    ment{index}   = results.validity;
    index          = index + 1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%
% Plotting



if param.val == 1 || param.val == 3
    PC=[]; CE=[]; SC=[]; S=[]; XB=[]; DI=[]; ADI=[];
    for i = 2:nkmax
        PC = [PC ment{i}.PC];
        CE = [CE ment{i}.CE];
        SC = [SC ment{i}.SC];
        XB = [XB ment{i}.XB];
    end
    figure(1)
    clf
    subplot(4,1,1), plot(2:nkmax,PC, 'k-o')
    title('Partition Coefficient (PC)')
    subplot(4,1,2), plot(2:nkmax,CE,'r-o')
    title('Classification Entropy (CE)')
    subplot(4,1,3), plot(2:nkmax,SC,'g-o')
    title('Separation Index (S)')
    subplot(4,1,4), plot(2:nkmax,XB, 'b-o')
    title('Xie and Beni Index (XB)')

end

if param.val == 2 || param.val == 3
    Gmean=[]; ARI=[]; RI=[];

    for i = 2:nkmax
        Gmean = [Gmean ment{i}.Gmean];
        ARI   = [ARI ment{i}.ARI];
        RI    = [RI ment{i}.RI];

    end
    figure(2)
    clf
    subplot(3,1,1), plot(2:nkmax, Gmean, 'r-o')
    title('Gmean')
    subplot(3,1,2), plot(2:nkmax, ARI,'g-o')
    title('ARI')
    subplot(3,1,3), plot(2:nkmax, RI,'b-o')
    title('RI')
end
