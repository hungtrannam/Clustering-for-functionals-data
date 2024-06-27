%calling function to search the optimal number of clusters
close all
clear all


%%%%%%%%%%%%%%%%%%%%%%%
% Setting

[Data, param.x, param.truelabels] = SimPDFAbnormal( ...
    { ...
    linspace(0, 1, 10*100), ...
    linspace(3, 4, 10)}, ...
    sqrt([.5, .5]));


% FCM
param.L = 0.1;
param.epsilon = 1e-5;
param.val = 2;
param.maxIter = 500;


ment=[];
nkmax = 30;
index = 1;


%%%%%%%%%%%%%%%%%%%%%%%
% Loot

for cln = 0.1:0.1:5
    param.L = cln;
    results = ERACF_(Data, param);
    results = validityClustering(results, param);
    % validation
    ment{index} = results.validity;
    index = index + 1; % Increment the index
end

%%%%%%%%%%%%%%%%%%%%%%%
% Plotting


if param.val == 1 || param.val == 3
    PC=[];CE=[];SC=[];S=[];XB=[];DI=[];ADI=[];
    for i=2:nkmax
        PC=[PC ment{i}.PC];
        CE=[CE ment{i}.CE];
        SC=[SC ment{i}.SC];
        XB=[XB ment{i}.XB];
    end
    figure(1)
    clf
    subplot(4,1,1), plot(2:nkmax,PC)
    title('Partition Coefficient (PC)')
    subplot(4,1,2), plot(2:nkmax,CE,'r')
    title('Classification Entropy (CE)')
    subplot(4,1,3), plot(2:nkmax,SC,'g')
    title('Separation Index (S)')
    subplot(4,1,4), plot(2:nkmax,XB)
    title('Xie and Beni Index (XB)')

end

if param.val == 2 || param.val == 3
    Gmean=[];ARI=[];RI=[];

    for i = 2:nkmax
        Gmean = [Gmean ment{i}.Gmean];
        ARI = [ARI ment{i}.ARI];
        RI = [RI ment{i}.RI];

    end
    figure(2)
    clf
    subplot(3,1,1), plot(2:nkmax,Gmean, 'r')
    title('Gmean')
    subplot(3,1,2), plot(2:nkmax,ARI,'r')
    title('ARI')
    subplot(3,1,3), plot(2:nkmax,RI,'r')
    title('RI')
end

h = PlotPDFeachIteration(Data, results.Cluster.IDX, param.x); legend("hide"); hold on;
plot(param.x, results.Data.fv, "LineWidth", 3, "DisplayName", "Representative PDF"); hold off;