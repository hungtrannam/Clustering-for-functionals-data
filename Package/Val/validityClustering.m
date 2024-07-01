function results = validityClustering(results, param)
% validation of the clustering
N = size(results.Data.Data, 2);

if param.val == 1 || param.val == 3
    m = param.mFuzzy;

    % partition coefficient (PC)
    fm = (results.Cluster.U).^m;
    PC = 1/N * sum(sum(fm));
    % classification entropy (CE)
    fm = (results.Cluster.U) .* log(results.Cluster.U);
    CE = -1/N * sum(sum(fm));

    % Xie and Beni's index (XB)
    XB = sum((sum(results.Dist.D .* results.Cluster.U.^2)) ./ (N * min(results.Dist.D)));

    % Silhouette Coefficient
    silhouetteVals = silhouette_(results.Data.Data, results.Cluster.IDX, param);
    SC = mean(silhouetteVals);

    results.validity.SC = SC;
    results.validity.XB = XB;
    results.validity.PC = PC;
    results.validity.CE = CE;
end

if param.val == 2 || param.val == 3
    [RI, ARI] = randindex(param.truelabels, results.Cluster.IDX);
    Gmean     = g_mean(param.truelabels, results.Cluster.IDX);


    results.validity.Gmean = Gmean;
    results.validity.ARI = ARI;
    results.validity.RI = RI;
end





end

