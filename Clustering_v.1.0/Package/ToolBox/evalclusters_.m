function [bestK, silhouetteValues] = evalclusters_(Data, param, Klist)

% Initialize variables to store the results
numK = length(Klist);
silhouetteValues = zeros(numK, 1);

% Loop over each value of K in Klist
for k = 1:numK
    param.kClust = Klist(k);

    % Perform k-means clustering
    IDXFCM = H_(Data, param);

    % Compute the silhouette values
    s = silhouette_(Data, IDXFCM.Cluster.IDX, param);

    % Store the mean silhouette value for the current K
    silhouetteValues(k) = s;

end

% Find the best K that maximizes the mean silhouette value
[~, bestIdx] = max(silhouetteValues);
bestK = Klist(bestIdx);

end

