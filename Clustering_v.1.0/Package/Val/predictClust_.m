function results = predictClust_(Data, training, param, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting
% rng(4);
f = Data;
fm = param.mFuzzy;
numSample = size(f, 2);
numCluster = param.kClust;

if isfield(param, 'x')
    x = param.x;
end

% from training

fv = training.Data.fv;

if isfield(param, 'Imba')
    fci = training.Imba;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering

% Calculate the distance between fv with fi PDFs
for j = 1:numSample
    for i = 1:numCluster
        Wf(i, j) = trapz(x, abs(fv(:, i) - f(:, j))) + 10^(-10); % L1
        % Wf(i,j) = trapz(x, (fv(:, i) - f(:, j)).^2 ./ max(fv(:, i), f(:, j))) + 10^(-10); % Vicissitude chi^2
        % Wf(i, j) = trapz(x, fv(:,i) .* log(fv(:,i) ./ f(:,j))) + 10^(-10); % Kullback-Leibler
        % Wf(i, j) = (2 - trapz(x, max(fv(:, i), f(:, j))))^2 + 10^(-10); % Overlap
        % Wf(i, j) = 2 * (1 - (1/2) * trapz(param.x, max(fv(:,i), f(:,j)))) + 10^(-10); % SCC
        % Wf(i,j) = 1 - trapz(x, min(fv(:, i), f(:, j)))+ 10^(-10); % Intersection
        % Wf(i,j) = trapz(x, (fv(:, i) - f(:, j)).^2 ./ fv(:, i)).^2 + 10^(-10); % Neyman chi^2
    end
end


if isfield(param, 'Imba') %IFCM
    %Update partition matrix
    for i = 1:numCluster
        for j = 1:numSample
            numerator = fci(i) / ((Wf(i, j)) ^ (2/(fm - 1)));
            denominator = 0;
            for k = 1:numCluster
                denominator = denominator + fci(k) / ((Wf(k, j)) ^ (2/(fm - 1)));
            end
            Unew(i, j) = numerator / denominator;
        end
    end
else %FCM

    % Update partition matrix
    Unew = ones(numCluster, numSample) / numCluster;
    for i = 1:numCluster
        for j = 1:numSample
            sumDist = 0;
            for k = 1:numCluster
                sumDist = sumDist + (Wf(i,j)^(2/(fm-1))) / (Wf(k,j)^(2/(fm-1)));
            end
            Unew(i, j) = 1 / sumDist;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results
[~,IDX] = max(Unew);

results.Cluster.U = Unew;
results.Cluster.IDX = IDX;
results.Dist.D = Wf;
results.Data.Data = f;

end
