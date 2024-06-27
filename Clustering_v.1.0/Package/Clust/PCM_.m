function results = PCM_(Data, param, varargin)

if param.kClust <= 1
    error('The number of clusters (kClust) must be greater than 1 for clustering to be performed.');
end

p = inputParser;
addParameter(p, 'Visualize', 'None');
parse(p, varargin{:});
Visualize = p.Results.Visualize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting
f = Data;
iter = 0;
max_iter = param.maxIter;
fm = param.mFuzzy;
epsilon = param.epsilon;
numSample = size(f, 2);
numCluster = param.kClust;
alpha = param.alphaCut;
if isfield(param, 'x')
    x = param.x;
end
K = param.K;
abnormal = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering
%% Initialize the partition matrix with FCM
if param.FvIni == 1
    fv = f(:, randperm(numSample, numCluster));
elseif param.FvIni == 2
    U0 = rand(numCluster, numSample);
    U0 = U0 ./ sum(U0);

    fv = (f * U0.^fm') ./ sum(U0.^fm, 2)'; 
elseif param.FvIni == 3

    fv = zeros(size(f, 1), numCluster);
    tt = randi([1, numSample], 1, 1);
    fv(:,1) = f(:,tt);

    min_distance_threshold = 0.000001; % Define a minimum distance threshold

    for i = 2:numCluster
        max_distance = -inf;
        best_tt = -1;
        attempts = 0;
        max_attempts = 1000;

        while attempts < max_attempts
            tt = randi([1, numSample], 1, 1);
            nearest = inf;
            valid = true;

            for j = 1:i-1
                distance = 2 - trapz(x, min(fv(:, j), f(:, tt)));
                if distance < min_distance_threshold
                    valid = false;
                    break;
                end
                nearest = min(nearest, distance);
            end

            if valid && nearest > max_distance
                max_distance = nearest;
                best_tt = tt;
            end

            attempts = attempts + 1;
        end

        if best_tt ~= -1
            fv(:, i) = f(:, best_tt);
        else
            error('Could not find a valid initial center. Try increasing max_attempts or adjusting min_distance_threshold.');
        end
    end

end

fv0 = fv;

figure;
plot(x, fv); 
hold off;

%% Repeat FCM until convergence or max iterations
while iter < max_iter
    iter = iter + 1;

    % Calculate the distance between fv with fi PDFs
    for j = 1:numSample
        for i = 1:numCluster
            Wf(i, j) = trapz(x, abs(fv0(:, i) - f(:, j))).^2 + 10^(-10); % L1  
            % Wf(i,j) = trapz(x, (fv(:, i) - f(:, j)).^2 ./ max(fv(:, i), f(:, j))) + 10^(-10); % Vicissitude chi^2
            % Wf(i, j) = trapz(x, fv(:,i) .* log(fv(:,i) ./ f(:,j))) + 10^(-10); % Kullback-Leibler
            % Wf(i, j) = (2 - trapz(x, max(fv(:, i), f(:, j))))^2 + 10^(-10); % Overlap
            % Wf(i, j) = 2 * (1 - (1/2) * trapz(param.x, max(fv(:,i), f(:,j)))) + 10^(-10); % SCC
            % Wf(i,j) = 1 - trapz(x, min(fv(:, i), f(:, j)))+ 10^(-10); % Intersection
            % Wf(i,j) = trapz(x, (fv(:, i) - f(:, j)).^2 ./ fv(:, i)).^2 + 10^(-10); % Neyman chi^2
        end
    end

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

    % Calculate the cluster centers
    fvnew = (f * Unew.^fm') ./ sum(Unew.^fm, 2)';

    % Check for convergence
    if norm(fv - fvnew, 1) < epsilon
        break;
    end

    fv = fvnew;


end


% Calculate the distance between fv with fi PDFs
for j = 1:numSample
    for i = 1:numCluster
        Wp(i, j) = trapz(x, abs(fv0(:, i) - f(:, j))).^2 + 10^(-10);
    end
end

% Estimate eta by FCM results
for i = 1:numCluster
    eta(i) = K * sum((Unew(i,:).^fm) .* Wp(i,:)) / sum(Unew(i,:).^fm);
end

iter = 0;
%% Repeat PCM until true condition nor max_ter
while iter < max_iter
    iter = iter + 1;

    % Calculate the distance between fv with fi PDFs
    for j = 1:numSample
        for i = 1:numCluster
            Wp(i, j) = trapz(x, abs(fv0(:, i) - f(:, j))).^2 + 10^(-10);
        end
    end


    % Update partition matrix
    for j = 1:numSample
        m = 0;
        for k = 1:numCluster
            if Wp(k, j) == 0
                m = m + 1;
            end
        end
        if m == 0
            Upcm = 1./(1 + (Wp./eta').^(1/(fm-1)));
        else
            for l = 1:numCluster
                if Wp(l, j) == 0
                    Upcm(l, j) = 1 / m;
                else
                    Upcm(l, j) = 0;
                end
            end
        end
    end
            
    % Update the representation PDF fv
    fvpnew = (f * (Upcm.^fm)') ./ sum(Upcm.^fm, 2)';

    % Calculate ObfFun by Krishnapuram 1993
    ObjFun = sum(sum(Upcm.^fm .* Wp)) + sum(eta).*sum((1 - Upcm).^fm, 'all');
    fprintf('Iteration count = %d, obj. pcm = %f\n', iter, ObjFun);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting each iteration

    if strcmp(Visualize, 'CDF')
        [~, IDX] = max(Upcm);
        for i = 1:numSample
            if all(Upcm(:, i) < alpha) == 1
                IDX(i) = 0;
                abnormal = [abnormal i];
            end
        end
        h = PlotPDFeachIteration(f, IDX, x);
        hold on;
        plot(x, fvpnew, "LineWidth", 2, "DisplayName", "representative PDFs");
        fig_filename = fullfile('Figure_Output/CDF_cont/', ...
            sprintf('CDF_Plot(%03d).fig', iter));
        saveas(h, fig_filename);
        pause(0.2); close;
    end


    if norm(fv0 - fvpnew, 1) < epsilon
        break
    end

    fv0 = fvpnew;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results
[~, IDX] = max(Upcm);
for i = 1:numSample
    if all(Upcm(:, i) < alpha) == 1
        IDX(i) = 0;
        abnormal = [abnormal i];
    end
end

results.Cluster.U = Upcm;
results.Data.fv = fvpnew;
results.iter = iter;
results.ObjFun = ObjFun;
results.Data.Data = f;
results.Cluster.IDX = IDX;
results.Dist.D = Wp;
results.isnoise = abnormal;

end
