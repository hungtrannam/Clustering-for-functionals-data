function results = FCM_(Data, param, varargin)

if param.kClust <= 1
    error('The number of clusters (kClust) must be greater than 1 for clustering to be performed.');
end

p = inputParser;
addParameter(p, 'Visualize', 'False');
addParameter(p, 'Saving', 'False');

parse(p, varargin{:});
Visualize = p.Results.Visualize;
Saving    = p.Results.Saving;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting
% rng(4);
f            = Data;
iter         = 0;
max_iter     = param.maxIter;
fm           = param.mFuzzy;
epsilon      = param.epsilon;
numSample    = size(f, 2);
numCluster   = param.kClust;
if isfield(param, 'x')
    x = param.x;
end


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

figure;
plot(x, fv);
hold off;


%% Repeat FCM until convergence or max iterations
while iter < max_iter
    iter = iter + 1;

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting each iteration

    if strcmp(Visualize, 'True')
        [~, labels] = max(Unew);
        h = figure(1);
        subplot(1,2,1);
        PlotPDFeachIteration(f, labels, x); hold on
        plot(x, fv, "LineWidth", 2, "DisplayName", "representative PDFs");
        hold off;
        pause(0.01);
        if strcmp(Saving, 'True')
            fig_filename = fullfile('Figure_Output/CDE_cont/', ...
                sprintf('CDE_Plot(%03d).png', iter));
            saveas(h, fig_filename);
        end
        subplot(1,2,2);
        PlotCDEsContour(fvnew, labels, param, iter);
    end



    ObjFun = sum(sum((Unew) .* Wf));
    fprintf('Iteration count = %d, obj. fcm = %f\n', iter, ObjFun);

    % Check for convergence
    if norm(fv - fvnew, 1) < epsilon
        break;
    end

    fv = fvnew;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results
[~,IDX] = max(Unew);

results.Cluster.U = Unew;
results.Data.fv = fv;
results.iter = iter;
results.ObjFun = ObjFun;
results.Data.Data = f;
results.Cluster.IDX = IDX;
results.Dist.D = Wf;

end
