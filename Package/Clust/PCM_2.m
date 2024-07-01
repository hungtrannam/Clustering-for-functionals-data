function results = PCM_2(Data, param, varargin)

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
% Clustering PCM



%% Initialize the partition matrix with FCM (U, fv, num_cluster)
options = fcmOptions(...
    NumClusters=numCluster,...
    Exponent=fm,...
    MaxNumIteration=max_iter, ...
    DistanceMetric = 'euclidean', ...
    Verbose = false);

[fv, U] = fcm(f', options);


fv = fv';
% Calculate the distance between fv with fi PDFs
for j = 1:numSample
    for i = 1:numCluster
        Wf(i, j) = trapz(x, abs(fv(:, i) - f(:, j))).^2;
    end
end

% Estimate eta by FCM results
for i = 1:numCluster
    eta(i) = K * sum((U(i,:).^fm) .* Wf(i,:)) / sum(U(i,:).^fm);
end

%% Repeat PCM until true condition nor max_ter
while iter < max_iter
    iter = iter + 1;

    % Calculate the distance between fv with fi PDFs
    for j = 1:numSample
        for i = 1:numCluster
            Wp(i, j) = trapz(x, abs(fv(:, i) - f(:, j))).^2;
            % Wp(i,j) = dtw(fv(:, i),f(:, j));
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

    % Calculate ObfFun by Krishnapuram 1993
    ObjFun = sum(sum(Upcm.^fm .* Wp)) + sum(eta).*sum((1 - Upcm).^fm, 'all');

    % Update the representation PDF fv
    fv = (f * (Upcm.^fm)') ./ sum(Upcm.^fm, 2)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting each iteration

    [~, IDX] = max(Upcm);
    for i = 1:numSample
        if all(Upcm(:, i) < alpha) == 1
            IDX(i) = 0;
            abnormal = [abnormal i];
        end
    end

    if strcmp(Visualize, 'CDF')
        h = PlotPDFeachIteration(f, IDX+1, x);
        hold on;
        plot(x, fv, "LineWidth", 2, "DisplayName", "representative PDFs");
        fig_filename = fullfile('Figure_Output/CDF_cont/', ...
            sprintf('CDF_Plot(%03d).fig', iter));
        saveas(h, fig_filename);
        pause(0.2); close;
    elseif strcmp(Visualize, 'CDE')
        h = PlotCDEseachIteration(f, fv',IDX, param);
        fig_filename = fullfile('Figure_Output/CDE_cont/', ...
            sprintf('CDE_Plot(%03d).fig', iter));
        saveas(h, fig_filename);
        pause(0.2); close;

    end

    Cond = norm(Upcm - U, 1);
    fprintf('Iteration count = %d, obj. pcm = %f\n', iter, ObjFun);

    if Cond < epsilon
        break
    end
    U = Upcm;

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
results.Data.fv = fv;
results.iter = iter;
results.ObjFun = ObjFun;
results.Data.Data = f;
results.Cluster.IDX = IDX;
results.Dist.D = Wf;
results.isnoise = abnormal;


end
