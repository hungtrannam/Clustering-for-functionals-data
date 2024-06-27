function results = GKFCM_(Data, param, varargin)
    
    if param.kClust <= 1
        error('The number of clusters (kClust) must be greater than 1 for clustering to be performed.');
    end

    p = inputParser;
    addParameter(p, 'Visualize', 'None');
    parse(p, varargin{:});
    Visualize = p.Results.Visualize;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setting
    % rng(4);
    f = Data;
    iter = 0;
    max_iter = param.maxIter;
    fm = param.mFuzzy;
    epsilon = param.epsilon;
    numSample = size(f, 2);
    numCluster = param.kClust;
    if isfield(param, 'x')
        x = param.x;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Clustering 

    %% Initialize the partition matrix with FCM
    U = rand(numCluster, numSample);
    U = U ./ sum(U, 1);

    %% Repeat FCM until convergence or max iterations
    while iter < max_iter
        iter = iter + 1;

        % Calculate the cluster centers
        Um = U.^fm;
        fv = (f * Um') ./ sum(Um, 2)';

        % Calculate the distance between fv with fi PDFs
        Wf = zeros(numCluster, numSample);
        for i = 1:numCluster
            covMatrix = cov(f'); 
            invCovMatrix = inv(covMatrix); 
            for j = 1:numSample
                diff = fv(:, i) - f(:, j);
                Wf(i, j) = sqrt(diff' * invCovMatrix * diff);
            end
        end


        % Update partition matrix
        Unew = zeros(numCluster, numSample);
        for i = 1:numCluster
            for j = 1:numSample
                sumDist = 0;
                for k = 1:numCluster
                    sumDist = sumDist + (Wf(i,j)^(2/(fm-1))) / (Wf(k,j)^(2/(fm-1)));
                end
                Unew(i, j) = 1 / sumDist;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting each iteration

        if strcmp(Visualize, 'CDF')
            [~, labels] = max(Unew);
            if length(unique(labels)) < numCluster
                labels = randi([1, numCluster], 1, numSample);
            else
                [~, labels] = max(Unew);
            end
            h = PlotPDFeachIteration(f, labels, x);
            hold on;
            plot(x, fv, "LineWidth", 2, "DisplayName", "representative PDFs");
            fig_filename = fullfile('Figure_Output/CDF_cont/', ...
                sprintf('CDF_Plot(%03d).fig', iter));
            saveas(h, fig_filename);
            pause(0.2); close;
        elseif strcmp(Visualize, 'CDE')
            [~, labels] = max(Unew);
            h = PlotCDEseachIteration(f, fv',labels, param);
            fig_filename = fullfile('Figure_Output/CDE_cont/', ...
                sprintf('CDE_Plot(%03d).fig', iter));
            saveas(h, fig_filename);
            pause(0.2); close;
            
        end



        ObjFun = sum(sum((Unew) .* Wf));

        fprintf('Iteration count = %d, obj. fcm = %f\n', iter, ObjFun);

        % Check for convergence
        if norm(U - Unew, 1) < epsilon
            break;
        end

        U = Unew;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Results
    [~,IDX] = max(U);

    results.Cluster.U = U;
    results.Data.fv = fv;
    results.iter = iter;
    results.ObjFun = ObjFun;
    results.Data.Data = f;
    results.Cluster.IDX = IDX;
    results.Dist.D = Wf; 

end
