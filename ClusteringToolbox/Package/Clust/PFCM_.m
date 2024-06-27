function results = PFCM_(Data, param, varargin)

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
x = param.x;
K = param.K;
cF = param.cF;
cP = param.cP;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering

%% Initialize the partition matrix with FCM (U, fv, num_cluster)
% load inima.mat
numCluster = 2;

ops = fcmOptions(...
    NumClusters=numCluster,...
    Exponent=fm,...
    MaxNumIteration=max_iter, ...
    Verbose = true);

[fv,Uf] = fcm(f', ops);

% fv = SUF
fv =fv';
% Calculate the distance between fv with fi PDFs
for j = 1:numSample
    for i = 1:numCluster
        Wf(i, j) = trapz(x, abs(fv(:, i) - f(:, j))).^2;
    end
end

% Estimate eta by FCM results
for i = 1:numCluster
    eta(i) = K * sum((Uf(i,:).^fm) .* Wf(i,:)) / sum(Uf(i,:).^fm);
end

%% Repeat FCM until convergence or max iterations
while iter < max_iter
    iter = iter + 1;

    % Calculate the distance between fv with fi PDFs
    for j = 1:numSample
        for i = 1:numCluster
            Wp(i, j) = trapz(x, abs(fv(:, i) - f(:, j))).^2;
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


    % Update partition matrix
    for j=1:numSample
        m=0;
        for k=1:numCluster
            if Wp(k,j)==0
                m=m+1;
            end
        end
        if m==0
            for i=1:numCluster
                tong=0;
                for k=1:numCluster
                    tong=tong+(Wp(i,j)^(2/(fm-1)))/(Wp(k,j)^(2/(fm-1)));
                end
                Ufcm(i,j)=1/tong;
            end
        else
            for l=1:numCluster
                if Wp(l,j)==0
                    Ufcm(l,j)=1/m;
                else
                    Ufcm(l,j)=0;
                end
            end
        end
    end

    ObjFun = sum(sum((cF.*Ufcm.^fm + cP.*Upcm.^fm) .* Wp)) + sum(eta).*sum((1 - Upcm).^fm, 'all');

    fv = (f * (cF.*Ufcm.^fm + cP.*Upcm.^fm)') ./ sum(cF.*Ufcm.^fm + cP.*Upcm.^fm, 2)';

    % Calculate the norm
    Cond = norm(Upcm - Uf, 1);
    fprintf('Iteration count = %d, obj. pcm = %f\n', iter, ObjFun);

    if Cond < epsilon
        break
    end
    Uf = Upcm;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results
[~,IDX] = max(Upcm);

results.Cluster.U = Upcm;
results.Data.fv = fv;
results.iter = iter;
results.ObjFun = ObjFun;
results.Data.Data = f;
results.Cluster.IDX = IDX;
results.Dist.D = Wf;

end
