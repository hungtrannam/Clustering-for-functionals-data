function results = ERACF_(Data, param, varargin)


p = inputParser;
addParameter(p, 'Visualize', 'False');
addParameter(p, 'Saving', 'False');

parse(p, varargin{:});
Visualize = p.Results.Visualize;
Saving    = p.Results.Saving;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting
% rng(4);
f           = Data;
iter        = 0;
numSample   = size(f, 2);
epsilon     = param.epsilon;
l           = param.L;
max_iter    = param.maxIter;


if isfield(param, 'x')
    x = param.x;
end

% Calculate the distance between each pair of PDFs
for i = 1:numSample
    for j = 1:numSample
        if i ~= j
            Wf(i, j) = Integration(param.h, abs(f(:, i) - f(:, j)), 1) + 10^(-10); % L1  
        else
            % Distance between a PDF and itself should be a small positive number
            Wf(i, j) = 10^(-10);
        end
    end
end
Wf = Wf;

% mu
d      = Wf(:);
mu     = median(d);
sigma  = iqr(d)/2;

%f
alp    = ones(numSample,numSample);
alpnew = alp;

%x_(t+1)
for i=1:numSample
    tu  = 0;
    mau = 0;
    idx = [i randsample(numSample, floor(numSample/2))'];
    for j = idx
        if Wf(j,i) > mu*alp(j,i)
            U = 0;
        else
            U = exp(-Wf(j,i)/(sigma/l));
        end
        alpnew(j,i) = alp(j,i) / (1 + alp(j,i)*U);
        tu = tu + f(:,j)*U;
        mau = mau + U;
    end
    fnew(:,i) = tu./mau;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering

while iter < max_iter
    iter = iter + 1;

        % Calculate the distance between each pair of PDFs
    for i = 1:numSample
        for j = 1:numSample
            if i ~= j
                % L1 distance between PDFs i and j
                Wf(i, j) = Integration(param.h, abs(f(:, i) - f(:, j)), 1) + 10^(-10); % L1  
            else
                % Distance between a PDF and itself should be a small positive number
                Wf(i, j) = 10^(-10);
            end
        end
    end

    Wf     = Wf;

    % mu
    d      = Wf(:);
    mu     = median(d);
    sigma  = iqr(d)/2;

    %f
    alp    = ones(numSample,numSample);
    alpnew = alp;

    %x_(t+1)
    for i = 1:numSample
        tu  = 0;
        mau = 0;
        idx = [i randsample(numSample, floor(numSample/2))'];
        for j=idx
            if Wf(j,i)>mu*alp(j,i)
                U = 0;
            else
                U = exp(-Wf(j,i)/(sigma/l));
            end
            alpnew(j,i) = alp(j,i) / (1 + alp(j,i)*U);
            tu          = tu + f(:,j)*U;
            mau         = mau + U;
        end
        fnew(:,i)       = tu./mau;
    end


    if strcmp(Visualize, 'True')
        res = roundn(fnew',-1);
        idx = unique(res,'rows');
        IDX = zeros(size(res,1),1);
        for i = 1:size(IDX,1)
            for j = 1:size(idx,1)
                if res(i,:) == idx(j,:)
                    IDX(i) = j;
                end
            end
        end
        h = figure(1);
        PlotPDFeachIteration(Data, IDX, x);
        hold off;
        pause(0.01);
        if strcmp(Saving, 'True')
            fig_filename = fullfile('Figure_Output/CDE_cont/', ...
                sprintf('CDE_Plot(%03d).png', iter));
            saveas(h, fig_filename);
        end
    end


    % Check for convergence
    if norm(fnew-f,1) < epsilon
        break;
    end

    f = fnew;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results

res = roundn(fnew',-1);
idx = unique(res,'rows');
IDX = zeros(size(res,1),1);
for i = 1:size(IDX,1)
    for j = 1:size(idx,1)
        if res(i,:) == idx(j,:)
            IDX(i) = j;
        end
    end
end


fvnew = round(fnew, 4);
fv    = [];

for i = 1:numSample
    fvnewsort = fvnew(:, i);
    
    Overlap = false;
    for j = 1:size(fv, 2)
        if isequal(fvnewsort, fv(:, j))
            Overlap = true;
            break;
        end
    end
    
    % Nếu không trùng lặp, thêm cột vào B
    if ~Overlap
        fv = [fv fvnewsort];
    end
end

% results.Cluster.U = U;
results.Data.fvnew = fnew;
results.Data.fv = fv;

results.iter = iter;
results.Data.Data = Data;
results.Cluster.IDX = IDX';
results.Dist.D = Wf;


end
