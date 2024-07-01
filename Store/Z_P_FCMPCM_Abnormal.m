close all; clear;
a = 0.01; 
ARI_N = [];F1 = [];time_N =[];iter_N=[];sen=[];Gmean=[];
% rng(3401);

Path_data = 'D:\FCM\Data\Images\TexturalImages\Clustering_Texture2';
% Data_annotations = 'D:\FCM\Data\Images\TexturalImages\PH2\PH2_lables.csv';

tic
% Truelables_table=readtable(5Data_annotations,'VariableNamingRule','preserve');
% Truelables = categorical(Truelables_table.Diagnosis);
% Truelables = [1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;3];

val=[];
for cv = 1:30
imds = imageDatastore(Path_data, ...
    'IncludeSubfolders', true, ...
    'FileExtensions', '.jpg', ...
    'LabelSource','foldernames');
% 'Labels', Truelables);

label   = imds.Labels;
numClasses = size(tabulate((imds.Labels)),1);
num_cluster = numClasses-1;

RotationRange = [-90 90];
% pixelRange = [-30 30];
imageAugmenter = imageDataAugmenter( ...
    'RandXReflection',true, ...
    'RandRotation',RotationRange);
% 'RandXTranslation',pixelRange, ...
% 'RandYTranslation',pixelRange, ...

for p = 1:length(imds)
    img = imread(imds.Files{p});
    augimds = augmentedImageDatastore(size(img, 1:2), imds);%, 'DataAugmentation', imageAugmenter);
end

pdf = ExtractPDFbyKernel(augimds, imds.Labels);

% ECG = readtable("ECG200_TRAIN.txt");
% pdf = table2array(ECG)';
% pdf = pdf(2:end,:);

for channels = 1:size(pdf, 2)

    fchannel = pdf{channels};
    f = fchannel(1:end-1,:);
    % Setup parameters
    iter = 0; Cond= Inf;
    max_iter = 300;
    fm = 2;
    epsilon = 0.00001;
    K=1;
    num_sample = size(f, 2);


    %% Initialize the partition matrix with FCM (U, fv, num_cluster)
    ops = fcmOptions(...
        NumClusters=num_cluster,...
        Exponent=fm,...
        MaxNumIteration=max_iter, ...
        DistanceMetric = 'euclidean');
    opt.Verbose = false;


    [fv,Uf] = fcm(f', ops);
    fvfcm = fv'; Ufcm = Uf;

    fv =fv';
    % Calculate the distance between fv with fi PDFs
    for j = 1:num_sample
        for i = 1:num_cluster
            Wf(i, j) = norm((fv(:, i) - f(:, j)),2).^2;
            % Wf(i, j) = dtw(fv(:, i),f(:, j));
            % Wf(i, j) = (1 - OverlapCoefficient(fv(:, i), f(:, j))).^2;
        end
    end

    % Estimate eta by FCM results
    for i = 1:num_cluster
        eta(i) = K * sum((Uf(i,:).^fm) .* Wf(i,:)) / sum(Uf(i,:).^fm);
    end

    %% Repeat PCM until true condition nor max_ter
    for e = 1:max_iter
        iter = iter + 1;

        % Calculate the distance between fv with fi PDFs
        for j = 1:num_sample
            for i = 1:num_cluster
                Wp(i, j) = norm((fv(:, i) - f(:, j)),2).^2;
                % Wp(i,j) = dtw(fv(:, i),f(:, j));
                % Wp(i, j) = (1 - OverlapCoefficient(fv(:, i), f(:, j))).^2;
            end
        end

        % Update partition matrix
        for j = 1:num_sample
            m = 0;
            for k = 1:num_cluster
                if Wp(k, j) == 0
                    m = m + 1;
                end
            end
            if m == 0
                Upcm = 1./(1 + (Wp./eta').^(1/(fm-1)));
            else
                for l = 1:num_cluster
                    if Wp(l, j) == 0
                        Upcm(l, j) = 1 / m;
                    else
                        Upcm(l, j) = 0;
                    end
                end
            end
        end

        % Calculate ObfFun by Krishnapuram 1993
        ObjFun(e) = sum(sum(Upcm.^fm .* Wp)) + sum(eta).*sum((1 - Upcm).^fm, 'all');

        % Update the representation PDF fv
        fv = (f * (Upcm.^fm)') ./ sum(Upcm.^fm, 2)';

        % scatter(f(1,:)', f(2,:)');
        % hold on
        % scatter(fv(1,:)',fv(2,:)',80,"filled");
        % hold off
        % pause(.5)
        %
        % figure
        % plot(f,'k-')
        % hold on
        % plot(fv,'LineWidth',2)
        % hold off

        % Calculate the norm
        Cond(e) = norm(Upcm - Uf, 1);
        fprintf('Iteration count = %d, obj. pcm = %f\n', e, ObjFun(e));

        if Cond(e) < epsilon
            break
        end
        Uf = Upcm;

    end

    % figure
    % subplot(4,1,1)
    % heatmap(Upcm);
    % subplot(4,1,2)
    % heatmap(Ufcm);


    memship{channels} = Upcm;
end

result.data.f=Upcm;
result.data.d=sqrt(Wp);
result.cluster.v=fv;
result.iter = iter;
result.cost = ObjFun;
result.param.fm = fm;

% threshold for dectection noising and outliers.
DectectNoise = zeros(1,num_sample);
for i = 1:num_sample
    if all(Upcm(:, i) <= a)
        DectectNoise(i) = 1;
    end
end
[~,NoiseIDX] = find(DectectNoise==1);

% XB  = validity(result)

[~,idx] = max(Ufcm);
truelabels = pdf{1,1}(end,:);
[RI, ARI] = randindex(idx, truelabels);


truelabels_binary = [zeros(1,size(pdf{1,1},2)-15) ones(1,15)];
idx_binary = zeros(size(truelabels_binary));
for i = 1:length(NoiseIDX)
    idx_binary(NoiseIDX(i)) = 1;
end

cfs = confusionmat(truelabels_binary, idx_binary);

performPCM = performance(cfs,1);


ARI_N = [ARI_N ARI];
F1 = [F1 performPCM.classes(9,2)]
sen = [sen performPCM.classes(6,2)];
time_N = [time_N toc];
iter_N = [iter_N iter];

Gmean = [Gmean performPCM.classes(6,2)*performPCM.classes(7,2)];
end


mean(ARI_N)
mean(F1)
mean(Gmean)
mean(sen)
mean(time_N)
mean(iter_N)