clear; close all;


addpath /home/hung-tran-nam/Hung_matlab/Clustering_v.1.0/Package/ToolBox
addpath /home/hung-tran-nam/Hung_matlab/Data


% I = [1 1 0 0 0;
%      1 1 0 0 0;
%      0 0 0 0 1];
%
% h = (4/3)^(1/5) * length(I(:)).^(-1/5) * std(I(:));
%
% [ker, x] = ksdensity(I(:), "Kernel","normal", NumPoints=1000);
% ker2 = ksdensity(I(:),"Bandwidth", h, "Kernel","normal", NumPoints=1000);
%
% plot(x, ker', 'r'); hold on;
% plot(x, ker2', 'b');


Path_data = '/home/hung-tran-nam/Hung_matlab/Data/IMBA1';

param.maxIter = 1000;
param.mFuzzy = 2;
param.epsilon = 1e-10;
param.kClust = 3;
param.val = 2;

imds = imageDatastore(Path_data, ...
    'IncludeSubfolders', true, ...
    'FileExtensions', '.png', ...
    'LabelSource','foldernames');

numClasses = length(unique(imds.Labels));
param.truelabels      = imds.Labels;
kfold      = 5;
fold       = cvpartition(param.truelabels,'kfold',kfold,'Stratify',true);
Afold      = zeros(kfold,1);
confmat    = 0;

RotationRange = [-90 90];
imageAugmenter = imageDataAugmenter( ...
    'RandXReflection',true, ...
    'RandRotation',RotationRange);

for i = 1
    Train_idx    = fold.training(i);
    Test_idx     = fold.test(i);
    imdsTrain    = imageDatastore(imds.Files(Train_idx),"Labels",param.truelabels(Train_idx));
    imdsTest     = imageDatastore(imds.Files(Test_idx),"Labels",param.truelabels(Test_idx));
    % [imdsVal, ~] = splitEachLabel(imdsTrain, .2);

    for p = 1:length(imdsTrain)
        img = imread(imdsTrain.Files{p});
        augimdsTrain = augmentedImageDatastore(size(img, 1:2), imdsTrain, 'DataAugmentation', imageAugmenter);
    end

    [pdfTrain,param.x]   = ExtractKernel(imdsTrain);
    % [pdfTest, x]   = ExtractKernel(imdsTest);
    % [pdfTest,x]   = ExtractI2PDF(imdsVal);


    results = FCM_(pdfTrain,param);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting
    % figure;
    heatmap(results.Cluster.U);
    
    
    h = PlotPDFeachIteration(pdfTrain, results.Cluster.IDX, param.x); hold on;
    plot(param.x,results.Data.fv, "LineWidth", 3, "DisplayName", "Representative PDF"); hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Validation
    % results = validityClustering(results, param);
    
end