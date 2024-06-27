clear; close all;


addpath /home/hung-tran-nam/Hung_matlab/Clustering_v.1.0/Package/ToolBox
addpath /home/hung-tran-nam//Hung_matlab/Clustering_v.1.0/Package/Clust/
addpath /home//hung-tran-nam//Hung_matlab/Clustering_v.1.0/Package/Vis/
addpath /home//hung-tran-nam//Hung_matlab/Clustering_v.1.0/Package/Val/

addpath /home/hung-tran-nam/Hung_matlab/Data


% Path_data = '/home/hung-tran-nam/Hung_matlab/Data/IMBA1';
Path_data = '/home/hung-tran-nam/Hung_matlab/Data/Brodatz/cropped_images';

param.maxIter    = 1000;
param.mFuzzy     = 2;
param.epsilon    = 1e-10;
param.val        = 2;
param.FvIni      = 3;
param.Klist      = 2:5;

imds = imageDatastore(Path_data, ...
    'IncludeSubfolders', true, ...
    'FileExtensions', '.gif', ...
    'LabelSource','foldernames');

rng('default') % For reproducibility

numClasses            = length(unique(imds.Labels));
kfold                 = 5;
fold                  = cvpartition(imds.Labels,'kfold',kfold,'Stratify',true);
Afold                 = zeros(kfold,1);

% RotationRange = [-90 90];
% imageAugmenter = imageDataAugmenter( ...
%     'RandXReflection',true, ...
%     'RandRotation',RotationRange);

for i = 1:1
    Train_idx    = fold.training(i);
    Test_idx     = fold.test(i);
    imdsTrain    = imageDatastore(imds.Files(Train_idx),"Labels", imds.Labels(Train_idx));
    imdsTest     = imageDatastore(imds.Files(Test_idx),"Labels", imds.Labels(Test_idx));

    % [imdsVal, ~] = splitEachLabel(imdsTrain, .2);

    % for p = 1:length(imdsTrain)
    %     img = imread(imdsTrain.Files{p});
    %     augimdsTrain = augmentedImageDatastore(size(img, 1:2), imdsTrain, 'DataAugmentation', imageAugmenter);
    % end

    [pdfTrain, param.x] = ExtractKernel(imdsTrain);
    [pdfTest,  ~]       = ExtractKernel(imdsTest);
    % [pdfTest, param.x]   = ExtractKernel(imdsTest);
    % [pdfTest, param.x]   = ExtractI2PDF(imdsVal);


    bestK = evalclusters_(pdfTrain, param, param.Klist);
    param.kClust = bestK;
    training = IFCM_(pdfTrain, param);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting
    % figure;
    % heatmap(training.Cluster.U);

    figure;
    PlotPDFeachIteration(pdfTrain, training.Cluster.IDX, param.x); hold on;
    plot(param.x,training.Data.fv, "LineWidth", 3, "DisplayName", "Representative PDF"); hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Validation
    % [param.truelabels, ~, label_index] = grp2idx(imdsTrain.Labels);
    % param.truelabels                   = param.truelabels';
    % training                           = validityClustering(training, param);

    [param.truelabels, ~, label_index] = grp2idx(imdsTest.Labels);
    param.truelabels                   = param.truelabels';

    testing                            = predictClust_(pdfTest, training, param);
    testing                            = validityClustering(testing, param);

    Afold(i,:)                         = testing.validity.Gmean;

end

mean(Afold)
