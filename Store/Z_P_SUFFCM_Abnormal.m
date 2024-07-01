close all; clear;
a = 0.05; 
ARI_N = [];F1 = [];time_N =[];iter_N=[];sen=[];Gmean=[];


Path_data = 'D:\FCM\Data\Images\TexturalImages\Clustering_Face';
% Data_annotations = 'D:\FCM\Data\Images\TexturalImages\PH2\PH2_lables.csv';

tic
% Truelables_table=readtable(5Data_annotations,'VariableNamingRule','preserve');
% Truelables = categorical(Truelables_table.Diagnosis);
% Truelables = [1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;3];

val=[];
imds = imageDatastore(Path_data, ...
    'IncludeSubfolders', true, ...
    'FileExtensions', '.png', ...
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
    augimds = augmentedImageDatastore(size(img, 1:2), imds, 'DataAugmentation', imageAugmenter);
end

pdf = ExtractPDFbyNorm(augimds, imds.Labels);


f = pdf{1,1}(1:end-1,:);
epsilon = 0.00001;
lamda = 0.5;
num_sample = size(pdf{1,1}, 2);
tic
%Lap ma tran khoang cach
W=[];
for i=1:num_sample-1
    for j=i+1:num_sample
        W(i,j)=norm((f(:, i)-f(:, j)),2);
        W(j,i)=W(i,j);
    end
end

%Tinh ds
Ws=sum(sum(W))/(length(W)*length(W)-length(W));
%Tinh ma tran K lamda
for i=1:num_sample
    for j=1:num_sample
        if W(i,j)<Ws
            k(i,j)=exp(-W(i,j)/(Ws/lamda));
        else
            k(i,j)=0;
        end
    end
end

%Tinh z(t+1)
fv=[];
for i=1:length(W)
    tu=[];
    for l=1:size(f,1)
        tu=[tu;0];
    end
    mau=0;
    for j=1:length(W)
        tu=tu+f(:,j)*k(i,j);
        mau=mau+k(i,j);
    end
    fv=[fv tu/mau];
end
%so sanh chuan va do lech toi da cho phep

cond=[];
for i=1:size(fv,2)
    cond=[cond norm(f(:,i)-fv(:,i),2)];
end
cond=max(cond);


iter=0;
while cond>epsilon
    iter=iter+1;
    f=fv;

    %Lap ma tran khoang cach
    W=[];
    for i=1:size(f,2)-1
        for j=i+1:size(f,2)
            W(i,j)=norm((f(:, i) - f(:, j)),2);
            W(j,i)=W(i,j);
        end
    end
    %Tinh ds
    Ws=sum(sum(W))/(length(W)*length(W)-length(W));
    %Tinh ma tran K lamda
    for i=1:length(W)
        for j=1:length(W)
            if W(i,j)<Ws
                k(i,j)=exp(-W(i,j)/(Ws/lamda));
            else
                k(i,j)=0;
            end
        end
    end
    %Tinh f(t+1)
    fv=[];
    for i=1:num_sample
        tu=[];
        for l=1:size(f,1)
            tu=[tu;0];
        end
        mau=0;
        for j=1:length(W)
            tu=tu+f(:,j)*k(i,j);
            mau=mau+k(i,j);
        end
        fv=[fv tu/mau];
    end
    
    cond=[];
    for i=1:size(fv,2)
        cond=[cond norm(f(:,i) -fv(:,i),2)];
    end
    cond=max(cond);
    fprintf('Iteration count = %d, Ws = %f\n', iter, Ws);

end
toc

hardSUF=fv';
hardSUF=roundn(hardSUF,-5);
labels=unique(hardSUF,'rows');
IDX=zeros(size(hardSUF,1),1);
for i=1:size(IDX,1)
    for j=1:size(labels,1)
        if norm(hardSUF(i,:)-labels(j,:))==0
            IDX(i)=j;
        end
    end
end
inima=zeros(size(labels,1),size(IDX,1));
for i=1:size(inima,2)
    inima(IDX(i),i)=1;
end

fm = 2;
max_iter = 300;
ops = fcmOptions(...
    NumClusters=num_cluster,...
    Exponent=fm,...
    MaxNumIteration=max_iter);

[~,Uf] = fcm(f', ops);

Usuf = Uf;
fv = hardSUF';