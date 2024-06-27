close all; clear;

%
% pdf = [60 150;65 150;70 150;55 150;50 150;60 145;60 155; 140 150;145 150;150 150; 135 150;130 150;140 145;140 155;100 150;100 200]';
% pdf = [30 80 140 150;40 50 160 170;60 70 130 190;80 90 170 180;50 90 120 160;80 100 110 130;180 220 160 190;190 200 110 180;180 220 140 150;210 230 130 150;210 240 170 180;240 250 180 190;10 20 10 20;10 20 280 290;270 300 270 300;280 300 10 20]';
% pdf = [-5 0; -3.34 1.67; -3.34 0; -3.34 -1.67; -1.67 0; 0 0; 1.67 0; 3.34 0; 3.34 1.67; 3.34 -1.67; 5 0; 0 27]';

f1 =[]; f2 = []; f3 = [];f4 = [];f5 = [];
x = linspace(-0.2, 1.5, 1000);

mu_1 = 0.25:.001:0.35;
mu_2 = 0.75:.001:0.8;
mu_3 = 0.55:0.01:0.56;
sigma_sq_3 = 0.01;

mu = size(mu_1,2) + size(mu_2,2);
sigma_sq = .01;
for i = 1:length(mu_1)
    f = normpdf(x, mu_1(i), sqrt(sigma_sq));
    f1 = [f1; f 1];
end
for i = 1:length(mu_2)
    f = normpdf(x, mu_2(i), sqrt(sigma_sq));
    f2 = [f2; f 2];
end

for i = 1:length(mu_3)
    f = normpdf(x, mu_3(i), sqrt(sigma_sq_3));
    f3 = [f3; f 3];
end

pdf = [f1; f2; f3]';

% ECG = readtable("ECG200_TRAIN.txt");
% pdf = table2array(ECG)';
% pdf = pdf(2:end,:);

% ECG = readtable("ECG5000_TRAIN.txt");
% pdf = table2array(ECG)';
% pdf = pdf(2:end,:);

% load D:\FCM\Data\pdf.mat
% for i = 1:130
%     kd(i,:) = ksdensity(f(i,:));
% end
% fa = kd(103:115,:)';
% fb = kd(81:102,:)';

f = pdf(1:end-1,:);

% Setup parameters
iter = 0; Cond= Inf;
max_iter = 300;
fm =2;
epsilon = 0.00001;
K=1;
num_sample = size(f, 2);
num_cluster =2;


%% Initialize the partition matrix with FCM (U, fv, num_cluster)
ops = fcmOptions(...
    NumClusters=num_cluster,...
    Exponent=fm,...
    MaxNumIteration=max_iter, ...
    DistanceMetric = 'euclidean');
tic
[fv,Uf] = fcm(f', ops);
tFCM = toc

fvfcm = fv'; Ufcm = Uf;

tic
fv =fv';
% Calculate the distance between fv with fi PDFs
for j = 1:num_sample
    for i = 1:num_cluster
        Wf(i, j) = norm((fv(:, i) - f(:, j)),2).^2;
        % Wf(i, j) = dtw(fv(:, i),f(:, j)).^2;
        % Wf(i, j) = mahal(fv(:, i),f(:, j)).^2;
        Wf(i, j) = (OverlapCoefficient(fv(:, i), f(:, j))).^2;
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
            % Wp(i, j) = norm((fv(:, i) - f(:, j)),2).^2;
            % Wp(i,j) = dtw(fv(:, i),f(:, j)).^2;
            % Wp(i, j) = mahal(fv(:, i),f(:, j)).^2;
            Wp(i, j) = (OverlapCoefficient(fv(:, i), f(:, j))).^2;
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

    % Calculate the norm
    Cond(e) = norm(Upcm - Uf, 1);
    fprintf('Iteration count = %d, obj. pcm = %f\n', e, ObjFun(e));

    if Cond(e) < epsilon
        break
    end
    Uf = Upcm;

end
tPCM = toc

figure
subplot(4,1,1)
heatmap(Upcm);
subplot(4,1,2)
heatmap(Ufcm);
subplot(4,1,3)
plot(f,'g-.')
hold on
plot(fv,'LineWidth',3)
hold off
subplot(4,1,4)
plot(f,'g-.')
hold on
plot(fvfcm,'LineWidth',3)
hold off

% for i = 1:num_cluster
%     for j=1:num_sample
%         rectangle('Position',[f(1,j) f(3,j)  f(2,j)-f(1,j)  f(4,j)-f(3,j)],'EdgeColor','k')
%         rectangle('Position',[fv(1,i) fv(3,i)  fv(2,i)-fv(1,i)  fv(4,i)-fv(3,i)],'LineWidth',3, 'EdgeColor', 'b')
%         rectangle('Position',[fvfcm(1,i) fvfcm(3,i)  fvfcm(2,i)-fvfcm(1,i)  fvfcm(4,i)-fvfcm(3,i)],'LineWidth',3,'EdgeColor', 'g')
%     end
% end
% axis([5 300 0 301]);
% hold off


% for k = 1:100
%     a = 0 + 0.005 * (k - 1); % Tính giá trị a tương ứng
%     DectectNoise = zeros(1, num_sample);
% 
%     % Tính DectectNoise cho giá trị a hiện tại
%     for i = 1:num_sample
%         if all(Upcm(:, i) <= a)
%             DectectNoise(i) = 1;
%         end
%     end
% 
%     % Lưu kết quả vào mảng results
%     Noise(k, :) = DectectNoise;
% end

% threshold for dectection noising and outliers.
a = 0.1; DectectNoise = zeros(1,num_sample);
for i = 1:num_sample
    if all(Upcm(:, i) <= a)
        DectectNoise(i) = 1;
    end
end
[~,NoiseIDX] = find(DectectNoise==1);


result.data.f=Upcm;
result.data.d=sqrt(Wp);
result.cluster.v=fv;
result.iter = iter;
result.cost = ObjFun;
result.param.fm = fm;

% XB  = validity(result);

[~,idx] = max(Ufcm);
truelabels = pdf(end,:);
% [RI, ARI] = randindex(idx, truelabels)


% truelabels_binary = [zeros(1,size(f1,1) + size(f2,1) + size(f3,1)) ones(1,size(f4,1)) ones(1,size(f5,1)) ones(1,size(f6,1))];
% idx_binary = zeros(size(truelabels_binary));
% for i = 1:length(NoiseIDX)
%     idx_binary(NoiseIDX(i)) = 1;
% end
% 
% cfs = confusionmat(truelabels_binary, idx_binary);
% 
% performPCM = performance(cfs,1)