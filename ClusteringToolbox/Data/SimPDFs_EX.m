clear; close all;

%%%%%%%%%%%%%%%%%%%%%
% CDF
% Data

mu_ranges = { ...
    linspace(0, 1, 5*100), ...
    % linspace(3, 6, 5), ...
    };

% abnormal_params = {
%     {[0.15, 2.5], sqrt([.4, .4])} 
%     {[0.3, 2.2], sqrt([.4, .4])}
%     };

[Data, param.x, param.truelabels] = SimPDFAbnormal( ...
    mu_ranges, ...
    sqrt([.5, .5, .5]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

figure;
hold on;
colors = parula(length(mu_ranges)); % Tạo bảng màu
for clust = 1:max(param.truelabels)
    clust_pdfs = Data(:, param.truelabels == clust);
    for i = 1:size(clust_pdfs, 2)
        plot(param.x, clust_pdfs(:,i), 'Color', colors(clust,:), 'LineWidth', 2);
    end
end
title(sprintf('Plotting of %d distributions', length(param.truelabels)));
xlabel('Value');
ylabel('Probability');
hold off;


INFO.SCC = (size(Data,2) / (size(Data,2) - 1)) * (1 - (1/size(Data,2)) *trapz(param.x, max(Data')'));
INFO.CW = (trapz(param.x, max(Data')')) - 1