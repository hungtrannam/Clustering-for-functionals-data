function [Data, labels] = SimPDFAbnormal(mu_ranges, sig_values, grid, abnormal_params)
% SimPDFAbnormal generates and plots PDFs for multiple sets of distributions
% with specified means and standard deviations, including optional
% "abnormal" distributions.
%
% Parameters:
% mu_ranges      : Cell array of ranges of means for each set of distributions
% sig_values     : Array of standard deviations for all distributions
% abnormal_params: Cell array of cells, each containing two cells:
%                  - Cell 1: Array of means for the abnormal distribution
%                  - Cell 2: Array of standard deviations for the abnormal distribution
%
% Returns:
% Data   : Matrix of generated PDFs without labels
% x      : X-axis values for plotting the PDFs
% labels : Labels indicating the set to which each PDF belongs
% f      : PDFs with labels

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE

% mu_ranges = {0.25:0.01:0.35, 0.75:0.01:0.8};
% sig_values = [0.2, 0.2];
% abnormal_params = {
%     {[mu1, mu2], [sig1, sig2]}
%     {[mu1, mu2], [sig1, sig2]}
% };
% [Data, x, labels, f] = SimPDFAbnormal(mu_ranges, sig_values, abnormal_params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation

num_groups = length(mu_ranges);
pdfs = [];
labels = [];

for group = 1:num_groups
    mu_range = mu_ranges{group};
    for i = 1:length(mu_range)
        f_single = normpdf(grid, mu_range(i), sig_values(group));
        pdfs = [pdfs; f_single];
        labels = [labels; group];
    end
end

if exist('abnormal_params', 'var') && iscell(abnormal_params)
    for a_group = 1:length(abnormal_params)
        abnormal_mu_sig = abnormal_params{a_group};
        if iscell(abnormal_mu_sig) && numel(abnormal_mu_sig) == 2
            mus = abnormal_mu_sig{1};
            sigmas = abnormal_mu_sig{2};
            abnormal_pdf = zeros(1, length(grid));
            for i = 1:length(mus)
                abnormal_pdf = abnormal_pdf + normpdf(grid, mus(i), sigmas(i));
            end
            pdfs = [pdfs; abnormal_pdf];
            labels = [labels; num_groups + 1]; % Gán nhãn cho nhóm "abnormal"
        end
    end
end

Data = pdfs';
labels = labels';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

% figure;
% hold on;
% colors = hsv(num_groups + length(abnormal_params)); % Tạo bảng màu
% for group = 1:max(labels)
%     group_pdfs = pdfs(labels == group, :);
%     for i = 1:size(group_pdfs, 1)
%         plot(x, group_pdfs(i, :), 'Color', colors(group, :), 'LineWidth', 2);
%     end
% end
% title(sprintf('Plotting of %d distributions', max(labels)));
% xlabel('Value');
% ylabel('Probability');
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving

% name = sprintf("Sim_%d_%d.mat", max(labels), size(pdf, 2));
% save(name, "pdfs", "x", "labels");

end
