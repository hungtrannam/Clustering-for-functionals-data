function [Data, labels] = SimPDFAbnormal_Uniform(mu_ranges, sig_values, grid, abnormal_params)
% SimPDFAbnormal_Uniform generates and plots PDFs for multiple sets of distributions
% with specified ranges for the uniform distribution, including optional
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE

% mu_ranges = {0.25:0.01:0.35, 0.75:0.01:0.8};
% sig_values = [0.2, 0.2];
% abnormal_params = {
%     {[0.4, 0.5], [0.1, 0.2]}
%     {[0.6, 0.7], [0.1, 0.2]}
% };
% [Data, x, labels] = SimPDFAbnormal_Uniform(mu_ranges, sig_values, abnormal_params);

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
        % Define the bounds for the uniform distribution
        a = mu_range(i) - sig_values(group);
        b = mu_range(i) + sig_values(group);
        f_single = unifpdf(grid, a, b);
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
                % Define the bounds for the uniform distribution
                a = mus(i) - sigmas(i);
                b = mus(i) + sigmas(i);
                abnormal_pdf = abnormal_pdf + unifpdf(grid, a, b);
            end
            pdfs = [pdfs; abnormal_pdf];
            labels = [labels; num_groups + 1]; % Label for "abnormal" group
        end
    end
end

Data = pdfs';
labels = labels';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting (uncomment if needed for visualization)

% figure;
% hold on;
% colors = hsv(num_groups + length(abnormal_params)); % Generate color map
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
% Saving (uncomment if needed for saving the data)

% name = sprintf("Sim_%d_%d.mat", max(labels), size(pdfs, 2));
% save(name, "pdfs", "x", "labels");

end
