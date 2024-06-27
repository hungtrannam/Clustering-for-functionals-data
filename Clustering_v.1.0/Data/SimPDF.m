function [Data, x, labels, f] = SimPDF(mu_ranges, sig_values, range)
% SimPDFAnormal generates and plots PDFs for multiple sets of distributions
% with specified means and standard deviations.
%
% Parameters:
% mu_ranges : Cell array of ranges of means for each set of distributions
% sig_values: Array of standard deviations for all distributions
%
% Returns:
% pdfs   : Matrix of generated PDFs without labels
% x      : X-axis values for plotting the PDFs
% labels : Labels indicating the set to which each PDF belongs
% f      : PDFs with labels

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE

% mu_ranges = {0.25:0.01:0.35, 0.75:0.01:0.8, 0.55:0.01:0.6};
% sig_values = [0.2, 0.2, 0.2];
% 
% [Data, x, labels] = SimPDFAnormal(mu_ranges, sig_values);


close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting
a = range{1};
b = range{2};

x = linspace(a, b, 1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation

num_groups = length(mu_ranges);
pdfs = [];

for group = 1:num_groups
    mu_range = mu_ranges{group};
    f = [];
    for i = 1:length(mu_range)
        f_single = normpdf(x, mu_range(i), sig_values(group));
        f = [f; f_single group];
    end
    pdfs = [pdfs; f];
end

Data = pdfs(:,1:end-1)';
labels = pdfs(:,end)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

% figure;
% hold on;
%
% % Plot the first set of PDFs
% plot(x', f1(:,1:end-1)', ...
%     'Color', 'g', ...
%     'LineWidth', 2);
%
% % Plot the second set of PDFs
% plot(x', f2(:,1:end-1)', ...
%     'Color', 'b', ...
%     'LineWidth', 2);
%
% % Plot the second set of PDFs
% plot(x', f3(:,1:end-1)', ...
%     'Color', 'k', ...
%     'LineWidth', 2);
%
% hold off;
% title(sprintf('Plotting of %d distributions', size(pdf, 2)));
% xlabel('Value');
% ylabel('Probability');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving

% name = sprintf("Sim_%d_%d.mat", max(labels), size(pdf, 2);
% save(name, ...
%     "pdfs", "x", "labels");

end
