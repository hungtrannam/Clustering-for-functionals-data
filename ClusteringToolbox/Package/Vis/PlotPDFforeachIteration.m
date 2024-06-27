function [h]=PlotPDFforeachIteration(data, labels)
% Plot probability density functions from data matrix
% Input:
%   - data: a matrix containing PDF data. Rows represent the x-axis values, and columns represent different PDFs.
%           The last row contains the label for each PDF.
%   - labels: a vector containing the labels for each PDF.
%   - x (optional): a vector of x-axis values. If not provided, a default vector is created.
% Output:
%   - A figure displaying the PDFs with a legend, axis labels, and title.

% Check if x is provided, if not, create a default vector
if nargin < 3
    x = linspace(-.5, 1.5, 500);
end

colors = spring(max(labels));

unique_labels = unique(labels); 
num_labels = numel(unique_labels); 
color_map = containers.Map(unique_labels, mat2cell(colors, ones(size(colors, 1), 1), size(colors, 2)));

h = figure;

for k = 1:4
    subplot(4, 1, k);
    hold on;

for i = 1:numel(labels)
    color = color_map(labels(i));
    plot(x', data{1, i}(:, k), 'Color', color);
end

title_name = {'Red', 'Green', 'Blue', 'Gray'};
title(sprintf('Probability Density Functions of %s channel', title_name{k}));

% Create legend based on color using graphics objects
legend_entries = cell(1, num_labels);
legend_objects = gobjects(1, num_labels);
for j = 1:num_labels
    legend_entries{j} = sprintf('Cluster %d', unique_labels(j));
    legend_objects(j) = plot(NaN, NaN, 'Color', colors(j,:), 'DisplayName', legend_entries{j}); % Plot dummy points for legend
end
legend(legend_objects);
end

end
