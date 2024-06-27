function h = PlotCDEsContour(fvnew, labels)
    % Calculate mean and standard deviation
    v_mean = mean(fvnew,1);
    v_std = std(fvnew,1);
    
    % Determine the number of clusters
    numCluster = max(labels);
    colors = lines(numCluster);
    colororder("meadow")
    
    % Create a figure
    hold on;

    % Plot centroids
    for k = 1:numCluster
    scatter(v_mean(k), v_std(k), 72, colors(k,:), 'filled', "DisplayName", "Centroids", "AlphaData", .7);  % 72 is the marker size
    end 

    % Add labels and title
    xlabel('Mean');
    ylabel('Std');
    title('Scatter Plot with Cluster Coloring');
end
