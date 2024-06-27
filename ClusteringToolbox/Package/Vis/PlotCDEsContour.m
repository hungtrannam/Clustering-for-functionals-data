function h = PlotCDEsContour(fvnew, labels, param, iter)
    % Calculate mean and standard deviation
    x = param.x;

    for i = 1:size(fvnew,2)
        v_mean(i) = trapz(x', x' .* fvnew(:,i));
    end

    % Determine the number of clusters
    numCluster = max(labels);
    colors = lines(numCluster);
    colororder("meadow")
    
    % Create a figure
    hold on;

    % Plot centroids
    for k = 1:numCluster
    scatter(iter, v_mean(k), 72, colors(k,:), 'filled', "DisplayName", "Centroids");  % 72 is the marker size
    end 

    % Add labels and title
    xlabel('Iteration');
    ylabel('Mean');
    ylim([-5,10]);
    title(sprintf('Mean of prototype PDFs over the %d iterations', iter));
end
