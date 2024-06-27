function results = H_(Data, param)
    
    threshold = .5;

    n = size(Data, 2); % Number of PDFs
    
    % Initialize the distance matrix
    D = zeros(n, n);
    
    % Calculate the distance between each pair of PDFs
    for i = 1:n
        for j = 1:n
            if i ~= j
                % L1 distance between PDFs i and j
                D(i, j) = trapz(param.x, abs(Data(:, i) - Data(:, j)));
            else
                % Distance between a PDF and itself should be a small positive number
                D(i, j) = 10^(-10);
            end
        end
    end
    
    % Perform hierarchical clustering using the distance matrix
    % 'average' linkage (UPGMA - Unweighted Pair Group Method with Arithmetic mean)
    Z = linkage(D, 'average');

    IDX = cluster(Z, 'cutoff', threshold, 'criterion', 'distance');



    results.Data.Data = Data;
    results.Cluster.IDX = IDX;
    results.Dist.D = D; 

end