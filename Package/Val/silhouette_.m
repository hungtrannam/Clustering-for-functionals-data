function cvi = silhouette_(Data, cluster_labels, param)
        uc = unique(cluster_labels);
        K = numel(uc);

        if K == 1 % one cluster only
            cvi = -1;
            return
        end

        n = size(Data,2);
        % precompute the distance matrix
        % Calculate the distance between each pair of PDFs
        for i = 1:n
            for j = 1:n
                if i ~= j
                    % L1 distance between PDFs i and j
                    d(i, j) = trapz(param.x, abs(Data(:, i) - Data(:, j)));
                else
                    % Distance between a PDF and itself should be a small positive number
                    d(i, j) = 10^(-10);
                end
            end
        end


        cvia = zeros(size(cluster_labels)); % array to store a's
        cvib = zeros(size(cluster_labels)); % array to store b's
        for i = 1:K
            this_cluster = cluster_labels == uc(i);
            nk = sum(this_cluster);
            if nk > 1 % guard against one-element cluster
                cvia(this_cluster) = mean(d(this_cluster,this_cluster),2) ...
                    * nk / (nk-1); % correction for self
                bb = [];
                for j = 1:K
                    if i ~= j % other cluster
                        other_cluster = cluster_labels == uc(j);
                        bb = [bb mean(d(this_cluster,other_cluster),2)];
                        % no need for correction for self
                    end
                end
                cvib(this_cluster) = min(bb,[],2);
            else
                cvia(this_cluster) = 1;
                cvib(this_cluster) = 1;
            end
        end
        cvi = mean((cvib - cvia)./max([cvia,cvib],[],2));

    end