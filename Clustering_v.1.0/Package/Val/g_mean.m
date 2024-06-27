function Gmean = g_mean(labels1, labels2)
    % The function calculates the Rand index (RI) and G-mean
    % between two label assignments: labels1 and labels2.
    
    N = numel(labels1);  % Get the number of elements in the label vector
    
    % Initialize the four quantities: TP (true positive), FN (false negative), FP (false positive), TN (true negative)
    TP = 0; FN = 0; FP = 0; TN = 0;
    
    % Calculate TP, FN, FP and TN
    for i = 1:N-1
        for j = i+1:N
            if (labels1(i) == labels1(j)) && (labels2(i) == labels2(j))  % TP: Both labels1 and labels2 have the same class for samples i and j
                TP = TP + 1;
            elseif (labels1(i) == labels1(j)) && (labels2(i) ~= labels2(j))  % FN: labels1 have the same class for samples i and j, but labels2 not
                FN = FN + 1;
            elseif (labels1(i) ~= labels1(j)) && (labels2(i) == labels2(j))  % FP: labels2 have the same class for samples i and j, but labels1 not
                FP = FP + 1;
            else   % TN: Both labels1 and labels2 have different classes for samples i and j
                TN = TN + 1;
            end
        end
    end
    
    % Calculate Sensitivity (Recall) and Specificity
    Sensitivity = TP / (TP + FN);  % True Positive Rate
    Specificity = TN / (TN + FP);  % True Negative Rate
    
    % Calculate G-mean
    Gmean = sqrt(Sensitivity * Specificity);
end
