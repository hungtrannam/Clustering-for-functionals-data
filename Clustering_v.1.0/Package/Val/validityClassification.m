function val = validityClassification(predicted, trueLabels)
    % Calculate classification performance metrics
    accuracy = sum(predicted == trueLabels) / numel(trueLabels);

    % Compute confusion matrix
    C = confusionmat(trueLabels, predicted);
    
    % Calculate sensitivity, specificity, precision, negative predictive value,
    % F1-score, and true positive rate for each label
    numLabels = size(C, 1);
    sensitivity = zeros(numLabels, 1);
    specificity = zeros(numLabels, 1);
    precision = zeros(numLabels, 1);
    npv = zeros(numLabels, 1);
    f1Score = zeros(numLabels, 1);
    
    for i = 1:numLabels
        TP = C(i, i);
        FN = sum(C(i, :)) - TP;
        FP = sum(C(:, i)) - TP;
        TN = sum(C(:)) - TP - FN - FP;
        
        sensitivity(i) = TP / (TP + FN);
        specificity(i) = TN / (TN + FP);
        precision(i) = TP / (TP + FP);
        npv(i) = TN / (TN + FN);
        f1Score(i) = 2 * TP / (2 * TP + FP + FN);
    end

    % Store the classification performance metrics in the output struct
    val.accuracy = accuracy;
    val.sensitivity = sensitivity;
    val.specificity = specificity;
    val.precision = precision;
    val.npv = npv;
    val.f1Score = f1Score;
end
