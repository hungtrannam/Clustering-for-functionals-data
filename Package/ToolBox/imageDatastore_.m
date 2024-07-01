function [imagePaths, labels] = imageDatastore_(Path_data, fileExtension, includeSubfolders)
    % Ensure the image package is loaded
    pkg load image;

    % Initialize cell arrays to store image paths and labels
    imagePaths = {};
    labels = {};

    % If includeSubfolders is true, get a list of all subfolders
    if includeSubfolders
        subfolders = dir(Path_data);
        subfolders = subfolders([subfolders.isdir] & ~strncmpi('.', {subfolders.name}, 1));
    else
        subfolders = dir(Path_data);
        subfolders = subfolders([subfolders.isdir]);
        subfolders = subfolders(1); % Only the main folder
    end

    % Loop through each folder (subfolder or main folder)
    for k = 1:length(subfolders)
        if includeSubfolders
            subfolderPath = fullfile(Path_data, subfolders(k).name);
        else
            subfolderPath = Path_data;
        end

        % Get all files with the specified extension in the current folder
        imageFiles = dir(fullfile(subfolderPath, ['*' fileExtension]));

        % Loop through each image file
        for j = 1:length(imageFiles)
            imagePath = fullfile(subfolderPath, imageFiles[j].name);

            % Store the image path and corresponding label (folder name)
            imagePaths{end+1} = imagePath;
            if includeSubfolders
                labels{end+1} = subfolders(k).name;
            else
                labels{end+1} = ''; % No label if not using subfolders
            end
        end

        % If not including subfolders, break after the first iteration
        if ~includeSubfolders
            break;
        end
    end
end

