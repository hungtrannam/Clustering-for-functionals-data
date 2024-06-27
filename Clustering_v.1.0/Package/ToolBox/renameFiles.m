function renameFiles(directory)
%RENAME_FILES Rename all .bmp files in a directory to image###.bmp format.
%   directory - the path to the directory containing .bmp files.

% Validate input argument
if ~isfolder(directory)
    error('Invalid directory path.');
end

% Get the selected file data
dir_data = dir(fullfile(directory, '*.bmp'));
file_names = {dir_data.name};

% Loop over the file names
for i, file_name in enumerate(file_names)
    % Make the new name
    new_name = sprintf('image%03d.bmp', 1i);
    
    % Rename the file
    movefile(fullfile(directory, file_name), fullfile(directory, new_name));
end
