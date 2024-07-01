function [pdf, x] = ExtractKernel(imds, varargin)
% Set default values for optional parameters
defaultX = 1000;
defaultExtensions = {'.png'};

% Parse the input arguments
p = inputParser;
addRequired(p, 'imds', @(x) isa(x, 'matlab.io.datastore.ImageDatastore'));
addParameter(p, 'numPoints', defaultX, @isnumeric);
addParameter(p, 'extensions', defaultExtensions, @iscellstr);
parse(p, imds, varargin{:});

% Get parameter values
imgExtensions = p.Results.extensions;

pdf = [];

% Find all images in the ImageDatastore
imgFiles = imds.Files;

% Check if any images were found
if isempty(imgFiles)
    error('No images found in the ImageDatastore with extension %s', imgExtensions{1});
end

% Process each image in the ImageDatastore
for i = 1:length(imgFiles)
    fullName = imgFiles{i};
    img = imread(fullName);

    if size(img, 3) == 3
        img_gray = rgb2gray(img);
        
    else
        img_gray = img;
    end
    img_gray = double(img_gray);

    
    % img_filtered = img_gray(img_gray > 10 & img_gray < 255);

    h = (4/3)^(1/5) * length(img_gray).^(-1/5) * std(img_gray(:));
    [pdf_img, x] = ksdensity(img_gray(:), 'Bandwidth', h, 'NumPoints', p.Results.numPoints);
    pdf(:, i) = pdf_img';
    textwaitbar(i,length(imgFiles));
end


end
