function imLab = simpleseg(im, threshVal)
%simpleseg This function performs a simple threshold segmentation procedure
%on a fluorescence image. 

% Median filter the image. 
imNorm = mat2gray(im);
imFilt = medfilt2(imNorm);
% Apply the threshold. 
imThresh = imFilt  > threshVal;

% Get rid of small objects. 
imLarge = bwareaopen(imThresh, 100);

% Label and return the image. 
imLab = bwlabel(imLarge);

end

