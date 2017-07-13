function segMask = LogSegmentation(im, A, E)
%LOGSEGMENTATION Perform Laplacian of Gaussian segmentation algorithm over a
%field of bacterial cells. 
% 	imSeg = LOGSEGMENTATION(im) segments im removing objects smaller than
% 		100 pixels
%	imSeg = LOGSEGMENTATION(im, A) segments im removing objects smaller
%		than N pixels.

%Ensure that the image is normalized.
if max(max(im)) > 1
	im = mat2gray(im);
end

%Instantiate the edgefinder function.
imEdge = edge(im, 'log', 0.001);

%Fill the holes in the edge image.
imFill = imfill(imEdge, 'holes');

%Remove the small stuff.
imFilt = bwareaopen(imFill, 100, 4);

%Dilate the objects by 2px.
selem = strel('diamond', 1);
imDil = imdilate(imFilt, selem);

%Remove objects touching the border and label the output image.
imSeg = imclearborder(imDil, 4); %Buffer of 4px.
imSeg = bwlabel(imSeg);

%Extract the region props of the mask.
props = regionprops(imSeg, 'Area', 'Eccentricity');

%Iterate through each object and apply an area and eccentricity filter.
areas = [props.Area];
eccs = [props.Eccentricity];
approvedAreas = (areas > A(1)) & (areas < A(2));
approvedEccs = (eccs > E(1)) & (eccs < E(2));

%Find the approved labels
approvedLabels = find((approvedAreas .* approvedEccs) > 0);

%Generate the final mask.
segMask = zeros(size(imSeg));
for i=1:length(approvedLabels)
	segMask = segMask + (imSeg==approvedLabels(i));
end

%Relabel the image.
segMask = bwlabel(segMask > 0);

end



