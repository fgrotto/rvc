clear all;
close all;

% Display the original image
I = imread('rice.png');
figure;
imshow(I);

% Create a morphological opening
se = strel('disk',15);
background = imopen(I,se);
figure;
imshow(background);

% Subtract the background approximation image from the 
% original image
I2 = I - background;
figure;
imshow(I2);

% Increate the contract of the processed image
I3 = imadjust(I2);
figure;
imshow(I3);

% Binarize the image obtained at the previous step
bw = imbinarize(I3);
bw = bwareaopen(bw,50);
figure;
imshow(bw);

% Find all connected components from the binarized image
cc = bwconncomp(bw,4);
cc.NumObjects;

% View the rice grain that is labeled 50 in the image
grain = false(size(bw));
grain(cc.PixelIdxList{50}) = true;
figure;
imshow(grain);

% Visualize all the connected components in the image
% by creating a label matrix and then diplaying it as a
% pseudocolor indexed image
labeled = labelmatrix(cc);
whos labeled
RGB_label = label2rgb(labeled,'spring','c','shuffle');
figure;
imshow(RGB_label);

% Compute the area of each objec in the image
graindata = regionprops(cc,'basic');
grain_areas = [graindata.Area];
grain_areas(50);

% Find and diplay the grain with the smallest area
[min_area, idx] = min(grain_areas);
grain = false(size(bw));
grain(cc.PixelIdxList{idx}) = true;
figure;
imshow(grain);

% Display an histogram of rice grains areas
figure;
histogram(grain_areas);
title('Histogram of Rice Grain Area');
