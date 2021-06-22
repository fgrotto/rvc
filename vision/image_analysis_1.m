clear all;
close all;

% Display the original image
I = imread('eight.tif');
figure;
imshow(I);

% Create a morphological opening
se = strel('disk',8);
background = imopen(I,se);
figure;
imshow(background);

% Increate the contract of the processed image
I2 = imadjust(I);
figure;
imshow(I2);

% Subtract the background approximation image from the 
% original image
% I3 = I2 - background;
% figure;
% imshow(I3);

% Binarize the image obtained at the previous step
bw = imbinarize(I2);
bw = bwareaopen(bw,10);
figure;
imshow(bw);

% Find all connected components from the binarized image
cc = bwconncomp(bw,4);
cc.NumObjects;

% View the rice grain that is labeled 50 in the image
grain = false(size(bw));
grain(cc.PixelIdxList{1}) = true;
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
grain_areas(1);

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
