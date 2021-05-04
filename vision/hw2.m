close all;
clear all;

img = imread('sofa1.png');
figure(1);
imagesc(img);

fx = 525;
cx = 319.5;
fy = 525;
cy = 239.5;
[row,column] = size(img);
cloud = zeros(row*column, 3);

X = zeros(row, column);
Y = zeros(row, column);

k = 0;
threshold = 3000;
for i = 1:1:row
    for j = 1:1:row
        X(j,i) = -(i-cx) * double(img(j,i)) /fx;
        Y(j,i) = -(j-cy) * double(img(j,i)) /fy;
        if (img(j,i)~=0 && img(j,i)<threshold)
            k = k+1;
            cloud(k,:)=[X(j,i) Y(j,i) double(img(j,i))];
        end
    end
    disp(i);
end

cloud = cloud(1:k,:);
figure(2); plot3(cloud(:,1),cloud(:,2),cloud(:,3));

Triangle = [];
exportMeshToPly(cloud, Triangle, ones(size(cloud,1),3), 'mesh');