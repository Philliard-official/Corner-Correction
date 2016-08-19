%%%clustertest
clear all
close all
clc

I = imread('web8.jpg'); %read in image

% convert image to grayscale if it isn't already
try
    I = rgb2gray(I);
catch ME
    disp('Image is already grayscale');
end

imshow((I)) %display image to select regoin of interest (ROI)5
hold on % keep image shown even when plotting points

% h = imrect; % grab a rectangular ROI from user
[x,y,z] = size(I); % store the dimensions of the image in x, y, and z
t0 = clock; % start time
% pos = getPosition(h);
max_distance = 5.7; % maximum radius for clusters
corners = detectHarrisFeatures(I);%, 'ROI',getPosition(h)); 
sz = length(corners);
amountOfPoints = sz;
subtotalPoints = [];
finalPoints = [];
dist_count = 0;
one_point = 0;
distances = zeros(sz, (sz-1));
% k = int64(length(corners)/3);

% [idx C] = kmeans(corners.Location, k);
X = corners.Location;
xbounds = [min(X(:,1)) min(X(:,2)); max(X(:,1)) max(X(:,2))];

[C,S] = subclust(X,0.007, xbounds);

imshow(I);
hold on
% plot(corners);
plot(C(:,1),C(:,2), 'go', 'markerfacecolor' , 'g');
% plot(cornerPoints(C));
% 
% newPOI = cornerPoints(C);
% plot(newPOI);

