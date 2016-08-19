%%% July 29th version
%%% designation: Release Candidate 1.6.3
%%% 

%clear workspace
% clear all
close all
clc

%read in image
or_im = imread('web10.jpg'); 

% store the dimensions of the image in x, y, and z
[x,y,z] = size(or_im); 

% convert image to grayscale if it isn't already
try
    I = rgb2gray(or_im);
catch ME
    disp('Image is already grayscale');
    I = or_im;
end

% settings
base_distance = 4.0;
scaling_distance = 6.5;
base_radius = 0.012;
scaling_radius = 0.018;
percentPoints = 0.75; 
scan_width = x;
scan_length = y;

% adjust contrast to improve results
I = imadjust(I); 
% DEBUG
disp('Adjusted contrast');

%display image
imshow((I));
% keep image shown even when plotting points
hold on

% grab a rectangular ROI from user
h = imrect; 
% h = imrect(gca, [1 1 y-1 x-1]);
pos = getPosition(h);

% start time
t0 = clock; 

% maximum radius for clusters
max_distance = base_distance + scaling_distance * ((pos(3)-pos(1))*(pos(4)-pos(2))/(2500*2500));

% establish initial corners
corners = detectHarrisFeatures(I, 'ROI', pos); 
plot(corners);

% save the orginal corners for later
original = corners;


r = 1;
repeat = true;
% repeat for as long as it's making a difference each run
while(repeat)
    last = corners;
    sz = length(corners);
    subtotalPoints = zeros(sz(1),2);
    finalPoints = [];
    dist_count = 0;
    one_point = 0;
    
    % record the list of points
    X = corners.Location;
    % data structure to compare the points
    MdlKDT = KDTreeSearcher(X);
    
    % searches for point indeces within the radius of max_distance
    Idx = rangesearch(MdlKDT,X,max_distance);
    
    
    for i = 1:length(Idx)
        
        % get the indeces of the nearby points
        indeces = cell2mat(Idx(i));
        % record the points at those indeces
        cluster = corners.Location(indeces(1,:),:);
        
        solution = cluster;
        % matlab's mean function will only average correctly 
        % if there's more than one row
        if(size(cluster,1) > 1)
            solution = (mean(cluster));
        end
        
        % record the new point
        dist_count = dist_count + 1;
        subtotalPoints(dist_count, 1) = int64(solution(1));
        subtotalPoints(dist_count, 2) = int64(solution(2));
    end
    
    
    imshow(I);
    hold on
    disp('The time difference was');
    % new clock time
    now = clock;
    
    disp(now - t0);
    disp('Now the points:');
    
    corners = cornerPoints(subtotalPoints);
    plot(corners);
    message = strcat('This iteration was number ', int2str(r));
    disp(message);
    
    % continue only if the new set is different from the last one
    set = setxor(corners.Location, last.Location, 'rows');
    repeat = ~isempty(set);
    r = r+1;
end

disp('The initial amount of points was');
disp(sz);

original = original.selectStrongest(int64(length(original)*percentPoints));
output = vertcat(corners.Location, original.Location);
results = length(output);

%%% Return output as the list of points.
xbounds = [min(output(:,1)) min(output(:,2)); max(output(:,1)) max(output(:,2))];

radius = base_radius + scaling_radius * ((pos(3)-pos(1))*(pos(4)-pos(2))/(2500*2500));
% radius = 0.030
[Z,F] = subclust(output, radius);
% [C, S] = subclust(corners.Location, radius);
imshow(or_im);
hold on
plot(Z(:,1),Z(:,2), 'mo', 'markerfacecolor' , 'm')


z_axis = double(zeros(length(Z),1));

xscale = x / scan_width;
yscale = y / scan_length;

fileID = fopen('exp.txt','w');
text_file = vertcat(vertcat(vertcat(1:length(Z), (xscale*rot90(Z(:,1)))), (yscale*rot90(Z(:,2)))), rot90(z_axis));
fprintf(fileID,'%i\t%6.4f\t%6.4f\t%6.4f\n', text_file);
fclose(fileID);


disp('The time difference was');
now = clock;% new clock time

disp(now-t0);
