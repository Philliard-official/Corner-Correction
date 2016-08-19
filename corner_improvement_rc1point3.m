%%% July 7th version
%%% designation: Release Candidate 1.3
%%% 

%clear workspace
clear all
close all
clc


I = imread('web13.jpg'); %read in image

% convert image to grayscale if it isn't already
try
    I = rgb2gray(I);
catch ME
    disp('Image is already grayscale');
end

imshow((I)) %display image to select regoin of interest (ROI)5
hold on % keep image shown even when plotting points

h = imrect; % grab a rectangular ROI from user
[x,y,z] = size(I); % store the dimensions of the image in x, y, and z
t0 = clock; % start time
pos = getPosition(h);
max_distance = 8.7; % maximum radius for clusters
percentPoints = 0.95; % the percentage of points to be selected

corners = detectHarrisFeatures(I, 'ROI', getPosition(h)); 
% corners = corners.selectStrongest(int64(length(corners) * percentPoints)); 
sz = length(corners);
amountOfPoints = sz;
subtotalPoints = [];
finalPoints = [];
dist_count = 0;
one_point = 0;
distances = zeros(sz, (sz-1));

for i = 1:sz(1)
    
    alone = 1;
    cluster = zeros(sz, 3);
    clusterPoints = 1;
    
    
    % establish distances as arrays to save procesessing time
    j = 1:sz(1);
    dx = corners.Location(i,1) - corners.Location(j,1);
    dy = corners.Location(i,2) - corners.Location(j,2);
    
    for p = 1:sz(1)
        
        %calculate distance
        dist = sqrt(dx(p)*dx(p) + dy(p)*dy(p));
        if(dist <= max_distance)
            alone = 0;
            cluster(clusterPoints, 1) = p;
            cluster(clusterPoints, 2) = corners.Location(p, 1);
            cluster(clusterPoints, 3) = corners.Location(p, 2);
            clusterPoints = clusterPoints + 1;
            
        end
    end
    
    if(alone == 1)
        dist_count = dist_count +1;
        subtotalPoints(dist_count, 1) = corners.Location(i,1);
        subtotalPoints(dist_count, 2) = corners.Location(i,2);
        
    else
        total = (sum(cluster))/(clusterPoints-1);
        dist_count = dist_count +1;
        subtotalPoints(dist_count, 1) = int64(total(2));
        subtotalPoints(dist_count, 2) = int64(total(3));
        
    end
    
end

imshow(I);
hold on
disp('The time difference was');
now = clock;% new clock time

disp(now-t0);
disp('Now the points:');

final = cornerPoints(subtotalPoints);
corners = final;% save the new points in corners so final can be used again
plot(corners);

r = 1;
while r<10
    
    r=r+1;
    sz = length(corners);
    subtotalPoints = [];
    finalPoints = [];
    dist_count = 0;
    one_point = 0;
    distances = zeros(sz, (sz-1));

    for i = 1:sz(1)
        alone = 1;
        cluster = zeros(sz, 3);
        clusterPoints = 1;
        j = 1:sz(1);
        dx = corners.Location(i,1) - corners.Location(j,1);
        dy = corners.Location(i,2) - corners.Location(j,2);
        for p = 1:sz(1)
            dist = sqrt(dx(p)*dx(p) + dy(p)*dy(p));
            if(dist <= max_distance)
                alone = 0;
                cluster(clusterPoints, 1) = p;
                cluster(clusterPoints, 2) = corners.Location(p, 1);
                cluster(clusterPoints, 3) = corners.Location(p, 2);
                clusterPoints = clusterPoints + 1;
            end
        end
        
        if(alone == 1)
%             if(sum(find(subtotalPoints(:,1)==corners.Location(i,1)))==0 && sum(find(subtotalPoints(:,2)==corners.Location(i,2)))==0)
                dist_count = dist_count +1;
                subtotalPoints(dist_count, 1) = corners.Location(i,1);
                subtotalPoints(dist_count, 2) = corners.Location(i,2);
%             end
            
        else
            total = (sum(cluster))/(clusterPoints-1);
%             if(isempty(subtotalPoints) || (sum(find(subtotalPoints(:,1)==total(2)))==0 &&sum(find(subtotalPoints(:,2)==total(3)))==0))
                dist_count = dist_count +1;
                subtotalPoints(dist_count, 1) = int64(total(2));
                subtotalPoints(dist_count, 2) = int64(total(3));
%             end
            
        end
        
    end
    
    message = strcat('This iteration was number ', int2str(r));
    disp(message);
    
    disp('The time difference was');
    now = clock;
    disp(now-t0);
    disp('Now the points:');
    
    final = cornerPoints(subtotalPoints);
    corners = final;
    disp(length(final));
    sz = length(corners);
    imshow(I);
    plot(corners);
    
end

disp('The initial amount of points was');
disp(amountOfPoints);

output = zeros(1, 2);
results = 0;
for i=1:length(final)
    if(sum(find(output(:,1) == corners.Location(i, 1))) == 0 && sum(find(output(:,2) == corners.Location(i, 2))) == 0)
        results = results + 1;
        output(results ,1) = corners.Location(i, 1);
        output(results, 2) = corners.Location(i, 2);
    end
    
end

%%% Return output as the list of points.
