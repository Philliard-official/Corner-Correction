%%% July 7th version
%%% designation: Release Candidate 1.3
%%% 

%clear workspace
% clear all
close all
clc


OG = imread('web13.jpg'); %read in image
imshow(OG)
% [x,y,z] = size(OG); % store the dimensions of the image in x, y, and z
% h = imrect(gca, [1 1 y-1 x-1]);
% coord = getPosition(h);
% sub_I = OG(coord(2):coord(2)+coord(4), coord(1):coord(1)+coord(3));
% pause(1)
% filename = strcat('web13tc.jpg');
% disp(filename);
% 
% saveas(gcf, filename);

% convert image to grayscale if it isn't already
try
    I = rgb2gray(OG);
catch ME
    disp('Image is already grayscale');
    I=OG;
end


uniqueI = unique(I(:));
n = histc(I, uniqueI);
[maxNum, maxInd] = max(n);
frequentVals = uniqueI(maxInd);

%check the lowest value
base = abs((double(min(frequentVals))/80));%frequentVals(int64(length(frequentVals)/4))

%value? 0.5528
if(base<0.7 && base > 0.3)
    I= imadjust(I, [base 0.80], []); % adjust contrast to improve results
    
    disp('Adjusted contrast');
end

imshow((I)) %display image to select regoin of interest (ROI)5
hold on % keep image shown even when plotting points
% h = imrect; % grab a rectangular ROI from user
[x,y,z] = size(I); % store the dimensions of the image in x, y, and z
h = imrect(gca, [1 1 y-1 x-1]);
% 
% pause(1)
% filename = strcat('web13tbw.jpg');
% disp(filename);
% imwrite(I, filename);

t0 = clock; % start time
pos = getPosition(h);

% maximum radius for clusters
max_distance = 4.5 + 6.5 * ((x*y)/(2500*2500));

% the percentage of points to be selected
% decrease to increase speed by a few seconds but lower accuracy
percentPoints = 0.85; 
corners = detectHarrisFeatures(I, 'ROI', getPosition(h)); 
% corners = corners.selectStrongest(int64(length(corners) * percentPoints)); 
plot(corners);
original = corners;
% filename = strcat('web13t0.jpg');
% disp(filename);
% imwrite(I, filename);

for r = 1:10
    sz = length(corners);
    amountOfPoints = sz;
    subtotalPoints = zeros(sz(1),2);
    finalPoints = [];
    dist_count = 0;
    one_point = 0;
    % distances = zeros(sz, (sz-1));
    
    X = corners.Location;
    MdlKDT = KDTreeSearcher(X);
    n = size(X, 1);
    idx = randsample(n,5);
    % Y = X;
    Idx = rangesearch(MdlKDT,X,max_distance);
    
    for i = 1:length(Idx)
        indeces = cell2mat(Idx(i));
        cluster = corners.Location(indeces(1,:),:);
        if(size(cluster,1) > 1)
            total = (sum(cluster));
        else
            total = cluster;
        
        end
        solution = total/length(indeces);
        if(i ==1)
            disp('cluster was')
            disp(cluster)
            disp('total was:')
            disp(total)
            disp('solution was:');
            disp(solution);
        end
        if(solution(1)>y ||solution(2)>x)
            disp(strcat('The flaw started at #',int2str(i)))
            disp(solution);
            disp(str(33));
        end
        dist_count = dist_count +1;
        subtotalPoints(dist_count, 1) = int64(solution(1));
        subtotalPoints(dist_count, 2) = int64(solution(2));
    end
%     disp(subtotalPoints)
%     disp(length(subtotalPoints))
%     disp(length(corners))
    
    imshow(I);
    hold on
    disp('The time difference was');
    now = clock;% new clock time
    
    disp(now-t0);
    disp('Now the points:');
    
    final = cornerPoints(subtotalPoints);
    corners = final;% save the new points in corners so final can be used again
    plot(corners);
    plot(corners.Location(:,1),corners.Location(:,2), 'gs', 'markerfacecolor', 'c');
    message = strcat('This iteration was number ', int2str(r));
    disp(message);
%     filename = strcat('web13t', int2str(r), '.jpg');
%     disp(filename);
%     pause(1)
%     saveas(gcf, filename);
    
end

disp('The initial amount of points was');
disp(amountOfPoints);

output = zeros(1, 2);
results = 0;
for i=1:length(final)
%     if(sum(find(output(:,1) == corners.Location(i, 1))) == 0 && sum(find(output(:,2) == corners.Location(i, 2))) == 0)
        results = results + 1;
        output(results ,1) = corners.Location(i, 1);
        output(results, 2) = corners.Location(i, 2);
%     end
    
end

% original = original.selectStrongest(int64(length(corners)/length(original)));
original = original.selectStrongest(int64(length(original)*percentPoints));
for i=1:length(original)
%     if(sum(find(output(:,1) == corners.Location(i, 1))) == 0 && sum(find(output(:,2) == corners.Location(i, 2))) == 0)
        results = results + 1;
        output(results ,1) = original.Location(i, 1);
        output(results, 2) = original.Location(i, 2);
%     end
    
end


% plot(cornerPoints(C));

%%% Return output as the list of points.
xbounds = [min(output(:,1)) min(output(:,2)); max(output(:,1)) max(output(:,2))];

radius = 0.012 + 0.018 * ((x*y)/(2500*2500));
% radius = 0.030
[Z,F] = subclust(output, radius);
% [C, S] = subclust(corners.Location, radius);
imshow(OG);
hold on
plot(Z(:,1),Z(:,2), 'co', 'markerfacecolor' , 'm')


% filename = strcat('web13tf.jpg');
% disp(filename);
% pause(5)
% saveas(gcf, filename);

disp('The time difference was');
now = clock;% new clock time

disp(now-t0);