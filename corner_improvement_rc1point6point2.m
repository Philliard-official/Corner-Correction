%%% July 27th version
%%% designation: Release Candidate 1.6.2
%%% 

%clear workspace
% clear all
close all
clc

%read in image
OG = imread('web13.jpg'); 

% convert image to grayscale if it isn't already
try
    I = rgb2gray(OG);
catch ME
    disp('Image is already grayscale');
    I=OG;
end

base_distance = 4.0;
scaling_distance = 6.5;
base_radius = 0.012;
scaling_radius = 0.018;
percentPoints = 0.75; 

% uniqueI = unique(I(:));
% n = histc(I, uniqueI);
% [maxNum, maxInd] = max(n);
% frequentVals = uniqueI(maxInd);

%check the lowest value
% low_base = abs((double(min(frequentVals))/100));%frequentVals(int64(length(frequentVals)/4))
% high_base = abs((double(max(frequentVals))/100));
%value? 0.5528
% if(base<0.7 && base > 0.3)

% adjust contrast to improve results
I = imadjust(I); 
% DEBUG
disp('Adjusted contrast');

% end
[x,y,z] = size(I); % store the dimensions of the image in x, y, and z
imshow((I)) %display image to select regoin of interest (ROI)5
hold on % keep image shown even when plotting points
% h = imrect; % grab a rectangular ROI from user
h = imrect(gca, [1 1 y-1 x-1]);

t0 = clock; % start time
pos = getPosition(h);

% maximum radius for clusters
max_distance = base_distance + scaling_distance * ((x*y)/(2500*2500));

% the percentage of points to be selected
% decrease to increase speed by a few seconds but lower accuracy

corners = detectHarrisFeatures(I, 'ROI', getPosition(h)); 
% corners2 = detectMinEigenFeatures(I, 'ROI', getPosition(h)); 
% corners3 = detectFASTFeatures(I, 'ROI', getPosition(h)); 
% corners = cornerPoints(vertcat(vertcat(corners1.Location, corners2.Location), corners3.Location));
% corners = corners.selectStrongest(int64(length(corners) * percentPoints)); 
plot(corners);
original = corners;

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
    % Y = X;
    Idx = rangesearch(MdlKDT,X,max_distance);
    
    for i = 1:length(Idx)
        indeces = cell2mat(Idx(i));
        cluster = corners.Location(indeces(1,:),:);
        if(size(cluster,1) > 1)
            solution = (mean(cluster));
        else
            total = cluster;
            solution = total/length(indeces);
        end
        
%         if(i ==1)
%             disp('cluster was')
%             disp(cluster)
%             disp('total was:')
%             disp(total)
%             disp('solution was:');
%             disp(solution);
% %         end
%         if(solution(1)>y ||solution(2)>x)
%             disp(strcat('The flaw started at #',int2str(i)))
%             disp(solution);
%             disp(str(33));
%         end
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
    message = strcat('This iteration was number ', int2str(r));
    disp(message);
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

radius = base_radius + scaling_radius * ((x*y)/(2500*2500));
% radius = 0.030
[Z,F] = subclust(output, radius);
% [C, S] = subclust(corners.Location, radius);
imshow(OG);
hold on
plot(Z(:,1),Z(:,2), 'mo', 'markerfacecolor' , 'm')


z_axis = double(zeros(length(Z),1));
%
scan_width = x;
scan_length = y;
%
xscale = x / scan_width;
yscale = y / scan_length;
% text_file = zeros(length(Z),4);
% text_file = '';
% text_file = zeros(length(Z),1);
% text_file = [];

% fileID = fopen('exp.txt','w');
% text_file = horzcat(horzcat(horzcat(rot90(1:length(Z),3), Z(:,1)), Z(:,2)), z_axis);
% fprintf(fileID,'%i\t%6.4f\t%6.4f\t%6.4\n', text_file);
% fclose(fileID);

fileID = fopen('exp.txt','w');
text_file = vertcat(vertcat(vertcat(1:length(Z), rot90(Z(:,1))), rot90(Z(:,2))), rot90(z_axis));
fprintf(fileID,'%i\t%6.4f\t%6.4f\t%6.4f\n', text_file);
fclose(fileID);

% for pnt_idx = 1:length(Z)
   
%     text_file(pnt_idx) = 'd';
    %strcat(int2str(pnt_idx), '\t', num2str(xscale*double(Z(pnt_idx,1))), '\t',num2str(yscale*double(Z(pnt_idx,2))),'\t',num2str(z_axis));
% end
    
% disp(text_file)

disp('The time difference was');
now = clock;% new clock time

disp(now-t0);