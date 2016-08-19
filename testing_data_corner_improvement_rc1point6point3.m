%%% July 29th version
%%% designation: Release Candidate 1.6.3
%%%

%clear workspace
% clear all
% close all
% clc

% start_time = clock;
total_reps = 1;
files = 14;
rep_num = 1;

% base_distance = 4.0;
% scaling_distance = 6.5;
% base_radius = 0.012;
% scaling_radius = 0.018;
% percentPoints = 0.75;


bd = linspace(3, 5, 9);

sd = linspace(5,7,9);

br = linspace((11.0/1000), (13.0/1000), 9);
% br = [11 11.5 12 12.5 13];
% br = double(br)/1000;

sr = linspace((16.0/1000), (18.0/1000), 9);
% sr = [16 17 17.5 18];
% sr = double(sr)/1000;

pp = linspace(0.70,0.9,5);

or_im = imread(strcat('web',int2str(files),'.jpg'));
try
    I = rgb2gray(or_im);
catch ME
    disp('Image is already grayscale');
    I = or_im;
end
[x,y,z] = size(or_im);
I = imadjust(I);
h = imrect(gca, [1 1 y-1 x-1]);
pos = getPosition(h);

% null = zeros(length(br)*length(sr));
% null = null(:,1);
% positives = null;
% false_negatives = null;
% false_positives = null;
% cpoints = null;
% zpoints = null;
% len_points = null;

%read in image
% or_im = imread(strcat('web',int2str(im),'.jpg'));
% imshow(or_im);
% h0 = imrect;
% pos = getPosition(h0);
% for im = files
im = 14
for percentPoints = pp
    for base_distance = bd
        
        for scaling_distance = sd
            
            for base_radius = br
                
                for scaling_radius = sr
                    
                    if (rep_num >= 32356)
%                     disp(pos)
%                     close all
                    if(mod(rep_num,2)==0)
                        clc
                    end
                    disp(strcat('starting run #',int2str(rep_num)))
                    or_im = imread(strcat('web',int2str(im),'.jpg'));
%                     imshow(or_im);
%                     h = h0;
                    % store the dimensions of the image in x, y, and z
%                     [x,y,z] = size(or_im);
                    
                    % convert image to grayscale if it isn't already
%                     try
%                         I = rgb2gray(or_im);
%                     catch ME
%                         disp('Image is already grayscale');
%                         I = or_im;
%                     end
                    
                    % settings
                    
                    scan_width = x;
                    scan_length = y;
                    
                    % adjust contrast to improve results
%                     I = imadjust(I);
                    % DEBUG
%                     disp('Adjusted contrast');
                    
                    %display image
%                     imshow((I));
                    % keep image shown even when plotting points
%                     hold on
                    
                    % grab a rectangular ROI from user
                    
%                     h = imrect(gca, [1 1 y-1 x-1]);
%                     pos = getPosition(h);
                    
                    % start time
                    t0 = clock;
                    
                    % maximum radius for clusters
                    max_distance = base_distance + scaling_distance * ((pos(3)-pos(1))*(pos(4)-pos(2))/(2500*2500));
                    
                    % establish initial corners
                    corners = detectHarrisFeatures(I, 'ROI', pos);
%                     plot(corners);
                    
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
                        
                        
%                         imshow(I);
%                         hold on
                        disp('The time difference was');
                        % new clock time
                        now = clock;
                        
                        disp(now - t0);
                        disp('Now the points:');
                        
                        corners = cornerPoints(subtotalPoints);
%                         plot(corners);
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
%                     imshow(or_im);
%                     hold on
%                     plot(Z(:,1),Z(:,2), 'mo', 'markerfacecolor' , 'm')
                    
                    
                    z_axis = double(zeros(length(Z),1));
                    
                    xscale = x / scan_width;
                    yscale = y / scan_length;
                    
%                     fileID = fopen('exp.txt','w');
%                     text_file = vertcat(vertcat(vertcat(1:length(Z), (xscale*rot90(Z(:,1)))), (yscale*rot90(Z(:,2)))), rot90(z_axis));
%                     fprintf(fileID,'%i\t%6.4f\t%6.4f\t%6.4f\n', text_file);
%                     fclose(fileID);
%                     
                    
                    disp('The time difference was');
                    now = clock;% new clock time
                    
                    disp(now-t0);
                    
                    % now I'll compare the results with the ones I labeled by hand
%                     load(strcat('web',int2str(im),'true_points.mat'));
%                     confirmed_positive = 0;
%                     confirmed_false_positive = 0;
%                     for point1 = 1:length(Z)
%                         pos = 0;
%                         for point2 = 1:length(true_points)
%                             dist = sqrt((Z(point1,1)-true_points(point2,1))*(Z(point1,1)-true_points(point2,1))+(Z(point1,2)-true_points(point2,2))*(Z(point1,2)-true_points(point2,2)));
%                             if(dist < max_distance*1.5)
%                                 confirmed_positive = confirmed_positive + 1;
%                                 pos = 1;
%                                 break;
%                             end
%                         end
%                         if(pos==0)
%                             confirmed_false_positive = confirmed_false_positive + 1;
%                         end
%                     end
%                     positives(rep_num) = confirmed_positive;
%                     false_positives(rep_num) = confirmed_false_positive;
%                     len_points(rep_num) = length(Z);
%                     %                         cpoints(rep_num) = corners;
%                     %                         zpoints(rep_num) = Z;
%                     
%                     confirmed_false_negative = 0;
%                     for point1 = 1:length(true_points)
%                         % assume false negative first and foremost
%                         fn = 1;
%                         for point2 = 1:length(Z)
%                             dist = sqrt((true_points(point1,1)-Z(point2,1))*(true_points(point1,1)-Z(point2,1))+(true_points(point1,2)-Z(point2,2))*(true_points(point1,2)-Z(point2,2)));
%                             
%                             if(dist < max_distance*1.5)
%                                 fn = 0;
%                                 break;
%                             end
%                         end
%                         confirmed_false_negative = confirmed_false_negative + 1;
%                     end
%                     false_negatives(rep_num) = confirmed_false_negative;
%                     
                    
                    disp('The time difference was');
                    now = clock;
                    dif = now-t0;
                    disp(dif);
                    
                    disp('estimated time to completion')
                    %                         eta = int64(dif(5)*60+dif(6))
                    eta = (int64(dif(5)*60+dif(6))* ((length(br)*length(sr)*length(bd)*length(pp)*length(sd)*length(im))-total_reps))/60;
                    disp(strcat(int2str(eta), '~ minutes'));
                    
                    filename = strcat('C:\Users\Public\Documents\Corner-Correction\web',int2str(im),'_run',int2str(rep_num),'points.mat');
                    save(filename, 'corners', 'original', 'Z', 'true_points');
                    disp('The time difference was');
                    now = clock;
                    dif = now-t0;
                    disp(dif);
                    end
                    rep_num = rep_num+1;
                    total_reps = total_reps+1;
%                     pause(5)
                end
            end
        end
    end
end


disp('The time elapsed for all runs was');
now = clock;
dif = now-start_time;
disp(dif);
filename = strcat('C:\Users\Public\Documents\Corner-Correction\web',int2str(im),'data_extended2.mat');
save(filename, 'true_points', 'pp', 'br', 'sr','bd','sd','len_points', 'positives', 'false_positives', 'false_negatives');
run('combination_looper.m')

% end