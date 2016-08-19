total_reps = 1;
files = 13;
rep_num = 1;
load('thisispointer.mat')
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

null = zeros(length(br)*length(sr)*length(bd)*length(sd)*length(pp),1);
null = null(:,1);
positives = null;
false_negatives = null;
false_positives = null;
% cpoints = null;
% zpoints = null;
len_points = null;

%read in image
% or_im = imread(strcat('web',int2str(im),'.jpg'));
% imshow(or_im);
% h0 = imrect;
% pos = getPosition(h0);
% for im = files
im = 13
for percentPoints = pp
    for base_distance = bd
        
        for scaling_distance = sd
            
            for base_radius = br
                
                for scaling_radius = sr
%                                         if (rep_num >= 2000)
                    
                    if(mod(total_reps,100)==0)
                        clc
                    end
                    
                    t0 = clock;
%                     or_im = imread(strcat('web',int2str(files),'.jpg'));
%                     try
%                         I = rgb2gray(or_im);
%                     catch ME
%                         disp('Image is already grayscale');
%                         I = or_im;
%                     end
%                     [x,y,z] = size(or_im);
%                     I = imadjust(I);
%                     h = imrect(gca, [1 1 y-1 x-1]);
%                     pos = getPosition(h);
                    
                    filename = strcat('C:\Users\Public\Documents\Corner-Correction\web ',int2str(im),' points','\web',int2str(im),'_run',int2str(rep_num),'points.mat');
                    load(filename);
                    max_distance = base_distance + scaling_distance * ((pos(3)-pos(1))*(pos(4)-pos(2))/(2500*2500));
                    
                    % now I'll compare the results with the ones I labeled by hand
%                     load(strcat('web',int2str(im),'true_points.mat'));
                    confirmed_positive = 0;
                    confirmed_false_positive = 0;
                    
                    % record the list of points
%                         X = corners.Location;
                        % data structure to compare the points
                        MdlKDT = KDTreeSearcher(true_points);
                        
                        
                        % searches for point indeces within the radius of max_distance
                        Idx = rangesearch(MdlKDT,Z,max_distance);
                        
                    
                    for points = 1:length(Idx)
                        if(~isequal(Idx(points),pointer))
                            confirmed_positive = confirmed_positive+1;
%                             clc
                        else
%                             disp('it didnt go through')
%                             pause(2)
                            confirmed_false_positive = confirmed_false_positive + 1;
                        
                        end
                    end
                    positives(rep_num) = confirmed_positive;
                    false_positives(rep_num) = confirmed_false_positive;
                    len_points(rep_num) = length(Z);
                    %                         cpoints(rep_num) = corners;
                    %                         zpoints(rep_num) = Z;
                    
                        MdlKDT = KDTreeSearcher(Z);
                        Idx = rangesearch(MdlKDT,true_points,max_distance);
                    
                        
                    confirmed_false_negative = 0;
                    for points = 1:length(Idx)
                        % assume false negative first and foremost
                        if(isequal(Idx(points),pointer))
                            confirmed_false_negative = confirmed_false_negative + 1;
                        end
                    end
                    false_negatives(rep_num) = confirmed_false_negative;
                    
                    
%                     disp('The time difference was');
                    now = clock;
                    dif = now-t0;
%                     disp(dif);
                    
                    clear MdlKDT
                        clear Idx
                        clear corners
                        clear Z
                    
                    %                         eta = int64(dif(5)*60+dif(6))
                    eta = int64((dif(5)*60+dif(6))* ((length(br)*length(sr)*length(bd)*length(pp)*length(sd)*length(im))-total_reps)/60);
                    if(mod(total_reps,10)==0)
                        disp('estimated time to completion')
                        disp(strcat(int2str(eta), '~ minutes'));
%                         pause(2);
                    end
                    
%                     end
                    rep_num = rep_num+1;
                    total_reps = total_reps+1;
                end
            end
        end
    end
end
filename = strcat('C:\Users\Public\Documents\Corner-Correction\web',int2str(im),'data_extended2.mat');
save(filename, 'true_points', 'pp', 'br', 'sr','bd','sd','len_points', 'positives', 'false_positives', 'false_negatives');
