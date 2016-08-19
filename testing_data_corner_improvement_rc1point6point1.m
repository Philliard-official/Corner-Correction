%%% July 7th version
%%% designation: Release Candidate 1.6.1
%%%

%clear workspace
clear all
close all
clc
start_time = clock;
total_reps = 1;

% im = 13;
% for im = [13 14]
% for im = [9 14]
im = 13
    rep_num = 1;
    disp(im)
%     dk = 70:85;
%     dk = dk(mod((dk),5) == 0);
    dk = linspace(70,90,3);
    % dark_value = 75;
    
%     bd = 30:40;
%     bd = double(bd(mod((bd),5) == 0))/10;
%     bd = linspace(3,4,5);
% bd = [3.5 4]
bd = 4.0;
    % base_distance = 30;
    
%     sd = 60:70;
%     sd = double(sd(mod((sd),5) == 0))/10;
    sd = linspace(6,7,5);
    % scaling_distance = 60;
    
    br = [11 11.5 12 12.5 13];
    br = double(br)/1000;
%     br = linspace(11,13,5)
    % base_radius = 0.011;
    
%     sr = [16 16.5 17 17.5 18];
    sr = [16 17 17.5 18];
    sr = double(sr)/1000;
    % scaling_radius = 0.016;
    
%     pp = [0.9];
    pp = linspace(0.85,0.95,3);
%     pp = 85:95;
%     pp = double(pp(mod((pp),5) == 0))/100;
    
%     pp = 84:96;
%     pp = double(pp(mod((pp), 3) == 0))/100;
%     
    
    null = zeros(length(br)*length(sr));
    null = null(:,1);
    positives = null;
    false_negatives = null;
    false_positives = null;
    cpoints = null;
    zpoints = null;
    % len_points = null;
    % for dark_value = dk
    %     for base_distance = bd
    %       for scaling_distance = sd
    %           for base_radius = br
    %             for scaling_radius = sr
    %                 for percent_points = pp
    %                     disp('why is this so slow')
    %                 end
    %             end
    %           end
    %       end
    %     end
    % end
    % disp(str(4));
    % dark_value = dk(1);
    % base_distance = bd(1);
    % scaling_distance = sd(1);
    % base_radius = br(1);
    % scaling_radius = sr(1);
    % for im=images
    
    for dark_value = dk
        
        for percentPoints = pp
            for base_distance = bd
                
                for scaling_distance = sd
                    
                    for base_radius = br
                        
                        for scaling_radius = sr
                            
                            close all
                            clc
                            
                            disp(strcat('starting run #',int2str(rep_num)))
                            OG = imread(strcat('web',int2str(im),'.jpg')); %read in image
                            
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
                            
                            %                 dark_value = 80;
                            %check the lowest value
                            base = abs((double(min(frequentVals))/dark_value));%frequentVals(int64(length(frequentVals)/4))
                            
                            %value? 0.5528
                            upper_contrast = 0.7;
                            lower_contrast = 0.3;
                            if(base<upper_contrast && base > lower_contrast)
                                I= imadjust(I, [base 0.80], []); % adjust contrast to improve results
                                
                                disp('Adjusted contrast');
                            end
                            
                            %         imshow((I)) %display image to select regoin of interest (ROI)5
                            %         hold on % keep image shown even when plotting points
                            
                            %         h = imrect; % grab a rectangular ROI from user
                            
                            [x,y,z] = size(I); % store the dimensions of the image in x, y, and z
                            t0 = clock; % start time
                            h = imrect(gca, [1 1 y-1 x-1]);
                            pos = getPosition(h);
                            
                            % maximum radius for clusters
                            %                 base_distance = 4.5;
                            %                 scaling_distance = 6.5;
                            max_distance = base_distance + scaling_distance * ((x*y)/(2500*2500));
                            
                            % the percentage of points to be selected
                            % decrease to increase speed by a few seconds but lower accuracy
                            %                 percentPoints = 0.85;
                            corners = detectHarrisFeatures(I, 'ROI', getPosition(h));
                            imshow(I)
                            hold on
                            plot(corners)
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
                                idx = randsample(n,5);
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
                                    
%                                             if(i ==1)
%                                                 disp('cluster was')
%                                                 disp(cluster)
%                                                 disp('total was:')
%                                                 disp(total)
%                                                 disp('solution was:');
%                                                 disp(solution);
%                                             end
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
                                message = strcat('This iteration was number ', int2str(r));
                                disp(message);
                            end
                            
                            % disp('The initial amount of points was');
                            % disp(amountOfPoints);
                            
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
                            
                            %         base_radius = 0.012;
                            %         scaling_radius = 0.018;
                            radius = base_radius + scaling_radius * ((x*y)/(2500*2500));
                            % radius = 0.030
                            [Z,F] = subclust(output, radius);
                            % [C, S] = subclust(corners.Location, radius);
                            %         imshow(OG);
                            %         hold on
                            
                            %         plot(Z(:,1),Z(:,2), 'mo', 'markerfacecolor' , 'm');
                            
                            
                            % now I'll compare the results with the ones I labeled by hand
                            load(strcat('web',int2str(im),'true_points.mat'));
                            %         plot(true_points(:,1), true_points(:,2), 'bo', 'markerfacecolor', 'b');
                            confirmed_positive = 0;
                            confirmed_false_positive = 0;
                            for point1 = 1:length(Z)
                                pos = 0;
                                for point2 = 1:length(true_points)
                                    dist = sqrt((Z(point1,1)-true_points(point2,1))*(Z(point1,1)-true_points(point2,1))+(Z(point1,2)-true_points(point2,2))*(Z(point1,2)-true_points(point2,2)));
                                    if(dist < max_distance*1.5)
                                        confirmed_positive = confirmed_positive + 1;
                                        pos = 1;
                                        break;
                                    end
                                end
                                if(pos==0)
                                    confirmed_false_positive = confirmed_false_positive + 1;
                                end
                            end
                            positives(rep_num) = confirmed_positive;
                            false_positives(rep_num) = confirmed_false_positive;
                            len_points(rep_num) = length(Z);
                            %                         cpoints(rep_num) = corners;
                            %                         zpoints(rep_num) = Z;
                            
                            confirmed_false_negative = 0;
                            for point1 = 1:length(true_points)
                                % assume false negative first and foremost
                                fn = 1;
                                for point2 = 1:length(Z)
                                    dist = sqrt((true_points(point1,1)-Z(point2,1))*(true_points(point1,1)-Z(point2,1))+(true_points(point1,2)-Z(point2,2))*(true_points(point1,2)-Z(point2,2)));
                                    
                                    if(dist < max_distance*1.5)
                                        fn = 0;
                                        break;
                                    end
                                end
                                confirmed_false_negative = confirmed_false_negative + 1;
                            end
                            false_negatives(rep_num) = confirmed_false_negative;
                            
                            
                            disp('The time difference was');
                            now = clock;
                            dif = now-t0;
                            disp(dif);
                            
                            disp('estimated time to completion')
                            %                         eta = int64(dif(5)*60+dif(6))
                            eta = (int64(dif(5)*60+dif(6))* ((length(dk)*length(br)*length(sr)*length(bd)*length(pp)*length(sd)*length(im))-total_reps))/60;
                            disp(strcat(int2str(eta), '~ minutes'));
                            
                            filename = strcat('web',int2str(im),'_run',int2str(rep_num),'points.mat');
                            save(filename, 'corners', 'original', 'Z', 'true_points');
                            rep_num = rep_num+1;
                            total_reps = total_reps+1;
                             pause(1)
                        end
                    end
                end
            end
        end
    end
    disp('The time elapsed for all runs was');
    now = clock;
    dif = now-t0;
    disp(dif);
    filename = strcat('web',int2str(im),'data_extended2.mat');
    save(filename, 'true_points', 'dk','pp', 'br', 'sr','bd','sd','len_points', 'positives', 'false_positives', 'false_negatives');
    run('combination_looper.m')
    
% end

