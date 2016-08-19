%%% July 5th version
%%% designation: Release Candidate 1.1
%%%


I = imread('web12.jpg'); %read in image

% convert image to grayscale if it isn't already
try
    I = rgb2gray(I);
catch ME
    disp('Image is already grayscale');
end

imshow((I)) %display image to select regoin of interest (ROI)5
hold on % keep image shown even when plotting points

h = imrect; % grab a rectangular ROI from user
coord = getPosition(h);
poi = [];
numPOI = 0;
[x,y,z] = size(I); % store the dimensions of the image in x, y, and z
t0 = clock; % start time
max_distance = 4.75; % maximum radius for clusters
x32 = coord(3)/64;
y32 = coord(4)/64;
for d1 = 1:x32
    for d2 = 1:y32
        sub_h = imrect(gca, [(1-mod((d2-1),2)+(d2-1)*64)+coord(1) (1-mod((d1-1),2)+(d1-1)*64)+coord(2) 64 64]);
        pos = getPosition(sub_h);
        corners = detectHarrisFeatures(I, 'ROI',getPosition(sub_h));
        
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
        disp(length(final));
        
        
        corners = final;% save the new points in corners so final can be used again
        plot(corners);
        
        %%
        % sub_I = imcrop(I, getPosition(h));
        % gcf = figure;
        % saveas(gcf, 'im_first_run.jpg');
        %%
        
        r = 1;
        while r<5
            
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
                    if(sum(find(subtotalPoints(:,1)==corners.Location(i,1)))==0 && sum(find(subtotalPoints(:,2)==corners.Location(i,2)))==0)
                        dist_count = dist_count +1;
                        subtotalPoints(dist_count, 1) = corners.Location(i,1);
                        subtotalPoints(dist_count, 2) = corners.Location(i,2);
                    end
                    
                else
                    total = (sum(cluster))/(clusterPoints-1);
                    if(isempty(subtotalPoints) || (sum(find(subtotalPoints(:,1)==total(2)))==0 &&sum(find(subtotalPoints(:,2)==total(3)))==0))
                        dist_count = dist_count +1;
                        subtotalPoints(dist_count, 1) = int64(total(2));
                        subtotalPoints(dist_count, 2) = int64(total(3));
                    end
                    
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
            hold on
            plot(corners);
            %     saveas(gcf, strcat('im_after',int2str(r),'.jpg'));
        end
        
        % sub_I = imcrop(I, getPosition(h));
        
        % saveas(gcf, 'im_end_real.jpg');
        
        disp('The initial amount of points was');
        disp(amountOfPoints);
        
        for iteration = 1:length(subtotalPoints)
            numPOI = numPOI + 1;
            poi(numPOI, 1) = subtotalPoints(iteration,1);
            poi(numPOI, 2) = subtotalPoints(iteration, 2);
        end
        
    end
end


imshow(I);
hold on
poi = cornerPoints(poi);
plot(poi);





% corners = detectHarrisFeatures(I, 'ROI',getPosition(h));
%
%
% % disp(potato);
%
% sz = length(corners);
% amountOfPoints = sz;
% subtotalPoints = [];
% finalPoints = [];
% dist_count = 0;
% one_point = 0;
% distances = zeros(sz, (sz-1));
%
% for i = 1:sz(1)
%
%     alone = 1;
%     cluster = zeros(sz, 3);
%     clusterPoints = 1;
%
%
%     % establish distances as arrays to save procesessing time
%     j = 1:sz(1);
%     dx = corners.Location(i,1) - corners.Location(j,1);
%     dy = corners.Location(i,2) - corners.Location(j,2);
%
%     for p = 1:sz(1)
%
%         %calculate distance
%         dist = sqrt(dx(p)*dx(p) + dy(p)*dy(p));
%         if(dist <= max_distance)
%             alone = 0;
%             cluster(clusterPoints, 1) = p;
%             cluster(clusterPoints, 2) = corners.Location(p, 1);
%             cluster(clusterPoints, 3) = corners.Location(p, 2);
%             clusterPoints = clusterPoints + 1;
%
%         end
%     end
%
%     if(alone == 1)
%         dist_count = dist_count +1;
%         subtotalPoints(dist_count, 1) = corners.Location(i,1);
%         subtotalPoints(dist_count, 2) = corners.Location(i,2);
%
%     else
%         total = (sum(cluster))/(clusterPoints-1);
%         dist_count = dist_count +1;
%         subtotalPoints(dist_count, 1) = int64(total(2));
%         subtotalPoints(dist_count, 2) = int64(total(3));
%
%     end
%
%
%
% end
%
%
%
% imshow(I);
% hold on
% disp('The time difference was');
% now = clock;% new clock time
%
% disp(now-t0);
% disp('Now the points:');
%
% final = cornerPoints(subtotalPoints);
% disp(length(final));
%
%
% corners = final;% save the new points in corners so final can be used again
% plot(corners);
%
% %%
% % sub_I = imcrop(I, getPosition(h));
% % gcf = figure;
% % saveas(gcf, 'im_first_run.jpg');
% %%
%
% r = 1;
% while r<30
%
%     r=r+1;
%     sz = length(corners);
%     subtotalPoints = [];
%     finalPoints = [];
%     dist_count = 0;
%     one_point = 0;
%     distances = zeros(sz, (sz-1));
%
%     for i = 1:sz(1)
%         alone = 1;
%         cluster = zeros(sz, 3);
%         clusterPoints = 1;
%         j = 1:sz(1);
%         dx = corners.Location(i,1) - corners.Location(j,1);
%         dy = corners.Location(i,2) - corners.Location(j,2);
%         for p = 1:sz(1)
%             dist = sqrt(dx(p)*dx(p) + dy(p)*dy(p));
%             if(dist <= max_distance)
%                 alone = 0;
%                 cluster(clusterPoints, 1) = p;
%                 cluster(clusterPoints, 2) = corners.Location(p, 1);
%                 cluster(clusterPoints, 3) = corners.Location(p, 2);
%                 clusterPoints = clusterPoints + 1;
%             end
%         end
%
%         if(alone == 1)
%             if(sum(find(subtotalPoints(:,1)==corners.Location(i,1)))==0 && sum(find(subtotalPoints(:,2)==corners.Location(i,2)))==0)
%                 dist_count = dist_count +1;
%                 subtotalPoints(dist_count, 1) = corners.Location(i,1);
%                 subtotalPoints(dist_count, 2) = corners.Location(i,2);
%             end
%
%         else
%             total = (sum(cluster))/(clusterPoints-1);
%             if(isempty(subtotalPoints) || (sum(find(subtotalPoints(:,1)==total(2)))==0 &&sum(find(subtotalPoints(:,2)==total(3)))==0))
%                 dist_count = dist_count +1;
%                 subtotalPoints(dist_count, 1) = int64(total(2));
%                 subtotalPoints(dist_count, 2) = int64(total(3));
%             end
%
%         end
%
%     end
%
%     message = strcat('This iteration was number ', int2str(r));
%     disp(message);
%
%     disp('The time difference was');
%     now = clock;
%     disp(now-t0);
%     disp('Now the points:');
%
%     final = cornerPoints(subtotalPoints);
%     corners = final;
%     disp(length(final));
%     sz = length(corners);
%     imshow(I);
%     plot(corners);
% %     saveas(gcf, strcat('im_after',int2str(r),'.jpg'));
% end
%
% % sub_I = imcrop(I, getPosition(h));
%
% % saveas(gcf, 'im_end_real.jpg');
%
% disp('The initial amount of points was');
% disp(amountOfPoints);
%
