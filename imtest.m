%clear workspace
clear all
close all
clc

I = imread('web8.jpg'); %read in image

% convert image to grayscale
try
    I = rgb2gray(I);
catch ME
    disp('Image is already grayscale');
end


% I = im2bw(I, 0.6); % convert to binary black and white
% I = imcomplement(I); % this inverts the image

imshow((I)) %display image to select regoin of interest (ROI)5
hold on % keep image shown even when plotting points

% This section is used to get a sample image
h = imrect; %grab a rectangular ROI from user


%-------------------------------------------------
% make a subimage of I
% coord = getPosition(h);
% disp(h);
% sub_I = imcrop(I, coord);
% % sub_I = I(coord(2):coord(2)+coord(4), coord(1):coord(1)+coord(3));
% 
% imshow(sub_I);
%  
% h = imrect; %grab a rectangular ROI from user
% coord = getPosition(h);

% % sub_I2 = sub_I(coord(2):coord(2)+coord(4), coord(1):coord(1)+coord(3));
% sub_I2 = imcrop(sub_I, coord);
% 
% 
% % imwrite(sub_I2, 'webtestfile2.jpg');
% % imwrite(sub_I2, 'intersection16.jpg');
% imshow(sub_I2);
% % imwrite(sub_I2, 'sample1.jpg');
% imwrite(sub_I, 'neg23.jpg');

% disp(getPosition(h)); %DEBUG

[x,y,z] = size(I);

% --------------------------------------------
% This segment

% sqSize = 128;
% xsect = int64(x/sqSize);
% ysect = int64(y/sqSize);
% for i = 1:xsect
%     for j = 1:ysect
%         
% %         h = imrect(gca, [(1-(j-1)+(j-1)*sqSize) (1-(i-1)+(i-1)*sqSize) sqSize sqSize]);
%         h = imrect(gca, [(1-mod((j-1),2)+(j-1)*sqSize) (1-mod((i-1),2)+(i-1)*sqSize) sqSize sqSize]);
%         coord = getPosition(h);
%         sub_I = I(coord(2):coord(2)+coord(4), coord(1):coord(1)+coord(3));
%         classification = int2str((1000*i)+(j));
%         filename = strcat('web9_', classification, '.jpg');
%         disp(filename);
%         imwrite(sub_I, filename);
% %         [filename, 0] = imsave(h);
%         
%     end
% end
%-----------------------------------------------


%--------------------------------------------------
% Select corners in smaller chunks
% disp(x);
% disp(y);
% for i = 1:5
%     for j = 1:5
% %         h = imrect(gca, [(1-(j-1)+(j-1)*int64(y/5)) (1-(i-1)+(i-1)*int64(x/5)) int64(y/5) int64(x/5)]);
%         pos = getPosition(h);
%         disp(pos);
%         hs = imrect(gca, [(pos(1)-(j-1)+(j-1)*int64(pos(3)/5)) (pos(2)-(i-1)+(i-1)*int64(pos(4)/5)) int64(pos(3)/5) int64(pos(4)/5)]);
%         
%         corners = detectHarrisFeatures(I,'ROI',getPosition(hs));
% %         corners = detectMinEigenFeatures(I,'ROI',int16(getPosition(h)),'minquality',.2);
%         hold on
%         blackcount = sum(I(1,:));
%         whitecount = sum(I(:,1));
%         
%         ratio = (pos(3)*pos(4))/(x*y);
%         
% %         plot(corners.selectStrongest(5));
%         plot(corners.selectStrongest(int64(ratio*80)));
% %         plot(corners.selectStrongest(int64(175*blackcount/whitecount)));
% %         disp(int64(175*blackcount/whitecount));
%         if i==3 && j==3
%             disp(corners.selectStrongest(5));
%         end
%        delete(hs);
% %         hs(i, j, 1) = 0;
% %         hs(i, j, 2) = 0;
%     end
% end
%---------------------------------------------------


% disp(hs);
% disp(datestr(now,'dd-mm-yyyy HH:MM:SS FFF'))
t0 = clock; % start time
pos = getPosition(h);
disp(pos);
amountOfPoints = int64(5000*(pos(3)*pos(4))/(x*y));
disp(amountOfPoints);

%set the threshold for two points being "near" each other
max_distance = 8.9;

%get corners within ROI, different methods are available
corners = detectHarrisFeatures(I,'ROI',getPosition(h)); 
% corners = detectFASTFeatures(I, 'ROI',getPosition(h));
% corners = detectMinEigenFeatures(I,'ROI',int16(getPosition(h)),'minquality', 0.2);

% corners = corners.selectStrongest(amountOfPoints);
amountOfPoints = length(corners);

plot(corners);
% stop = imrect; % DEBUG to make the code pause

sz = length(corners);
disp('The size is:'); % DEBUG
disp(sz); % DEBUG

test = zeros(sz, sz, 2);

% disp('stop location: '); % DEBUG
% disp(getPosition(stop));
subtotalPoints = []; % zeros(sz,2);
finalPoints = [];
dist_count = 0;
one_point = 0;
distances = zeros(sz, (sz-1));

% disp(corners.Location(1,:)); % DEBUG

% check every point against every other point
for i = 1:sz(1)
    disp('point in question:');
    disp(strcat(int2str(corners.Location(i,1)), ',', int2str(corners.Location(i,2))));
    alone = 1;
    cluster = zeros(sz, 3);
    clusterPoints = 1;
    j = 1:sz(1);
    dx = corners.Location(i,1) - corners.Location(j,1);% 
%     disp(dx); % DEBUG
    
    dy = corners.Location(i,2) - corners.Location(j,2);
    for p = 1:sz(1)
        dist = sqrt(dx(p)*dx(p) + dy(p)*dy(p));
%         disp(dist);
        if(dist <= max_distance)
            alone = 0;
            cluster(clusterPoints, 1) = p;
            cluster(clusterPoints, 2) = corners.Location(p, 1);
            cluster(clusterPoints, 3) = corners.Location(p, 2);
            clusterPoints = clusterPoints + 1;
%             disp(strcat('distance was good: ',int2str(p))); % DEBUG
            
        end
    end
    
%     dist = sqrt(dx*dx + dy*dy);
%     goodDist = distances < max_distance;
%     disp('cluster');
%     disp(cluster);
    
%     dx = corners.Location(i,1) - corners.Location(j,1);
%     dy = corners.Location(i,2) - corners.Location(j,2);
    
%     for j = 1:sz(1)
%         
%         % don't check a point against itself
%         if(i~=j)
%             dx = corners.Location(i,1) - corners.Location(j,1);
%             dy = corners.Location(i,2) - corners.Location(j,2);
%             dist = sqrt(dx*dx + dy*dy);
%             distances(i,j) = dist;
%             
%             if(max_distance >= dist)
% %                  disp(dist);
%                  alone = 0;
%                  dist_count = dist_count + 1;
%                  subtotalPoints(dist_count, 1) = int64((corners.Location(i,1) + corners.Location(j,1))/2);
%                  subtotalPoints(dist_count, 2) = int64((corners.Location(i,2) + corners.Location(j,2))/2);
%                  one_point = one_point+1;
% %                     disp(strcat(int2str(corners.Location(i,1)), ',',
% %                     int2str(corners.Location(i,2)), '_and_',
% %                     int2str(corners.Location(j,1)), ',',
% %                     int2str(corners.Location(j,2)), '__average to__',
% %                     int2str(subtotalPoints(dist_count, 1)), ', ',
% %                     int2str(subtotalPoints(dist_count, 2)))); %DEBUG
% % %                 test(i,j,1) = (corners.Location(i,1) + corners.Location(j,1))/2;
% % %                 test(i,j,2) = (corners.Location(i,2) + corners.Location(j,2))/2;
% 
%             end
%         
%         
%         end
%     end
    
    
    if(alone == 1)
        dist_count = dist_count +1;
        subtotalPoints(dist_count, 1) = corners.Location(i,1);
        subtotalPoints(dist_count, 2) = corners.Location(i,2);
        disp('created point:');
        disp(strcat(int2str(subtotalPoints(dist_count, 1)), ',', int2str(subtotalPoints(dist_count, 2))));
%         disp(strcat('plotting__', int2str(subtotalPoints(dist_count, 1)), ', ', int2str(subtotalPoints(dist_count, 2))));
    else
        total = (sum(cluster))/(clusterPoints-1);
        disp(total);
        dist_count = dist_count +1;
        subtotalPoints(dist_count, 1) = int64(total(2));
        subtotalPoints(dist_count, 2) = int64(total(3));
        disp('created point:');
        disp(strcat(int2str(subtotalPoints(dist_count, 1)), ',', int2str(subtotalPoints(dist_count, 2))));
    end
%     disp(strcat('Printing all similar points ',num2str(i)));
%     disp(subtotalPoints);
end



% disp('Take the picture.');
% pause(15);


% z=imrect;

% N=hist(distances,length(distances));
                    

% datestr(now,'dd-mm-yyyy HH:MM:SS FFF')
imshow(I);
hold on
disp('The time difference was');
now = clock;% new clock time
% disp(now);
% disp(now(6));
disp(now-t0);
disp('Now the points:');

final = cornerPoints(subtotalPoints);


% disp(corners);
% disp(subtotalPoints);
% disp(final);
% disp(final.Location);
disp(length(final));
% disp(final.Location(1,2));
corners = final;



plot(corners);



% plot(subtotalPoints);
% plot(corners);

% DEBUG
% disp(test);
% disp(test(21));
% test = [100 100
%         1000 100
%         100 1000];



% test = (cornerPoints(test));
% disp(test);
% plot(test.selectStrongest(2));

% hold on %"hold" image to plot detected corner points
% plot(corners.selectStrongest(amountOfPoints), 'Color', 'blue'); %plot corners
% h=imrect;

% disp(corners);
% plot(corners);

% plot(corners.selectStrongest(1000));
% d = corners.selectStrongest(5);
% d1 = d(1);
% 
% disp(d1);


%-----------------------------------------------------
%extra loop (old code is below)
% 
x = 1;
while x<20
    
    x=x+1;
%     if (x == 100)
%         x = input('repeat?');
%         if(x~=0)
%             x=x-1;
%         end
%     end
    
    
    sz = length(corners);
%     disp(sz);
    test = zeros(sz, sz, 2);
%     disp('stop location: ');
    % disp(getPosition(stop));
    subtotalPoints = [];% zeros(sz,2);
    finalPoints = [];
    dist_count = 0;
    one_point = 0;
    distances = zeros(sz, (sz-1));

    for i = 1:sz(1)
%         disp('point in question:');
%         disp(strcat(int2str(corners.Location(i,1)), ',', int2str(corners.Location(i,2))));
        alone = 1;
        cluster = zeros(sz, 3);
        clusterPoints = 1;
        j = 1:sz(1);
        dx = corners.Location(i,1) - corners.Location(j,1);% 
    %     disp(dx); % DEBUG
%         if(i==3)
%             disp(strcat('point 3: (',int2str(corners.Location(i,1)),',',int2str(corners.Location(i,2)),')'));
%         end
        
        dy = corners.Location(i,2) - corners.Location(j,2);
        for p = 1:sz(1)
            dist = sqrt(dx(p)*dx(p) + dy(p)*dy(p));
            if(dist <= max_distance)
%                 disp(dist);
                alone = 0;
                cluster(clusterPoints, 1) = p;
                cluster(clusterPoints, 2) = corners.Location(p, 1);
                cluster(clusterPoints, 3) = corners.Location(p, 2);
                clusterPoints = clusterPoints + 1;
    %             disp(strcat('distance was good: ',int2str(p))); % DEBUG

            end
        end
        
        %This part was debugging by checking the distances for a single
        %point
%         if(i == 3)
%             disp(cluster);
%         end

        if(alone == 1)
            %if it's not within the array
            if(sum(find(subtotalPoints(:,1)==corners.Location(i,1)))==0 && sum(find(subtotalPoints(:,2)==corners.Location(i,2)))==0)
                dist_count = dist_count +1;
                subtotalPoints(dist_count, 1) = corners.Location(i,1);
                subtotalPoints(dist_count, 2) = corners.Location(i,2);
%                 disp('created point:');
%                 disp(strcat(int2str(subtotalPoints(dist_count, 1)), ',', int2str(subtotalPoints(dist_count, 2))));

            end
%             dist_count = dist_count +1;
%             subtotalPoints(dist_count, 1) = corners.Location(i,1);
%             subtotalPoints(dist_count, 2) = corners.Location(i,2);
%             disp('created point:');
%             disp(strcat(int2str(subtotalPoints(dist_count, 1)), ',', int2str(subtotalPoints(dist_count, 2))));
    %         disp(strcat('plotting__', int2str(subtotalPoints(dist_count, 1)), ', ', int2str(subtotalPoints(dist_count, 2))));
        else
            total = (sum(cluster))/(clusterPoints-1);
%             disp(total);
            %if it's not within the array 
            if(isempty(subtotalPoints) || (sum(find(subtotalPoints(:,1)==total(2)))==0 &&sum(find(subtotalPoints(:,2)==total(3)))==0))
                dist_count = dist_count +1;
                subtotalPoints(dist_count, 1) = int64(total(2));
                subtotalPoints(dist_count, 2) = int64(total(3));
%                 disp('created point:');
%                 disp(strcat(int2str(subtotalPoints(dist_count, 1)), ',', int2str(subtotalPoints(dist_count, 2))));
            end
%             dist_count = dist_count +1;
%             subtotalPoints(dist_count, 1) = int64(total(2));
%             subtotalPoints(dist_count, 2) = int64(total(3));
%             disp('created point:');
%             disp(strcat(int2str(subtotalPoints(dist_count, 1)), ',', int2str(subtotalPoints(dist_count, 2))));
        end
    %     disp(strcat('Printing all similar points ',num2str(i)));
    %     disp(subtotalPoints);
    end
    
    
    % for i = 1:sz(1)
    %     alone = 1;
    %     for j = 1:sz(1)
    %         if(i~=j)
    %             dx = corners.Location(i,1) - corners.Location(j,1);
    %             dy = corners.Location(i,2) - corners.Location(j,2);
    %             dist = sqrt(dx*dx + dy*dy);
    %             
    %             distances(i,j) = dist;
    %             
    %             if(max_distance >= dist)
    % %                  disp(dist);
    %                  alone = 0;
    % %                 if(one_point <= amountOfPoints*2)
    %                     dist_count = dist_count + 1;
    %                     subtotalPoints(dist_count, 1) = int64((corners.Location(i,1) + corners.Location(j,1))/2);
    %                     subtotalPoints(dist_count, 2) = int64((corners.Location(i,2) + corners.Location(j,2))/2);
    % %                     disp(strcat(int2str(corners.Location(i,1)), ',', int2str(corners.Location(i,2)), '_and_', int2str(corners.Location(j,1)), ',', int2str(corners.Location(j,2)), '__average to__', int2str(subtotalPoints(dist_count, 1)), ', ', int2str(subtotalPoints(dist_count, 2))));
    %                     one_point = one_point+1;
    % %                 end
    %                 
    % % %                 test(i,j,1) = (corners.Location(i,1) + corners.Location(j,1))/2;
    % % %                 test(i,j,2) = (corners.Location(i,2) + corners.Location(j,2))/2;
    %                 
    %                 
    % 
    %             end
    %             
    %         end
    %         
    %     end
    %     
    %     if(alone == 1)
    %         dist_count = dist_count +1;
    %         subtotalPoints(dist_count, 1) = corners.Location(i,1);
    %         subtotalPoints(dist_count, 2) = corners.Location(i,2);
    % %         disp(strcat('plotting__', int2str(subtotalPoints(dist_count, 1)), ', ', int2str(subtotalPoints(dist_count, 2))));
    %     end
    %     
    % %     disp(strcat('Printing all similar points ',num2str(i)));
    % %     disp(subtotalPoints);
    % end
    % 
    %                     

    % datestr(now,'dd-mm-yyyy HH:MM:SS FFF')
    
    message = strcat('This iteration was number ', int2str(x));
    disp(message);
    
    disp('The time difference was');
    now = clock;
    disp(now-t0);
    disp('Now the points:');
    
    final = cornerPoints(subtotalPoints);
    corners = final;
    disp(length(final));
    % disp(final.Location);
    % disp(final.Location(1,2));
    
%     corners = final.selectStrongest(int64(length(subtotalPoints)*0.96));
%     %attempt to make the selection more accurate. failure. removed
     pause(1); % DEBUG
    sz = length(corners);
    imshow(I);
%     hold on
    plot(corners);
    
end
disp('The initial amount of points was');
disp(amountOfPoints);

