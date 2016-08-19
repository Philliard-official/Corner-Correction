% clear all
% clc
close all
im = 14

filename = strcat('web',int2str(im),'true_points.mat');
OG = imread(strcat('web', int2str(im), '.jpg')); %read in image

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
% if(base<0.7 && base > 0.3)
%     I= imadjust(I, [base 0.80], []); % adjust contrast to improve results
%     
%     disp('Adjusted contrast');
% end
try
    load(filename);
catch ME
    true_points = [];
end
                        
imshow((I)) %display image to select regoin of interest (ROI)5
hold on % keep image shown even when plotting points
% zoom;
if(isempty(true_points)~=true)
    close all
    imshow((I)) %display image to select regoin of interest (ROI)5
    hold on % keep image shown even when plotting points
    plot(true_points(:,1), true_points(:,2), 'gx', 'markerfacecolor' , 'g');
end
% true_points = points;
start = length(true_points)+1;
i=start;
while i<start+10
    
    
    try
        [true_points(i,1), true_points(i,2)] = ginput(1);
        i = i+1;
        close all
        imshow((I)) %display image to select regoin of interest (ROI)5
        hold on % keep image shown even when plotting points
        plot(true_points(:,1), true_points(:,2), 'gx', 'markerfacecolor' , 'g');
    catch
        true_points = true_points(1:(length(true_points)-1),:);
        disp('undo');
        close all
        imshow((I)) %display image to select regoin of interest (ROI)5
        hold on % keep image shown even when plotting points
        plot(true_points(:,1), true_points(:,2), 'gx', 'markerfacecolor' , 'g');
    end

end
true_points = true_points((true_points(:,1)~=0),:);
% filename = strcat('web',int2str(im),'true_points.mat');
save(filename, 'true_points');


