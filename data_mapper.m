

% load('web9_run199points.mat');



% I = imread('web9.jpg');
% imshow(I)
% hold on
% plot(true_points(:,1),true_points(:,2), 'c+','markerfacecolor','c');
% plot(Z(:,1),Z(:,2), 'mo', 'markerfacecolor' , 'm')
% plot(Z(:,1),Z(:,2), 'mo', 'markerfacecolor' , 'm')
% plot(Z(:,1),Z(:,2), 'mx', 'markerfacecolor' , 'm');


times6p1 = zeros(15,1)
times6 = zeros(15,1)
for spot = 1:15
    run('corner_improvement_rc1point6point1.m');
    dif = now-t0;
    times6p1(spot) = (dif(length(dif))) + 60*mod((dif(length(dif)-1)),60);
    run('corner_improvement_rc1point6.m');
    dif = now-t0;
    times6(spot) =  (dif(length(dif))) + 60*mod((dif(length(dif)-1)),60);%dif(length(dif))+(60*dif(length(dif)-1));% dif(length(dif));
    
end

avg6p1 = mean(times6p1)
avg6 = mean(times6)