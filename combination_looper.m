counter = 1;
% dk = 75:85;
% dk = dk(mod((dk),5) == 0);
% 
% bd = 30:40;
% bd = double(bd(mod((bd),5) == 0))/10;
% sd = 60:70;
% sd = double(sd(mod((sd),5) == 0))/10;
% br = 11:13;
% br = double(br)/1000;
% sr = 16:18;
% sr = double(sr)/1000;
% pp = 85:95;
% pp = double(pp(mod((pp),5) == 0))/100;
% null = zeros(length(br)*length(sr));
% null = null(:,1);
% positives = null;
% false_negatives = null;
% false_positives = null;
% len_points = null;
% im = 13;
% load(strcat('web', int2str(im), 'data_extended2.mat'));
% close all

% if dk exists replace 1 with length(dk)


total_data = zeros(length(1*length(br)*length(sr)*length(bd)*length(pp)*length(sd)), 10);
% for dark_value = dk
    for percentPoints = pp
        for base_distance = bd
            for scaling_distance = sd
                for base_radius = br
                    for scaling_radius = sr
%                         total_data(counter, 1) = dark_value;
                        total_data(counter, 2) = percentPoints;
                        total_data(counter, 3) = base_distance;
                        total_data(counter, 4) = scaling_distance;
                        total_data(counter, 5) = base_radius;
                        total_data(counter, 6) = scaling_radius;
                        total_data(counter, 7) = positives(counter);
                        total_data(counter, 8) = false_positives(counter);
                        total_data(counter, 9) = false_negatives(counter);
                        total_data(counter, 10) = len_points(counter);
%                         total_data(counter, 11) = 
                        counter = counter+1;
                    end
                end
            end
        end
    end
% end
disp(total_data)
avg_positives = zeros(length(pp),1);
avg_precision = zeros(length(pp),1);
avg_sensitivity = zeros(length(pp),1);
avg_false_positives = zeros(length(pp),1);
avg_false_negatives = zeros(length(pp),1);

max_positives = find(positives==max(positives))
min_positives = find(positives==min(positives))


max_false_positives = find(false_positives==max(false_positives))
min_false_positives = find(false_positives==min(false_positives))


max_false_negatives = find(false_negatives==max(false_negatives))
min_false_negatives = find(false_negatives==min(false_negatives))

precision = positives;
for j = 1:length(precision)
    precision(j) = positives(j)/(positives(j)+false_positives(j));
end
max_precision = find(precision==max(precision))
min_precision =  find(precision==min(precision))

% sensitivity = positives;
% for j = 1:length(sensitivity)
%     sensitivity(j) = positives(j)/length(Z);
% end
% max_sensitivity = find(sensitivity==max(sensitivity))
% min_sensitivity = find(sensitivity==min(sensitivity))


false_discovery_rate = 1-precision;
max_FDR = find(false_discovery_rate==max(false_discovery_rate))
min_FDR = find(false_discovery_rate==min(false_discovery_rate))


% % false_negative_rate = 1-sensitivity;
% max_FNR = find(false_negative_rate==max(false_negative_rate))
% min_FNR = find(false_negative_rate==min(false_negative_rate))

best_rows = union(union(union(max_precision,max_positives),min_false_positives),min_false_negatives)
% best_rows = (union(union(union(union(min_FDR,min_FNR), union(max_precision,max_sensitivity)), max_positives), min_false_positives))
best_data = total_data(best_rows,:);

best_dk = best_data == mode(best_data(:,1));
best_dk = best_dk(:,1);
best_dk = zeros(size(best_dk));
% best_dk = best_data(best_dk(:,1),:);
best_pp = best_data == mode(best_data(:,2));
best_pp = best_pp(:,2);
% best_pp = best_data(best_data == mode(best_data(:,2)),:)
best_bd = best_data == mode(best_data(:,3));
best_bd = best_bd(:,3);
% best_bd = best_data(best_data == mode(best_data(:,3)),:)
best_sd = best_data == mode(best_data(:,4));
best_sd = best_sd(:,4);
% best_sd = best_data(best_data == mode(best_data(:,4)),:)
best_br = best_data == mode(best_data(:,5));
best_br = best_br(:,5);
% best_br = best_data(best_data == mode(best_data(:,5)),:)
best_sr = best_data == mode(best_data(:,6));
best_sr = best_sr(:,6);
% best_sr = best_data(best_data == mode(best_data(:,6)),:)

bestidx = best_dk + best_pp + best_bd + best_sd + best_br + best_sr
best = best_data((bestidx == max(bestidx)),:);

for b=1:length(sr)
    avg_positives(b) = mean(total_data(total_data(:,6)==sr(b),7))
    avg_precision(b) = mean(precision(total_data(:,6)==sr(b)))
    avg_false_positives(b) = mean(total_data(total_data(:,6)==sr(b),8))
    avg_false_negatives(b) = mean(total_data(total_data(:,6)==sr(b),9));
%     avg_sensitivity(b) = mean(sensitivity(total_data(:,2)==pp(b)))
%     avg(b) = sum(double(total_data(total_data(:,5)==br(b),7)))/sum(total_data(total_data(:,1)==br(1)));
end

hold on
plot(sr,avg_positives, 'bo', 'markerfacecolor', 'c')
plot(sr,avg_positives, 'b--')
xlabel('Scaling Radius')
ylabel('Average Positive Values')
title('Scaling Radius of Point Grouping vs Confirmed Positive Points (Release Candidate 1.6.3, Web 13)')

pause(5)
close all
hold on

% plot(pp,avg_sensitivity, 'ko', 'markerfacecolor', 'r')
% plot(pp,avg_sensitivity, 'k--')
% xlabel('Base Distance')
% ylabel('Average Sensitivity')
% title('Base Distance of Subtractive Clustering vs Average Precision (Release Candidate 1.6)')


plot(sr,avg_precision, 'bo', 'markerfacecolor', 'm')
plot(sr,avg_precision, 'b--')
xlabel('Base Radius')
ylabel('Average Precision')
title('Base Radius for Subtractive Clustering vs Average Precision (Release Candidate 1.6.3, Web 13)')

pause(5)
close all
hold on

plot(sr,avg_false_positives, 'bo', 'markerfacecolor', 'm')
plot(sr,avg_false_positives, 'b--')
xlabel('Scaling Radius')
ylabel('Average False Positives')
title('Scaling Radius for Subtractive Clustering vs Average False Positives (Release Candidate 1.6.3, Web 13)')

pause(5)
close all
hold on

plot(sr,avg_false_negatives, 'bo', 'markerfacecolor', 'm')
plot(sr,avg_false_negatives, 'b--')
xlabel('Base Radius')
ylabel('Average False Negatives')
title('Base Radius for Subtractive Clustering vs Average False Negatives (Release Candidate 1.6.3, Web 13)')

