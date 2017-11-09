%% THIS SCRIPT LOADS DATA AND GENERATE ROC CURVES...

% LOAD DATA
clear all;close all;clc
% numRuns = 15;
% for l1 = 1:numRuns
%     save('num.mat','l1');
%     main
%     
%     if exist('correct1')
%         load('num.mat');
%         save(strcat('output9_',num2str(l1),'.mat'),'correct1','correct2','correct3');
%     end
% end

weight = 0.1; % 0.1, 0.5, 0.9

maxRuns = 100;
q1 = 0;
for q0 = 1:maxRuns
    if (exist(strcat('./result/wt_',num2str(weight),'_',num2str(q0),'.mat'),'file'))
        load(strcat('./result/wt_',num2str(weight),'_',num2str(q0),'.mat'));
        q1 = q1 + 1;
        for q2 = 1:length(correct1)
            gTruth{q1}(q2) = round(length(find(correct3{q2}))/length(correct3{q2}));
            data1{q1}(q2) = (nanmean(correct1{q2}));
            data2{q1}(q2) = (nanmean(correct2{q2}));
        end        
        
    end

end
genROC(data1, data2, gTruth, weight)



weight = 0.5; % 0.1, 0.5, 0.9

maxRuns = 100;
q1 = 0;
for q0 = 1:maxRuns
    if (exist(strcat('./result/wt_',num2str(weight),'_',num2str(q0),'.mat'),'file'))
        load(strcat('./result/wt_',num2str(weight),'_',num2str(q0),'.mat'));
        q1 = q1 + 1;
        for q2 = 1:length(correct1)
            gTruth{q1}(q2) = round(length(find(correct3{q2}))/length(correct3{q2}));
            data1{q1}(q2) = (nanmean(correct1{q2}));
            data2{q1}(q2) = (nanmean(correct2{q2}));
        end        
        
    end

end
genROC(data1, data2, gTruth, weight)

weight = 0.9; % 0.1, 0.5, 0.9

maxRuns = 100;
q1 = 0;
for q0 = 1:maxRuns
    if (exist(strcat('./result/wt_',num2str(weight),'_',num2str(q0),'.mat'),'file'))
        load(strcat('./result/wt_',num2str(weight),'_',num2str(q0),'.mat'));
        q1 = q1 + 1;
        for q2 = 1:length(correct1)
            gTruth{q1}(q2) = round(length(find(correct3{q2}))/length(correct3{q2}));
            data1{q1}(q2) = (nanmean(correct1{q2}));
            data2{q1}(q2) = (nanmean(correct2{q2}));
        end        
        
    end

end
genROC(data1, data2, gTruth, weight)
