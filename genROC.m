function f0 = genROC(data1, data2, gTruth, weight)

R1 = [];R2 = [];R3 = [];
nData = length(data1);
for q1 = 1:nData
    R1 = [R1 data1{q1}];
    R2 = [R2 data2{q1}];
    R3 = [R3 gTruth{q1}];
end
R1 = (R1-min(R1))./(max(R1)-min(R1));
R2 = (R2-min(R2))./(max(R2)-min(R2));

GT = R3;
P = nnz(GT); % number of positive responses in ground truth
N = nnz(1-GT);
% responses
% R = imread('../../Downloads/matlab/data.tif');
% your thresholds
% thresholds = [0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8];
% alternatively, use 100 thresholds between min(R) and max(R)
thresholds = linspace(min(R1(:)), max(R1(:)));
% pre-allocate for speed
tp = nan(1, length(thresholds));
fp = nan(1, length(thresholds));
for i = 1:numel(thresholds)
  t = thresholds(end-i+1); % thresholds from high to low as i increases
  Rt = R1 > t; % thresholded response
  tp(i) = nnz(Rt & GT);
  fp(i) = nnz(Rt & ~GT);
end
% convert to rates
TPR = tp/P;
FPR = fp/N;
f0 = figure('position',[100 100 800 800],'Color',[1 1 1]);
plot(FPR, TPR,'k--') % ROC
hold on;
thresholds = linspace(min(R2(:)), max(R2(:)));
% pre-allocate for speed
tp = nan(1, length(thresholds));
fp = nan(1, length(thresholds));
for i = 1:numel(thresholds)
  t = thresholds(end-i+1); % thresholds from high to low as i increases
  Rt = R2 > t; % thresholded response
  tp(i) = nnz(Rt & GT);
  fp(i) = nnz(Rt & ~GT);
end
% convert to rates
clear TPR FPR
TPR = tp/P;
FPR = fp/N;
plot(FPR, TPR,'k-'); hold on; % ROC
line([0 1],[0 1],'LineStyle','-.','Color','r');
legend('radar only','ww+radar');
ylabel('TPR (sensitivity)');
xlabel('FPR (1-specificity)');
axis('equal');
axis([0 1 0 1]);
title(sprintf('prior weight on radar: %.2f',weight));
set(gca,'FontSize',12);