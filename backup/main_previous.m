% Spatial distribution estimation (SPE) using Bayesian Approach (MMLE + FIR filter)
% Hyongju Park
clear all;close all;clc
%% 

nsampInfo = 100;    % number of information sample
dSamp = 8;          % down-sample, e.g., 1/8

varInfo = 0.01;     % how noisy the sensor is (e.g., 0: perfect) default 0.5....
varPos = 0.01;     % decay as the distance between the windshield wiper measurement and the source of rain, increases...

vehicleData = csvread('./data/20140811_vehicle_filtered.csv');

addpath(genpath('./gpml-matlab-v3.6-2015-07-07/'))

%% load data
ncdisp('./data/data_20140811.nc');
gageData = ncread('./data/data_20140811.nc','gage');
radarData = ncread('./data/data_20140811.nc','radar');
wiperData = ncread('./data/data_20140811.nc','wiper');

%set radar data to zero if it's NaN
radarData(isnan(radarData)) = 0;
wiperData(isnan(wiperData)) = 0;
lonNet = ncread('./data/data_20140811.nc','longitude');
latNet = ncread('./data/data_20140811.nc','latitude');

%radius of detection
radD = 0.1;

% change scale (GPS locations -> [0 1 0 1]) % for the sake of convenience
lonNetScaled=(lonNet-min(lonNet))/(max(lonNet) - min(lonNet));
latNetScaled=(latNet-min(latNet))/(max(latNet) - min(latNet));

% generate timeseries data {radar, windshield wiper}
for i = 1:size(radarData,3)
    radarTSeries{i} = radarData(:,:,i);
    wiperTSeries{i} = wiperData(:,:,i);    
end

radarIdx = [];
for i = 1:length(radarTSeries)
    if ~isempty(find(radarTSeries{i}, 1))
        radarIdx(end+1) = i;
    end
end
for i = 1:length(radarIdx)
    radar_nz(i) = length(find(radarTSeries{radarIdx(i)}));
end


% some sample time where there is radar measurements..., e.g., 127
k1 = 127;

% find non-empty wiperData
wiperIdx = [];
for i = 1:length(wiperTSeries)
    if ~isempty(find(wiperTSeries{i}, 1))  
        if i >= k1
            wiperIdx(end+1) = i;
        end
        
    end
end
nRuns = length(wiperIdx);  % total number of runs

% generate matrix M to plot the radar measurements
M(:,1) = repmat(lonNetScaled,size(radarTSeries{k1},2),1);
lvec = [];
for i= 1:size(radarTSeries{k1},2)
    lvec = [lvec;repmat(latNetScaled(i),size(radarTSeries{k1},1),1)];
end
M(:,2) = lvec;
% normalize
M(:,3) = radarTSeries{k1}(:)/max(radarTSeries{k1}(:));

% ++++++++++++PLOT++++++++++++
[qx,qy,qz] = drawFigure(M,radarTSeries{k1},dSamp);
% ++++++++++++PLOT++++++++++++

% clear M for later use...
clear M;
M = [qx(:) qy(:) qz(:)];
% let ground truth be 1 for all regions...
gTruth = ones(size(M,1),1);

% generate samples for information vector
hSet2 = haltonset(1,'Skip',1e3,'Leap',1e2);
hScrambled2 = scramble(hSet2,'RR2');
sampInfo = net(hScrambled2,nsampInfo)*max(gTruth);


% apply different weights to priors:
%   - prvWgt0: prior from radar measurements
%   - prvWgt1: uniform prior (no information)
prvWgt0 = [];
for i = 1:size(M,1)
    p0 = mvnpdf(sampInfo,M(i,3),varInfo)/sum(mvnpdf(sampInfo,M(i,3),varInfo)); 
    prvWgt0 = [prvWgt0 p0];
end
prvWgt1 = ones(size(prvWgt0))/size(prvWgt0,1);
prvWgt = 0.9*prvWgt0 + 0.1*prvWgt1;

% generate samples for locations (use M)
sampPos = M(:,1:2);

% expected value for the information, e.g., rain intensity...
particleWgt = sampInfo'*prvWgt;

% update M using the weighted prior
M(:,3) =particleWgt'/max(sampInfo);


% ++++++++++++PLOT++++++++++++
[qx,qy,qz] = drawFigure(M,radarTSeries{k1},dSamp);
% ++++++++++++PLOT++++++++++++

% save current data
savData{1,1} = sampPos;
savData{1,2} = particleWgt;

for curStep = 2:nRuns
    curStep
    % clear previous M to reuse the variable
    clear M;
    M(:,1) = repmat(lonNetScaled,size(wiperTSeries{wiperIdx(curStep-1)},2),1);
    lvec = [];
    for i= 1:size(wiperTSeries{wiperIdx(curStep-1)},2)
        lvec = [lvec;repmat(latNetScaled(i),size(wiperTSeries{wiperIdx(curStep-1)},1),1)];
    end
    M(:,2) = lvec;
    M(:,3) = wiperTSeries{wiperIdx(curStep-1)}(:);
    
    rainIdx = find(M(:,3));
    rainIdxNew = [];
    
    % extend precipitation data (to overcome possible data loss due to
    % downsampling)
    if ~isempty(rainIdx)
        for i = 1:size(M,1)
            if nearestPntDist(M(i,1:2),M(rainIdx,1:2))< 0.01 && ~ismember(i,rainIdx)
                rainIdxNew = [rainIdxNew i];
            end
        end
        rainIdxNew = [rainIdxNew rainIdx'];
    else
        rainIdxNew = rainIdx;
    end
    clear rainIdx; rainIdx = rainIdxNew';
    for i = 1:size(M,1)
        if ismember(i, rainIdx)
            M(i,3) = 1;
        end
    end
    F = scatteredInterpolant(M(:,1), M(:,2), M(:,3));
    [qx, qy] = meshgrid(linspace(min(M(:,1)), max(M(:,1)), ceil(size(wiperTSeries{wiperIdx(curStep-1)},2)/dSamp)), ...
                    linspace(min(M(:,2)), max(M(:,2)), ceil(size(wiperTSeries{wiperIdx(curStep-1)},1)/dSamp)));
    qz = F(qx, qy);
    detectLhd = qz(:); 
    clear M;
    M(:,1) = qx(:);M(:,2) = qy(:);M(:,3) = qz(:);
    clear rainIdx; rainIdx= find(M(:,3));
    % find detections...
%     tmp = find(detectLhd);
    for i_2 = 1:size(M,1)
        misDetectLhd = 1;
        [out1, out2] = nearestPntDist(M(i_2,1:2),M(rainIdx,1:2));
        if out1 <= radD
            misDetectLhd = misDetectLhd * (1-mvnpdf(M(i_2,1:2),M(rainIdx(out2),1:2),eye(2)*varPos)/mvnpdf([0 0],[0 0],eye(2)*varPos));
        end      
        detectLhd(i_2) = 1-misDetectLhd;
    end
%     clear qx qy qz M

    for i = 1:size(sampPos,1)
        for j = 1:length(sampInfo)
              infoLhd{i}(j)=normpdf(sampInfo(j),gTruth(i,:),max(gTruth)*varInfo);
        end
 
        infoLhd{i} = (infoLhd{i}-1/max(gTruth)) * detectLhd(i)+ 1/max(gTruth);
    end

    for i = 1:size(sampPos,1)
        wgt2(:,i) = sirFilter(infoLhd{i}',prvWgt(:,i));
        prvWgt(:,i) = wgt2(:,i);                
    end

    for i = 1:size(sampPos,1)
        particleWgt(:,i) = sampInfo' * wgt2(:,i);
    end
    
%     particleWgt = wgt2Avg/sum(wgt2Avg);
    clear M;
    M(:,1) = qx(:);
    M(:,2) = qy(:);
    % normalize
    M(:,3) =particleWgt'/max(sampInfo);
    % ++++++++++++PLOT++++++++++++
    [qx,qy,qz] = drawFigure(M,radarTSeries{k1},dSamp);
    % ++++++++++++PLOT++++++++++++
    
    % save current data    
    savData{curStep,1} = sampPos;
    savData{curStep,2} = particleWgt;
end



