% Spatial distribution estimation (SPE) using Bayesian Approach (MMLE + SIR filter)
% Hyongju Park
clear all;close all;clc
%% 

nsampInfo = 100;    % number of information sample
dSamp = 16;          % down-sample, e.g., 1/8, smaller the value, the more accurate

varInfo = 0.01;     % how noisy the sensor is (e.g., 0: perfect) default 0.01....
varPos = 0.01;     % decay as the distance between the windshield wiper measurement and the source of rain, increases...
weight = 0.9;		% weight on prior (e.g., if weight=1, it forgets about the prior measurements and only consider the new radar measurement, default: 0.9)

vehicleData = csvread('./data/20140811_vehicle_filtered.csv');
allVehicleID = unique(vehicleData(:,1))';
numExcludedVehicles = 1;
excludedVehicleID = allVehicleID(1,randsample(1:length(allVehicleID),numExcludedVehicles));


% addpath(genpath('./gpml-matlab-v3.6-2015-07-07/'))

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

% change scale (GPS locations -> [0, 1]x[0, 1]) % for the sake of convenience
lonNetScaled=(lonNet-min(lonNet))/(max(lonNet) - min(lonNet));
latNetScaled=(latNet-min(latNet))/(max(latNet) - min(latNet));

%radius of rain detection for each vehicle (scaled based upon our scaled area [0, 1]x[0, 1])
radD = 0.1;


% generate timeseries data {radar, windshield wiper}
for i = 1:size(radarData,3)
    radarTSeries{i} = radarData(:,:,i);
    wiperTSeries{i} = wiperData(:,:,i);    
end

% extract effective radar measurements (with non-NaNs)
radarIdx = [];
for i = 1:length(radarTSeries)
    if ~isempty(find(radarTSeries{i}, 1))
        radarIdx(end+1) = i;
    end
end
for i = 1:length(radarIdx)
    radar_nz(i) = length(find(radarTSeries{radarIdx(i)}));
end
for i = 1:length(radar_nz)
    radar_nz2(i) = radarIdx(i);
end
%===============NO NEED TO CHANGE ABOVE===========================

% number of runs for a finite time horizeon
nRuns = max(unique(vehicleData(:,end-2)));

% generate matrix M to plot the radar measurements
M(:,1) = repmat(lonNetScaled,size(radarTSeries{1},2),1);
lvec = [];
for i= 1:size(radarTSeries{1},2)
    lvec = [lvec;repmat(latNetScaled(i),size(radarTSeries{1},1),1)];
end
M(:,2) = lvec;
% % normalize
M(:,3) = radarTSeries{1}(:);
% 
% % ++++++++++++PLOT++++++++++++
[qx,qy,qz] = drawNoFigure(M,radarTSeries{1},dSamp);
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
% prvWgt0 = [];
% for i = 1:size(M,1)
%     p0 = mvnpdf(sampInfo,M(i,3),varInfo)/sum(mvnpdf(sampInfo,M(i,3),varInfo)); 
%     prvWgt0 = [prvWgt0 p0];
% end
prvWgt1 = ones(size(sampInfo,1),size(M,1))/size(sampInfo,1);
prvWgt = prvWgt1;

% generate samples for locations (use M)
sampPos = M(:,1:2);

% expected value for the information, e.g., rain intensity...
particleWgt = sampInfo'*prvWgt;

% update M using the weighted prior
M(:,3) =particleWgt'/max(sampInfo);


% ++++++++++++PLOT++++++++++++
[qx,qy,qz] = drawNoFigure(M,radarTSeries{1},dSamp);
% ++++++++++++PLOT++++++++++++

% save current data
savData{1,1} = sampPos;
savData{1,2} = particleWgt;

for curStep = 1:nRuns
    curStep
    
	% if there is a non-zero, non-NaN radar measurement, we will update the previous belief...
    if ismember(curStep,radar_nz2)
        clear M;
        % DO SENSOR FUSION...
        % generate matrix M to plot the radar measurements
        M(:,1) = repmat(lonNetScaled,size(radarTSeries{curStep},2),1);
        lvec = [];
        for i= 1:size(radarTSeries{curStep},2)
            lvec = [lvec;repmat(latNetScaled(i),size(radarTSeries{curStep},1),1)];
        end
        M(:,2) = lvec;
        % normalize
        M(:,3) = radarTSeries{curStep}(:)/max(radarTSeries{curStep}(:));
        [qx,qy,qz] = drawNoFigure(M,radarTSeries{curStep},dSamp);
        clear M;
        M = [qx(:) qy(:) qz(:)];
        prvWgt0 = [];
        for i = 1:size(M,1)
            p0 = mvnpdf(sampInfo,M(i,3),varInfo)/sum(mvnpdf(sampInfo,M(i,3),varInfo)); 
            prvWgt0 = [prvWgt0 p0];
        end
		% current belief = weight * (PDF from radar measurement) * (1 - weight) * prior belief
        prvWgt = weight*prvWgt0 + (1-weight)*prvWgt;         
    end
	% find time frame where there is windshield-wiper measurement 
    wiperOnIdx = find(vehicleData(:,7) == curStep);
    if ~isempty(wiperOnIdx)
		% find unique vehicle IDs
        vehicleID = unique(vehicleData(:,1));
        wiperInfo = [];        
        vehicleDataEff = vehicleData(wiperOnIdx,:);
		% scale GPS locations to fit [0, 1] x [0, 1]
        vehicleDataEff(:,4)=(vehicleDataEff(:,4)-min(lonNet))/(max(lonNet) - min(lonNet));
        vehicleDataEff(:,3)=(vehicleDataEff(:,3)-min(latNet))/(max(latNet) - min(latNet));        
		% find the closest points from the downsampled map to vehicle locations
        clear nearPos;
        for i = 1:size(vehicleDataEff,1)  
            [~,nearPos(i)] = nearestPntDist([vehicleDataEff(i,4) vehicleDataEff(i,3)],[qx(:),qy(:)]);
        end

		% [winshield wiper data, vehicle positions on the map]
        vehicleDataEff = [vehicleDataEff nearPos'];
		
		% find unique vehicle positions
        [nearPosUniq, idv] = unique(nearPos);
        clear nearPosMed;
        k = 0;
		
        for i = 1:length(nearPosUniq)
			% find windshield wiper value(5), and vehicleID(1) at unique vehicle positions...
            nearPosUniqWiper{i,1}=vehicleDataEff(find(vehicleDataEff(:,end) == nearPosUniq(i)),5);
            nearPosUniqWiper{i,2}=vehicleDataEff(find(vehicleDataEff(:,end) == nearPosUniq(i)),1);

			% for each vehicle (with different IDs)
            for j = 1:length(vehicleID)
                % if the vehicle is not excluded from the experiment...
                if ~ismember(vehicleID(j),excludedVehicleID)
                    % take median of the windshield-wiper measurements made by each vehicle
                    if ~isempty(find(nearPosUniqWiper{i,2}==vehicleID(j)))
                        k = k+1;
                        nearPosMed(k,:) = [nearPosUniq(i) vehicleID(j) min(median(nearPosUniqWiper{i,1}(find(nearPosUniqWiper{i,2}==vehicleID(j)),:)),2)];
                    end
                end
            end
        end
        clear M;

        M(:,1) = qx(:);M(:,2) = qy(:);
		
		% find indice for unique vehicle positions with windshield-wiper measurements
        if exist('nearPosMed') %#ok<EXIST>
            rainIdx = nearPosMed(:,1);
        end

		maxRainIntensity = 1;		% average rain intensity that is captured by windshield wiper 	(level 2)
		% obtain rain detection likelihoods from the wind-shield wiper measurements
        for i_2 = 1:size(M,1)
            misDetectLhd = 1;
            if exist('nearPosMed') %#ok<EXIST>
                [out1, out2] = nearestPntDist(M(i_2,1:2),M(rainIdx,1:2));
                if out1 <= radD
    %                 if ismember(i_2,nearPosMed(:,1))
                        i_1 = find(nearPosMed(:,1)==i_2);
                        if nearPosMed(i_1,3) == 1
                            gTruth(i_2) = maxRainIntensity/2;
                        elseif nearPosMed(i_1,3) == 2
                            gTruth(i_2) = maxRainIntensity;
                        else
                            gTruth(i_2) = 0;
                        end
    %                 end
                    misDetectLhd = misDetectLhd * (1-mvnpdf(M(i_2,1:2),M(rainIdx(out2),1:2),eye(2)*varPos)/mvnpdf([0 0],[0 0],eye(2)*varPos));
                end
            end
            detectLhd(i_2) = 1-misDetectLhd;
            
        end

		% obtain observation likelihoods
        for i = 1:size(sampPos,1)
            for j = 1:length(sampInfo)
                  infoLhd{i}(j)=normpdf(sampInfo(j),gTruth(i,:),max(gTruth)*varInfo)/normpdf(0,0,max(gTruth)*varInfo);
            end

            infoLhd{i} = (infoLhd{i}-1/max(gTruth)) * detectLhd(i)+ 1/max(gTruth);
        end
		
		% particle filtering
        for i = 1:size(sampPos,1)
            wgt2(:,i) = sirFilter(infoLhd{i}',prvWgt(:,i));
            prvWgt(:,i) = wgt2(:,i);                
        end

        for i = 1:size(sampPos,1)
            particleWgt(:,i) = sampInfo' * wgt2(:,i);
        end
        clear M;        
        M(:,1) = qx(:);
        M(:,2) = qy(:);        
        M(:,3) =particleWgt'/max(sampInfo);        
        [qx,qy,qz] = drawNoFigure(M,radarTSeries{curStep},dSamp);
%         hold on;
%         plot(M(nearPosUniq,1),M(nearPosUniq,2),'x');        
    else
        clear M;
        M(:,1) = qx(:);
        M(:,2) = qy(:);        
        M(:,3) =particleWgt'/max(sampInfo);
        [qx,qy,qz] = drawNoFigure(M,radarTSeries{curStep},dSamp);
%         hold on;        
    end   
  
	% save data...
    savData{curStep,1} = sampPos;
    savData{curStep,2} = particleWgt;
end


time = 25;
fig2 = figure('position',[100 100 1600 400],'Color',[1 1 1]);

load('radarOnly.mat');
M(:,1) = qx(:);
M(:,2) = qy(:);    
particleWgt = savData{time,2};
M(:,3) =particleWgt'/max(sampInfo);
subplot(1,3,1)
[qx,qy,qz] = drawsubFigure(M,radarTSeries{time},dSamp);
hold on; title('radar only');
load('radarWiper.mat');
M(:,1) = qx(:);
M(:,2) = qy(:);        
particleWgt = savData{time,2};
M(:,3) =particleWgt'/max(sampInfo);
subplot(1,3,2)
[qx,qy,qz] = drawsubFigure(M,radarTSeries{time},dSamp);
hold on; title('radar + windshield wiper');
load('radarWiperEx3.mat');
M(:,1) = qx(:);
M(:,2) = qy(:);        
particleWgt = savData{time,2};
M(:,3) =particleWgt'/max(sampInfo);
subplot(1,3,3)
[qx,qy,qz] = drawsubFigure(M,radarTSeries{time},dSamp);
hold on; title('radar + windshield wiper - 3 vehicles');
suptitle(sprintf('time = %d',time))
% colorbar;






