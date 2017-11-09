% Spatial distribution estimation (SPE) using Bayesian Approach (MMLE + SIR filter)
% Hyongju Park
clear all;close all;clc
%% 

nsampInfo = 100;    % number of information sample
dSamp = 16;          % down-sample, e.g., 1/8, smaller the value, the more accurate

varInfo = 0.01;     % how noisy the sensor is (e.g., 0: perfect) default 0.01....
varPos = 0.01;     % decay as the distance between the windshield wiper measurement and the source of rain, increases...
weight = 0.5;		% weight on prior (e.g., if weight=1, it forgets about the prior measurements and only consider the new radar measurement, default: 0.9)

vehicleData = csvread('./data/20140811_vehicle_filtered.csv');
allVehicleID = unique(vehicleData(:,1))';
numExcludedVehicles = 1;
excludedVehicleID = allVehicleID(1,randsample(1:length(allVehicleID),numExcludedVehicles));
% excludedVehicleID = 10153;

% addpath(genpath('./gpml-matlab-v3.6-2015-07-07/'))


ncdisp('./data/data_20140811.nc');
gageData = ncread('./data/data_20140811.nc','gage');
radarData = ncread('./data/data_20140811.nc','radar');
wiperData = ncread('./data/data_20140811.nc','wiper');

%set radar data to zero if it's NaN
radarData(isnan(radarData)) = 0;
wiperData(isnan(wiperData)) = 0;
lonNet = ncread('./data/data_20140811.nc','longitude');
latNet = ncread('./data/data_20140811.nc','latitude');


%% load data
% ncdisp('/Users/mdbartos/Data/combined_nc/data_20140811.nc');
% gageData = ncread('/Users/mdbartos/Data/combined_nc/data_20140811.nc','gage');
% radarData = ncread('/Users/mdbartos/Data/combined_nc/data_20140811.nc','radar');
% wiperData = ncread('/Users/mdbartos/Data/combined_nc/data_20140811.nc','wiper');
% 
% %set radar data to zero if it's NaN
% radarData(isnan(radarData)) = 0;
% wiperData(isnan(wiperData)) = 0;
% lonNet = ncread('/Users/mdbartos/Data/combined_nc/data_20140811.nc','longitude');
% latNet = ncread('/Users/mdbartos/Data/combined_nc/data_20140811.nc','latitude');

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

%%
% clear M for later use...
clear M;
M = [qx(:) qy(:) qz(:)];
% let ground truth be 1 for all regions...
gTruth = ones(size(M,1),1);

% generate samples for information vector
hSet2 = haltonset(1,'Skip',1e3,'Leap',1e2);
hScrambled2 = scramble(hSet2,'RR2');

maxes = []

for i=1:nRuns
    maxes(i) = max(max(radarTSeries{i}));    
end

grandmax = max(maxes);
sampInfo = net(hScrambled2,nsampInfo) ;


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



%%
% ++++++++++++PLOT++++++++++++
[qx,qy,qz] = drawFigure(M,radarTSeries{1},dSamp);
% ++++++++++++PLOT++++++++++++

% save current data
savData{1,1} = sampPos;
savData{1,2} = particleWgt;
nRuns = max(unique(vehicleData(:,end-2)));


% nRuns = 100;
k5 = 0;
for curStep =1:nRuns
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
        M(:,3) = radarTSeries{curStep}(:)/grandmax;
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
    if exist('prvWgt0')
        for i = 1:size(sampPos,1)
            particleWgtRad{curStep}(:,i) = sampInfo' * prvWgt0(:,i);
        end
    else
        particleWgtRad{curStep} = savData{1,2};
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
		k1 = 0;
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
                        t_val = median(nearPosUniqWiper{i,1}(find(nearPosUniqWiper{i,2}==vehicleID(j)),:));
                        if isnan(t_val)
                            t_val = randsample([0 1 2],1);
                        end
                        nearPosMed(k,:) = [nearPosUniq(i) vehicleID(j) min(t_val,2)];
                        
                    end
                else
                    if ~isempty(find(nearPosUniqWiper{i,2}==vehicleID(j)))
                        invalid = 0;
                        k1 = k1+1;
                        t_val = median(nearPosUniqWiper{i,1}(find(nearPosUniqWiper{i,2}==vehicleID(j)),:));
                        if isnan(t_val)
                            t_val = randsample([0 1 2],1);
                        end
                        
                        nearPosMedEx(k1,:) = [nearPosUniq(i) vehicleID(j) min(t_val,2)];
                        
                    else
                        invalid = 1;
                    end                    
                end
            end
        end
        
        clear M;

        M(:,1) = qx(:);M(:,2) = qy(:);
		
		% find indice for unique vehicle positions with windshield-wiper measurements
        if exist('nearPosMed')
            rainIdx = nearPosMed(:,1);
        end

		maxRainIntensity = 1;		% average rain intensity that is captured by windshield wiper 	(level 2)
		% obtain rain detection likelihoods from the wind-shield wiper measurements
        for i_2 = 1:size(M,1)
            misDetectLhd = 1;
            if exist('nearPosMed')
            [out1, out2] = nearestPntDist(M(i_2,1:2),M(rainIdx,1:2));
            if out1 <= radD
%                 if ismember(i_2,nearPosMed(:,1))
                    i_1 = find(nearPosMed(:,1)==i_2);
                    if nearPosMed(i_1,3) == 1
                        gTruth(i_2) = max(radarTSeries{curStep}(i_2)/grandmax, 0.1/grandmax);
                    elseif nearPosMed(i_1,3) == 2
                        gTruth(i_2) = max(radarTSeries{curStep}(i_2)/grandmax, 0.2/grandmax);
                    elseif nearPosMed(i_1,3) == 3
                        gTruth(i_2) = max(radarTSeries{curStep}(i_2)/grandmax, 0.1/grandmax);
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
        M(:,3) =particleWgt'./max(sampInfo);    
        set(gca,'color','none');
        [qx,qy,qz] = drawFigure(M,radarTSeries{curStep},dSamp);
        hold on;
        plot(M(nearPosUniq,1),M(nearPosUniq,2),'cd');
        filestr = sprintf('./img/%03d.png', curStep);
%         saveas(gcf, filestr);
        clf;
        close;
    else
        clear M;
        M(:,1) = qx(:);
        M(:,2) = qy(:);        
        M(:,3) =particleWgt'./max(sampInfo);
        set(gca,'color','none');
        [qx,qy,qz] = drawFigure(M,radarTSeries{curStep},dSamp);
        hold on;
        plot(M(nearPosUniq,1),M(nearPosUniq,2),'x');
        filestr = sprintf('./img/%03d.png', curStep);
%         saveas(gcf, filestr);
        clf;
        close;
    end   
  
	% save data...
    savData{curStep+1,1} = sampPos;
    savData{curStep+1,2} = particleWgt;
    if exist('nearPosMedEx')
        tmpComp1 = particleWgt./max(sampInfo);
        tmpComp5 = particleWgtRad{curStep}./max(sampInfo);
        tmpComp3 = tmpComp1(:,nearPosMedEx(:,1));
        tmpComp6 = tmpComp5(:,nearPosMedEx(:,1));
%         if all(tmpComp6 > 0)
%             tmpComp4 = ones(size(tmpComp6));
%         else
%             tmpComp4 = zeros(size(tmpComp6));
%         end
%         if all(tmpComp3 > 0)
%             tmpComp7 = ones(size(tmpComp3));
%         else
%             tmpComp7 = zeros(size(tmpComp3));
%         end    
        tmpComp2 = nearPosMedEx(:,3)';

        if invalid == 0
            k5 = k5 + 1;
            correct1{k5} = tmpComp6; % only radar
            correct2{k5} = tmpComp3; % fused
            correct3{k5} = tmpComp2; % windshield wiper measurement of the car left out...
%             if (length(find(tmpComp2)) == length(find(tmpComp4)))
%                 correct1(k5) = 1;
%             else
%                 correct1(k5) = 0;
%             end
%             if (length(find(tmpComp2)) == length(find(tmpComp7)))
%                 correct2(k5) = 1;
%             else
%                 correct2(k5) = 0;
%             end        
            correct1
            correct2
            correct3
        end
    end
end