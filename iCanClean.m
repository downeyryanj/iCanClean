% iCanClean() - perform EEG artifact removal using iCanClean algorithm.
%
% Usage:
%   >>  cleanEEG = iCanClean( EEG, xChan, yChan, visualizeFinalResults, params)
%
% Inputs:
%   EEG     - input EEG dataset
%   xChan   - indices of X channels (Data)
%   yChan   - indices of Y channels (Reference Noise)
%   visualizeFinalResults - option to turn on/off final results
%                           visualization
%   params  - parameter structure containing additional options
%
% Outputs:
%   cleanEEG     - cleaned EEG dataset
%
% See also:
%   pop_iCanClean, EEGLAB

% Copyright (C) <2020-2023>  <Ryan Downey>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% Please note that this software implements patented methods as described in:
%
% "Removing Latent Noise Components from Data Signals," patent pending, 
% non-provisional patent application submitted to United States Patent and 
% Trademark Office on August 25, 2021, Application No. PCT/US21/71283. 
% Published under International Publication No. WO 2022/061322 on March 24,
%  2022. https://patents.google.com/patent/WO2022061322A1/
%
% and in 
%
% “Using Pseudo Reference Noise Signals to Remove Latent Noise From 
% Data Signals and Identify Data Sources,” sole inventor, patent pending, 
% provisional patent application submitted to United States Patent and 
% Trademark Office on September 5, 2023, Serial No. 63/580,664.
%
% For licensing information, please contact the Univeristy of Florida's
% Tech Licensing center. Contact information below.
%
%
% Physical Address:
% UF Innovate | Tech Licensing
% 747 Southwest 2nd Avenue, Suite 108
% Gainesville, FL 32601
%
% Mailing Address:
% UF Innovate | Tech Licensing
% P.O. Box 115575
% Gainesville, FL 32611-5575
%
% Phone Number:
% (352) 392-8929
%
% Website
% https://innovate.research.ufl.edu/tech-licensing/
%
%
% To cite iCanClean, please reference the following paper:
%
% Downey, R.J.; Ferris, D.P. iCanClean Removes Motion, Muscle, Eye, 
% and Line-Noise Artifacts from Phantom EEG. Sensors 2023, 23, 8214. 
% https://doi.org/10.3390/s23198214
%

function cleanEEG = iCanClean( EEG, xChan, yChan, visualizeFinalResults, params)

if nargin < 1
    help iCanClean;
    return;
end

%% check channel indices validity
if     isempty(xChan), error('No data channels selected.');
elseif isempty(yChan), error('Cannot clean data without reference noise channels.');
end

if any(xChan < 1) || any(xChan > size(EEG.data,1))
    error('Invalid data channels (indices are out of range).');
elseif any(yChan < 1) || any(yChan > size(EEG.data,1))
    error('Invalid reference noise channels (indices are out of range).');
end

%% (temporary) use 3DEpochedData if EEG.data is 3D
if(ndims(EEG.data)==3)
    warning('The data appears to be epoched, so we are running a different script.')
    %     cleanEEG = iCanClean_3DEpochedData( EEG, xChan, yChan, visualizeFinalResults, params);
    cleanEEG = iCanClean_3DEpochedDataParallel( EEG, xChan, yChan, visualizeFinalResults, params);
    return;
end

%% initialize EEG structure that will eventually be our clean data and keep track of time
cleanEEG = EEG;
tic;

%% check/set parameters
params = iCanClean_setParams(params);
disp(params);

%% define temporary re-referenced and filtered data
%ryan new verion 2022-05-06
X_orig = EEG.data(xChan,:); %original data
Y_orig = EEG.data(yChan,:);

%X
switch params.filtXtype
    case 'no',   tempFiltRerefEEG = EEG; %no filtering
    case 'LP',   tempFiltRerefEEG = pop_eegfiltnew(EEG,'hicutoff',params.filtXfreq,'revfilt',0,'channels',xChan); %LP (motion)
    case 'HP',   tempFiltRerefEEG = pop_eegfiltnew(EEG,'locutoff',params.filtXfreq,'revfilt',0,'channels',xChan); %HP (muscle)
    case 'BP',   tempFiltRerefEEG = pop_eegfiltnew(EEG,'locutoff',params.filtXfreq(1),'hicutoff',params.filtXfreq(2),'revfilt',0,'channels',xChan); %BP (line noise)
    case 'Notch',tempFiltRerefEEG = pop_eegfiltnew(EEG,'locutoff',params.filtXfreq(1),'hicutoff',params.filtXfreq(2),'revfilt',1,'channels',xChan); %notch (all artifacts)
    otherwise,   error('Detected invalid filtering option for X.');
end
switch params.rerefX
    case 'no'
        X_temp = tempFiltRerefEEG.data(xChan,:); %no re-ref
    case 'yes-temp'
        X_temp = iCanClean_reref(tempFiltRerefEEG.data(xChan,:),params.rankX); %temp define a rereferenced version of X (to help CCA later)
    case {'yes-fullrank','yes-loserank'} %old, for backwards compatibility
        X_temp = iCanClean_reref(tempFiltRerefEEG.data(xChan,:),params.rerefX(5:end)); %temp define a rereferenced version of X (to help CCA later)
    otherwise
        error('Detected invalid rereferencing option for X.');    
end

%Y
switch params.filtYtype
    case 'no',   tempFiltRerefEEG = EEG; %no filtering
    case 'LP',   tempFiltRerefEEG = pop_eegfiltnew(EEG,'hicutoff',params.filtYfreq,'revfilt',0,'channels',yChan); %LP (motion)
    case 'HP',   tempFiltRerefEEG = pop_eegfiltnew(EEG,'locutoff',params.filtYfreq,'revfilt',0,'channels',yChan); %HP (muscle)
    case 'BP',   tempFiltRerefEEG = pop_eegfiltnew(EEG,'locutoff',params.filtYfreq(1),'hicutoff',params.filtYfreq(2),'revfilt',0,'channels',yChan); %BP (line noise)
    case 'Notch',tempFiltRerefEEG = pop_eegfiltnew(EEG,'locutoff',params.filtYfreq(1),'hicutoff',params.filtYfreq(2),'revfilt',1,'channels',yChan); %notch (all artifacts)
    otherwise,   error('Detected invalid filtering option for Y.');
end
switch params.rerefY
    case 'no'
        Y_temp = tempFiltRerefEEG.data(yChan,:); %no reref
    case 'yes-temp'
        Y_temp = iCanClean_reref(tempFiltRerefEEG.data(yChan,:),params.rankY); %temp define a rereferenced version of Y (to help CCA later)
    case {'yes-fullrank','yes-loserank'} %old, for backwards compatibility
        Y_temp = iCanClean_reref(tempFiltRerefEEG.data(yChan,:),params.rerefY(5:end)); %temp define a rereferenced version of Y (to help CCA later)
    otherwise
        error('Detected invalid rereferencing option for Y.');    
end

%% define time windows (allows use of moving window)
t = EEG.times/1000; %seconds
tTotal = EEG.times(end)/1000;

if(params.useBoundary) %optionally use boundary events to define windows
    %define a vector of boundaries
    latencies = round([EEG.event.latency])/EEG.srate;
    types = {EEG.event.type};
    boundaries = [latencies(strcmp(types,'boundary')) tTotal];
    
    tStart = [];
    tEnd = [];
    for i=1:length(boundaries)-1
        tStart = horzcat(tStart,boundaries(i):params.stepSize:boundaries(i+1)); %#ok<AGROW> %array of times defining start of window (moving or static)
        if isempty(tStart), tStart = 0; end %helping to avoid a bug momentarily
        tEnd = horzcat(tEnd,(boundaries(i):params.stepSize:boundaries(i+1))+params.cleanWindow); %#ok<AGROW> %array of end times
        tEnd(end) = boundaries(i+1);
    end
else %ignore boundary events (default)
    tStart = 0:params.stepSize:tTotal; %array of times defining start of window (moving or static)
    if isempty(tStart), tStart = 0; end %helping to avoid a bug momentarily
    tEnd = tStart+params.cleanWindow; %array of end times
    tEnd(end) = tTotal;
end

%% start starts figure (fill in with data later)
if params.plotStatsOn
    myStatsFig = figure(); %#ok<NASGU>
    myStatsFig_R = gca;
end

%% added section: external calibration data checking
if params.useExternCalibData
    %get global ALLEEG structure (this may not be kosher)
    eeg_global;
    
    EEG_calib = ALLEEG(params.indExternCalibData);
    
    %check number of channels, channel names, channel order
    if(EEG_calib.nbchan ~= EEG.nbchan)
        error('External calibration data has the wrong number of channels.')
    end
    calibChanNames = {EEG_calib.chanlocs.labels};
    targetChanNames = {EEG.chanlocs.labels};
    if(~isempty(setdiff(calibChanNames,targetChanNames)))
        error('External calibration data has different channel names. ')
    end
    for i=1:length(targetChanNames)
        if ~strcmp(calibChanNames{i},targetChanNames{i})
            error('External calibration data has the wrong channel order.')
        end
    end
    disp('External calibration data passed all checks.')
end

%% optionally use "external" data within each window by first computing CCA on the entire dataset
if params.calcCCAonWholeData == true
    if params.useExternCalibData
        X_calib = EEG_calib.data(xChan,:)';
        Y_calib = EEG_calib.data(yChan,:)';
    else
        X_calib = EEG.data(xChan,:)';
        Y_calib = EEG.data(yChan,:)';
    end
    
    %RYAN! ROEHL! NEED TO ACCOUNT FOR TEMP FILTERING AND REREF IN FUTURE?? APPLY SAME
    %STEPS TO CALIB DATASET
    %X_calib_temp = ...
    
    %calculate noise sources using CCA
    if any(isnan(X_calib(:))) || any(isnan(Y_calib(:)))
        disp('I coud not do calib with entire dataset b/c you have NaN')
        A_calib = [];
        B_calib = [];
    else %it's safe to calc
        [A_calib,B_calib,~,U_calib,~,~] = canoncorr(X_calib, Y_calib );
        %          [A_calib,B_calib,R_calib,U_calib,V_calib,STATS_calib] = canoncorr(X_calib_temp, Y_calib_temp) );
        X_calib_MC = (X_calib - mean(X_calib,1));
        %             Y_calib_MC = (Y_calib - mean(Y_calib,1));
        A_inv_calib = pinv(A_calib(:,:)); %#ok<NASGU>
        fakeWinv_calib = mrdivide(X_calib_MC(:,:)',U_calib(:,:)');
        %NOTE 2022-05-06: code above does not currently allow user to specify U
        %versus V!
        clear X Y X_MC Y_MC U V
    end
else
    A_calib = [];
    B_calib = [];
end

%% calculate intermediates (added 2022-31-8)
if(~isinf(params.cleanWindow) && ~isinf(params.statsWindow))
    if(params.RTBool)
        extraTime_pre = params.statsWindow-params.cleanWindow;
        extraTime_post = 0;
    else
        extraTime_pre = (params.statsWindow-params.cleanWindow)/2;
        extraTime_post = extraTime_pre;
    end
else %case where window is infinite (added 2022-08-15)
    extraTime_pre = 0;
    extraTime_post = 0;
end

%% go thru each window, grab the data, calculate potential noise sources, and clean x and or y with method of choice
% msg_n = 0; %used for trying to display output in one line only
numWindows = length(tStart);
numNoiseCompsRemoved = zeros(1,numWindows);
for window_i = 1:numWindows
    
    window = find( t>=tStart(window_i) & t<=tEnd(window_i) );
    
    %% finding adjacent boundary events to define stats window
    if(params.useBoundary)
        bStartIndex = find(boundaries <= tStart(window_i),1,'last');
        bEndIndex = find(boundaries >= tEnd(window_i),1,'first');
        if(bEndIndex-bStartIndex > 1)
            warning('The selected boundary events are not consecutive.');
        end
        boundaryStart = boundaries(bStartIndex);
        boundaryEnd = boundaries(bEndIndex);
        if(boundaryEnd-boundaryStart < params.statsWindow)
            warning('Boundary-defined segment is shorter than stats window.');
        end
    else
        boundaryStart = 0;
        boundaryEnd = tTotal;
    end
    
    %% optionally mimic real-time windowing (pseduo real-time)
    if(params.RTBool && tEnd(window_i)-boundaryStart < params.statsWindow)
        disp('Not enough real-time data yet for this window, so we are skipping it.')
        continue;
    end
    
    %% broad window edge cases (formerly just the else case)
    BW_leftEdge = tStart(window_i)-extraTime_pre; %if start is < 0, we need to add post time
    BW_rightEdge = tEnd(window_i)+extraTime_post; %if end is > tTotal, we need to add pre time
    if(boundaryEnd-boundaryStart < params.statsWindow)
        broadWindow = find( t>=boundaryStart & t<=boundaryEnd );
    elseif(BW_leftEdge < boundaryStart)
        broadWindow = find( t>=boundaryStart & t<=tEnd(window_i)+extraTime_post+(boundaryStart-BW_leftEdge) );
    elseif(BW_rightEdge > boundaryEnd)
        addMorePre = params.cleanWindow - (tEnd(window_i)-tStart(window_i)); %last window in a segment is short
        broadWindow = find( t>=tStart(window_i)-extraTime_pre-(BW_rightEdge-boundaryEnd)-addMorePre & t<=boundaryEnd );
    else %otherwise, define normally
        broadWindow = find( t>=tStart(window_i)-extraTime_pre & t<=tEnd(window_i)+extraTime_post );
    end
    
    %% grab data in local window
    [~, locb1]  = ismember(window,broadWindow);
    
    X_localWindow = X_orig(:, broadWindow)';
    Y_localWindow = Y_orig(:, broadWindow)';
    X_localWindow_temp = X_temp(:,broadWindow)';
    Y_localWindow_temp = Y_temp(:,broadWindow)';
    
    %% calculate noise sources using CCA
    if any(isnan(X_localWindow_temp(:))) || any(isnan(Y_localWindow_temp(:))) %if missing data (NaN)
        disp('I coud not clean this section because you have NaN')
        nT = size(X_localWindow_temp,1); nX = size(X_localWindow_temp,2); nY = size(Y_localWindow_temp,2); nS = min(nX,nY);
        R = zeros(nS,1);
        U = zeros(nT,nS);
        V = U;
    else %we have data available so it's safe to run
        if isempty(A_calib) || isempty(B_calib) %if we don't use to use external data and instead want to calc cca on this smaller window
            [~,~,R,U,V] = iCanClean_calcSources(X_localWindow_temp,Y_localWindow_temp);
            %             [A,B,R,U,V,STATS] = iCanClean_calcSources( X_temp , X_temp ,params); %hey Roehl look here 2022-05-06
        else %we were given external A and B matrices to use
            A = A_calib;
            B = B_calib;
            
            %note from Ryan to Roehl: may want to do temp filt/reref here too 2022-05-06
            %X = X_temp
            X_MC = (X_localWindow_temp - mean(X_localWindow_temp,1));
            Y_MC = (Y_localWindow_temp - mean(Y_localWindow_temp,1));
            U = X_MC*A;
            V = Y_MC*B;
            [tempRHO] = corr(U,V);
            R = diag(tempRHO);
        end
    end
    
    %% Select bad noise sources based on R_squared correlation
    badComps = find(R.^2>params.rhoSqThres_source);
    numNoiseCompsRemoved(window_i) = length(badComps);
    
    %% update plots
    if params.plotStatsOn
        %plot stats
        hold(myStatsFig_R,'off');stem(myStatsFig_R, R.^2); hold(myStatsFig_R,'on'); stem(myStatsFig_R, badComps,R(badComps).^2,'r'); set(myStatsFig_R,'ylim',[0 1]); xlabel(myStatsFig_R,'Component #'); ylabel(myStatsFig_R,'R^2'); %legend('Kept','Rejected'); %hide legend plotting because adds too much time?
        drawnow;
        
        %% EXPERIMENTAL: variance accounted for (warning: slow to calculate)
        if params.plotVAF %2022-03-11 Ryan added (Roehl please add clickable button to turn VAF fig on/off)
            if window_i == 1
                vafFig = figure(); %#ok<NASGU>
                subplot(2,2,1);
                vafFig_U2X = gca;
                subplot(2,2,2);
                vafFig_U2Y = gca;
                subplot(2,2,3);
                vafFig_V2X = gca;
                subplot(2,2,4);
                vafFig_V2Y = gca;
            end

            vaf_U2X = cumsum(iCanClean_calcVarAccFor(X_localWindow_temp,U));
            vaf_V2X = cumsum(iCanClean_calcVarAccFor(X_localWindow_temp,V));
            vaf_U2Y = cumsum(iCanClean_calcVarAccFor(Y_localWindow_temp,U));
            vaf_V2Y = cumsum(iCanClean_calcVarAccFor(Y_localWindow_temp,V));
            
            
            hold(vafFig_U2X,'off'); stem(vafFig_U2X, vaf_U2X); hold(vafFig_U2X,'on'); set(vafFig_U2X,'ylim',[0 100]); set(vafFig_U2X,'xlim',[0 120]); title(vafFig_U2X,'U to X');
            hold(vafFig_U2Y,'off'); stem(vafFig_U2Y, vaf_U2Y); hold(vafFig_U2Y,'on'); set(vafFig_U2Y,'ylim',[0 100]); set(vafFig_U2Y,'xlim',[0 120]); title(vafFig_U2Y,'U to Y');
            hold(vafFig_V2X,'off'); stem(vafFig_V2X, vaf_V2X); hold(vafFig_V2X,'on'); set(vafFig_V2X,'ylim',[0 100]); set(vafFig_V2X,'xlim',[0 120]); title(vafFig_V2X,'V to X');
            hold(vafFig_V2Y,'off'); stem(vafFig_V2Y, vaf_V2Y); hold(vafFig_V2Y,'on'); set(vafFig_V2Y,'ylim',[0 100]); set(vafFig_V2Y,'xlim',[0 120]); title(vafFig_V2Y,'V to Y');
            drawnow;
              
        end
    end
    %% define which noise sources to use to clean each of the original channels
    %i.e. do you want to clean using mixtures of X or mixtures of Y?
    switch params.cleanXwith
        case {'X','U'},     xNoiseSources = U(:,badComps);
        case {'Y','V'},     xNoiseSources = V(:,badComps);
        case {'XY','UV'},   xNoiseSources = ( U(:,badComps) + V(:,badComps) )./2;
        otherwise, error('you did not define params.cleanXwith properly');
    end
    
    switch params.cleanYwith
        case {'X','U'},     yNoiseSources = U(:,badComps);
        case {'Y','V','no'},yNoiseSources = V(:,badComps);
        case {'XY','UV'},   yNoiseSources = ( U(:,badComps)+V(:,badComps) )./2;
        otherwise, error('you did not define params.cleanYwith properly');
    end
    
    %% clean channels
    if isempty(A_calib) || isempty(B_calib)%ryan new 11/30/2020 quick test to see if we want to clean another way
        tempClean = iCanClean_cleanChansWithNoiseSources(X_localWindow,xNoiseSources)';
    else %we want to use calib data to clean
        disp('THIS IS EXPERIMENTAL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        X_MC = (X_localWindow - mean(X_localWindow,1));
        
        %         %option 1: reject with psuedo inverse
        %         X_est = X_MC*A_calib(:,badComps)*A_inv_calib(badComps,:);
        
        %option 2: reject with matrix division (least square?)
        %            fakeWinv_calib = mrdivide(X_MC(:,:)',U(:,:)'); disp('doing local scaling'); %calc locally if line executed. comment out if you want to use what was calc from larger dataset
        X_est = (fakeWinv_calib(:,badComps)*(X_MC*A_calib(:,badComps))')';
        %         X_est = X_MC*A_calib(:,badComps)*fakeWinv_calib(:,badComps)'; %equiv expression to line above? check and see
        
        tempClean = (X_localWindow-X_est)';
    end
    cleanEEG.data(xChan,window) = tempClean(:,locb1); clear tempClean;
    
    %clean y if desired
    if ~strcmpi(params.cleanYwith,'no') && ~isempty(yChan)
        if ~isempty(intersect(xChan,yChan))
            error('I am not yet programmed to prevent accidentally double cleaning in the case that xChan and yChan share one or more channels. Sorry but we can program this in the future');
        end
        tempClean = iCanClean_cleanChansWithNoiseSources(Y_localWindow,yNoiseSources)';%2022-07-15 Chang caught typo (X_localWindow versus Y_localWindow)
        cleanEEG.data(yChan,window) = tempClean(:,locb1); clear tempClean;
    end
    %NOT SET UP FOR CLEANING Y CHANNELS USING CALIB DATA YET!
            %         Y_MC = (Y_localWindow - mean(Y_localWindow,1));
    
    %% update user on progress (assuming moving window)
    %     fprintf(repmat('\b',1,msg_n+1));
    if params.giveCleaningUpdates
        msg = ['iCanClean cleaned [',num2str([tStart(window_i) tEnd(window_i)]),']s. Removed ',num2str(length(badComps)),' bad sources in this window'];
        disp(msg);
        
        %     msg_n = length(msg);
    end
end

%% finish keeping track of time
disp('Finished iCanClean cleaning');
toc

%% update EEG structure on the way out
cleanEEG.etc.iCanClean.numNoiseCompsRemovedPerWindow = numNoiseCompsRemoved;
cleanEEG.etc.iCanClean.numNoiseCompsRemovedOnAvg = mean(numNoiseCompsRemoved);
cleanEEG.etc.iCanClean.params = params;
disp(['Removed ',num2str(cleanEEG.etc.iCanClean.numNoiseCompsRemovedOnAvg),' noise components on average (across all time windows cleaned)']);

%% visualize final results
if visualizeFinalResults
    disp('Visualizing results. This may take a while.');
    %     iCanClean_visualizeResults(cleanEEG,EEG,xChan(1:4:end),yChan(1:8:end));
    iCanClean_visualizeResults(cleanEEG,EEG,xChan(1:4:end),[]); 
end
end
