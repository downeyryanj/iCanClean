% iCanClean_3DEpochedData() - perform EEG artifact removal using iCanClean 
%                             algorithm. Modified for an epoched dataset.
%
% Usage:
%   >>  cleanEEG = iCanClean_3DEpochedData( EEG, xChan, yChan, visualizeFinalResults, params)
%
% Inputs:
%   EEG     - input EEG dataset that is 3-Dimensional (epoched data)
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

function cleanEEG = iCanClean_3DEpochedData( EEG, xChan, yChan, visualizeFinalResults, params)

if nargin < 1
    help iCanClean;
    return;
end

%% initialize EEG structure that will eventually be our clean data and keep track of time
cleanEEG = EEG;
tic;

%% CLEANUP check/set parameters
params = iCanClean_setParams(params)

% if size(EEG.data,3)>1 %if epoched data
%     %set window length to epoch length
%     %
% end



%start starts figure (fill in with data later)
if params.plotStatsOn
    myStatsFig = figure();
    myStatsFig_R = gca;
end

%% optionally use "external" data within each window by first computing CCA on the entire dataset
if params.calcCCAonWholeData == true
    X = EEG.data(xChan,:)';
    Y = EEG.data(yChan,:)';
    
    
    %calculate noise sources using CCA
    if any(isnan(X(:))) | any(isnan(Y(:)))
        disp('I coud not do calib with entire dataset b/c you have NaN')
        A_calib = [];
        B_calib = [];
    else %it's safe to calc
        [A_calib,B_calib,R_calib,U_calib,V_calib,STATS_calib] = canoncorr(X, Y );
                    X_MC = (X - mean(X,1));
%             Y_MC = (Y - mean(Y,1)); 
            U = X_MC*A_calib;
%             V = Y_MC*B_calib;
            A_inv_calib = pinv(A_calib(:,:));
            fakeWinv_calib = mrdivide(X_MC(:,:)',U(:,:)');
            
        clear X Y X_MC Y_MC U V
    
    end
    
else
    A_calib = [];
    B_calib = [];
end



%% go thru each window, grab the data, calculate potential noise sources, and clean x and or y with method of choice
% msg_n = 0; %used for trying to display output in one line only
numWindows = size(EEG.data,3); %length(tStart);
numNoiseCompsRemoved = zeros(1,numWindows);
for window_i = 1:numWindows  
    locb1 = 1:size(EEG.data,2);
    
    %% grab data
    X = EEG.data(xChan,:,window_i)';
    Y = EEG.data(yChan,:,window_i)';
	
    
    %calculate noise sources using CCA
    if any(isnan(X(:))) | any(isnan(Y(:))) %if missing data (NaN)
        disp('I coud not clean this section because you have NaN')
        nT = size(X,1); nX = size(X,2); nY = size(Y,2); nS = min(nX,nY);
        R = zeros(nS,1);
        U = zeros(nT,nS);
        V = U;
    else %we have data available so it's safe to run
        if isempty(A_calib) || isempty(B_calib) %if we don't use to use external data and instead want to calc cca on this smaller window
            [A,B,R,U,V,STATS] = iCanClean_calcSources(X,Y,params);
        else %we were given external A and B matrices to use
            A = A_calib; 
            B = B_calib;
             
            X_MC = (X - mean(X,1));
            Y_MC = (Y - mean(Y,1)); 
            U = X_MC*A;
            V = Y_MC*B;
            [tempRHO,tempPVAL] = corr(U,V);
            R = diag(tempRHO);
            STATS.p = diag(tempPVAL);
        end
    end
    


    % define bad noise sources at source level (cancor stats)
    badComps = find(R.^2>params.rhoSqThres_source & STATS.p < params.pThres_source);
%     badComps = find(R.^2<params.rhoSqThres_source); disp('Ryan flipped Rsq criterion just as a quick test!'); %2022-03-31 ryan testing out something "wrong" on purpose
    numNoiseCompsRemoved(window_i) = length(badComps);


    if params.plotStatsOn
        %plot stats
        hold(myStatsFig_R,'off');stem(myStatsFig_R, R.^2); hold(myStatsFig_R,'on'); stem(myStatsFig_R, badComps,R(badComps).^2,'r'); set(myStatsFig_R,'ylim',[0 1]); xlabel(myStatsFig_R,'Component #'); ylabel(myStatsFig_R,'R^2'); %legend('Kept','Rejected'); %hide legend plotting because adds too much time?
        drawnow;
		
        %% added section: variance accounted for
        if params.plotVAF %2022-03-11 Ryan added (Roehl please add clickable button to turn VAF fig on/off)
            if window_i == 1
                vafFig = figure();
                subplot(2,2,1);
                vafFig_U2X = gca;
                subplot(2,2,2);
                vafFig_U2Y = gca;
                subplot(2,2,3);
                vafFig_V2X = gca;
                subplot(2,2,4);
                vafFig_V2Y = gca;
            end
            
             vaf_U2X = iCanClean_calcVarAccFor(X,U);
             vaf_V2X = iCanClean_calcVarAccFor(X,V);
             vaf_U2Y = iCanClean_calcVarAccFor(Y,U);
             vaf_V2Y = iCanClean_calcVarAccFor(Y,V);
            
            
            hold(vafFig_U2X,'off'); stem(vafFig_U2X, vaf_U2X); hold(vafFig_U2X,'on'); set(vafFig_U2X,'ylim',[0 100]); set(vafFig_U2X,'xlim',[0 120]); title(vafFig_U2X,'U to X');
            hold(vafFig_U2Y,'off'); stem(vafFig_U2Y, vaf_U2Y); hold(vafFig_U2Y,'on'); set(vafFig_U2Y,'ylim',[0 100]); set(vafFig_U2Y,'xlim',[0 120]); title(vafFig_U2Y,'U to Y');
            hold(vafFig_V2X,'off'); stem(vafFig_V2X, vaf_V2X); hold(vafFig_V2X,'on'); set(vafFig_V2X,'ylim',[0 100]); set(vafFig_V2X,'xlim',[0 120]); title(vafFig_V2X,'V to X');
            hold(vafFig_V2Y,'off'); stem(vafFig_V2Y, vaf_V2Y); hold(vafFig_V2Y,'on'); set(vafFig_V2Y,'ylim',[0 100]); set(vafFig_V2Y,'xlim',[0 120]); title(vafFig_V2Y,'V to Y');
            drawnow;
                    

        end
    end
    %% define which noise sources to use to clean each of the original channels
    %i.e. do you want to clean using mixtures of X or mixtures of Y?
    if strcmpi(params.cleanXwith,'X') | strcmpi(params.cleanXwith,'U')
        xNoiseSources = U(:,badComps);
    elseif strcmpi(params.cleanXwith,'Y') | strcmpi(params.cleanXwith,'V')
        xNoiseSources = V(:,badComps);
    elseif strcmpi(params.cleanXwith,'XY') | strcmpi(params.cleanXwith,'UV')
        xNoiseSources = ( U(:,badComps) + V(:,badComps) )./2;
    else
        error('you did not define params.cleanXwith properly');
    end
    
    
    if strcmpi(params.cleanYwith,'X') | strcmpi(params.cleanYwith,'U')
        yNoiseSources = U(:,badComps);
    elseif strcmpi(params.cleanYwith,'Y') | strcmpi(params.cleanYwith,'V')
        yNoiseSources = V(:,badComps);
    elseif strcmpi(params.cleanYwith,'XY') | strcmpi(params.cleanYwith,'UV') %typo prior to 2021-09-23 (cleanXwith vs cleanYwith) 
        yNoiseSources = ( U(:,badComps)+V(:,badComps) )./2;
    else
        error('you did not define params.cleanYwith properly');
    end
    
    %% clean channels
    %clean x
    if isempty(A_calib) || isempty(B_calib)%ryan new 11/30/2020 quick test to see if we want to clean another way
        tempClean = iCanClean_cleanChansWithNoiseSources(X,xNoiseSources,params.noiseRemovalMethod,params)';
    else %we want to use calib data to clean
        disp('THIS IS EXPERIMENTAL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        X_MC = (X - mean(X,1));
%         Y_MC = (Y - mean(Y,1)); 
        U = X_MC*A_calib;

%         %option 1: reject with psuedo inverse
%         X_est = X_MC*A_calib(:,badComps)*A_inv_calib(badComps,:);
        
        %option 2: reject with matrix division (least square?)
%            fakeWinv_calib = mrdivide(X_MC(:,:)',U(:,:)'); disp('doing local scaling'); %calc locally if line executed. comment out if you want to use what was calc from larger dataset
        X_est = (fakeWinv_calib(:,badComps)*(X_MC*A_calib(:,badComps))')';
%         X_est = X_MC*A_calib(:,badComps)*fakeWinv_calib(:,badComps)'; %equiv expression to line above? check and see
     
        tempClean = (X-X_est)';
    end
    cleanEEG.data(xChan,:,window_i) = tempClean(:,locb1); clear tempClean;
    
    
    %clean y if desired
    if params.cleanYBool && ~isempty(yChan)
        tempClean = iCanClean_cleanChansWithNoiseSources(Y,yNoiseSources,params.noiseRemovalMethod,params)';
        cleanEEG.data(yChan,:,window_i) = tempClean(:,locb1); clear tempClean;
    end
    
    %update user on progress (assuming moving window)
    %     fprintf(repmat('\b',1,msg_n+1));
    if params.giveCleaningUpdates
        msg = ['iCanClean cleaned epoch[',num2str(window_i),']. Removed ',num2str(length(badComps)),' bad sources in this window'];
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
