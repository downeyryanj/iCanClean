% iCanClean_3DEpochedDataParallel() - perform EEG artifact removal using
%                                     iCanClean algorithm. Modified for an 
%                                     epoched dataset. Includes
%                                     parallelization.
%
% Usage:
%   >>  cleanEEG = iCanClean3DEpochedDataParallel( EEG, xChan, yChan, visualizeFinalResults, params)
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

function cleanEEG = iCanClean_3DEpochedDataParallel( EEG, xChan, yChan, visualizeFinalResults, params)

if nargin < 1
    help iCanClean;
    return;
end

%% initialize EEG structure that will eventually be our clean data and keep track of time
cleanEEG = EEG;
tic;

%% check/set parameters
params = iCanClean_setParams(params);


%% optionally use "external" data within each window by first computing CCA on the entire dataset
if params.calcCCAonWholeData == true
    X = EEG.data(xChan,:)';
    Y = EEG.data(yChan,:)';

    
    %calculate noise sources using CCA
    if any(isnan(X(:))) || any(isnan(Y(:)))
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

%parfor test
tempData = cleanEEG.data;   %dumb stuff we have to initialize for parfor
tempDataxChan = tempData(xChan,:,:);
tempDatayChan = tempData(yChan,:,:);
parfor window_i = 1:numWindows

    locb1 = 1:size(tempData,2);
    
    %% grab data
%     X = EEG.data(xChan,:,window_i)';
%     Y = EEG.data(yChan,:,window_i)';
    X = tempDataxChan(:,:,window_i)';
    Y = tempDatayChan(:,:,window_i)';
    
    
    %calculate noise sources using CCA
    if any(isnan(X(:))) || any(isnan(Y(:))) %if missing data (NaN)
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
    badComps = find(R.^2<params.rhoSqThres_source); %disp('Ryan flipped Rsq criterion just as a quick test!'); %2022-03-31 ryan testing out something "wrong" on purpose
    numNoiseCompsRemoved(window_i) = length(badComps);


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
    %clean x
        tempClean = iCanClean_cleanChansWithNoiseSources(X,xNoiseSources)';

    
    %parfor test
    tempDataxChan(:,:,window_i) = tempClean(:,locb1); %clear tempClean;

    
    
    %clean y if desired
    if ~strcmpi(params.cleanYwith,'no') && ~isempty(yChan)
        tempClean = iCanClean_cleanChansWithNoiseSources(Y,yNoiseSources)';
        %parfor test
        tempDatayChan(:,:,window_i) = tempClean(:,locb1); %clear tempClean;

    end
    
    %update user on progress (assuming moving window)
    %     fprintf(repmat('\b',1,msg_n+1));
    if params.giveCleaningUpdates
        msg = ['iCanClean cleaned epoch[',num2str(window_i),']. Removed ',num2str(length(badComps)),' bad sources in this window'];
        disp(msg);
        
        %     msg_n = length(msg);
    end
end
%parfor test
tempData(xChan,:,:) = tempDataxChan;
tempData(yChan,:,:) = tempDatayChan;
cleanEEG.data = tempData;

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
