% pop_findBrainComps() - perform EEG artifact removal using iCanClean algorithm.
%                If less than three arguments are given, a window pops up
%                to ask for the value of the additional parameters.   
%
% Usage:
%   >>  OUTEEG = pop_iCanClean(EEG); % pop-up window asking for parameters
% 
%
% Inputs:
%   EEG     - input EEG dataset
%
%    
% Outputs:
%   OUTEEG  - output dataset
%
% See also:
%   iCanClean, EEGLAB 

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
% To cite iCanClean, please reference the following paper:
%
% Downey, R.J.; Ferris, D.P. iCanClean Removes Motion, Muscle, Eye, 
% and Line-Noise Artifacts from Phantom EEG. Sensors 2023, 23, 8214. 
% https://doi.org/10.3390/s23198214
%
function [ EEG ] = pop_findBrainComps( EEG )
%POP_FINDBRAINCOMPS GUI menu for finding brain components
%   Defines GUI parameters and parses inputs

%% define geometry

    % Rows sum to 60 (unless blank)
    geomhoriz = {[20 20 20] ...
                    1 ...
                    [20 15 5 8 12] ...
                    1 ...
                    [20 6 6 4 24] ...
                    [20 40] ...
                    1 ...
                    [4 16 6 34] ...
                    1 ...
                    [4 16 6 10 24] ...
                    1 ...
                    [20 10 30]
                    };
    
    geomvert = [1 0.5 1 0.2 1 1 0.5 1 0.2 1 1 1];

%% pre-fill conditions
    chanTypes = unique({EEG.chanlocs.type});
    if ismember('EEG',chanTypes)
        chanTypePrefill = 'EEG';
    elseif length(chanTypes)==1
        chanTypePrefill = chanTypes(1);
    else
        chanTypePrefill = '';
    end
    
%% define callbacks (if any)
    
    cb_chan = ['tmpchanlocs = EEG(1).chanlocs;' ...
            '[tmp tmpval] = pop_chansel({tmpchanlocs.labels},''withindex'',''on'');' ...
            'set(findobj(gcbf, ''tag'', ''chans''), ''string'',tmpval);' ...
            'clear tmp tmpchanlocs tmpval'];
    cb_chanType = strrep(cb_chan,'labels},''withindex'',''on','type},''field'',''type');
    cb_chanSpec = ['set(findobj(gcbf, ''tag'', ''chanButton''), ''callback'',' ... 
            'fastif(get(findobj(gcf,''tag'',''chanSpec''),''value'')==2,''' ...
            strrep(cb_chan,'''','''''') ''',''' strrep(cb_chanType,'''','''''') '''));'];
    cb_filtType = ['freqs = {'''',''5'',''25'',''55 65'',''5 25''};',...
                    'set(findobj(gcbf, ''tag'', ''filterFreq''),''string'',freqs{get(findobj(gcbf,''tag'',''filter''),''value'')});'];
    cb_plotCheck = ['set(findobj(gcbf,''tag'',''maxCompsCoeff''),''enable'',fastif(get(findobj(gcf,''tag'',''plotCheck''),''value''),''on'',''off''));',...
                   'set(findobj(gcbf,''tag'',''maxComps''),''enable'',fastif(get(findobj(gcf,''tag'',''plotCheck''),''value''),''on'',''off''));'];
    cb_usePrevInputs = ['homefolder = deblank(evalc(''!echo %USERPROFILE%''));',...
                     'if ~exist(fullfile(homefolder,''findBrainComps_options.mat'')), warndlg(''Previous run not found.'',''Cannot load previous inputs''); return; end;',...
                     'load(fullfile(homefolder,''findBrainComps_options.mat''));',...
                     'fields = fieldnames(structout);',...
                     'try,',...
                     'for ii=1:length(fields), set(findobj(gcbf,''tag'',fields{ii}),fastif(isnumeric(getfield(structout,fields{ii})), ''value'', ''string''),getfield(structout,fields{ii})); end;',...
                     'catch, warndlg(''The options file is from a previous version of iCanClean.'',''Cannot load previous inputs''); end;'];
    cb_resetDefaults = ['chanTypes = unique({EEG.chanlocs.type});',...
                        'if ismember(''EEG'',chanTypes), chanTypePrefill = ''EEG'';',...
                        'elseif length(chanTypes==1), chanTypePrefill = chanTypes(1);',...
                        'else, chanTypePrefill = ''''; end;',...
                        'defaultStruct = struct(''chans'',chanTypePrefill,''chanSpec'',1,',...
                        ' ''filter'',1,''filterFreq'','''',''maxCompsCoeff'',''0.5'',',...
                        ' ''maxComps'',''35'',''plotCheck'',0);',...
                        'fields = fieldnames(defaultStruct);',...
                        'for ii=1:length(fields), set(findobj(gcbf,''tag'',fields{ii}),fastif(isnumeric(getfield(defaultStruct,fields{ii})), ''value'', ''string''),getfield(defaultStruct,fields{ii})); end;',...
                        cb_plotCheck];

%% create uilist
    
    textFontSize = 12;
    
    uilist = {...
        {'style' 'pushbutton' 'string' 'Use previous inputs' 'callback' cb_usePrevInputs},...
        {'style' 'pushbutton' 'string' 'Reset to defaults' 'callback' cb_resetDefaults},...
        {},...
        ...
        {},...
        ...
        {'style' 'text' 'string' 'Define EEG channels' 'fontsize' textFontSize},...
        {'style' 'edit' 'string' chanTypePrefill 'tag' 'chans'},...
        {'style' 'pushbutton' 'string' '...' 'tag' 'chanButton' 'enable' ...
            fastif(isempty(EEG(1).chanlocs), 'off', 'on') 'callback' cb_chanType},...
        {'style' 'text' 'string' 'Specify by:'},...
        {'style' 'popupmenu' 'string' 'chan type|chan name' 'tag' 'chanSpec' 'callback' cb_chanSpec},...
        ...
        {},...
        ...
        {'style' 'text' 'string' 'Temp filter to remove brain activity' 'fontsize' textFontSize},...
        {'style' 'popupmenu' 'string' 'no|LP|HP|BP|Notch' 'tag' 'filter' 'callback' cb_filtType 'value' 5},... %roehl edited default value 2022-9-6
        {'style' 'edit' 'string' '5 25' 'tag' 'filterFreq'},... %roehl edited default value 2022-9-6
        {'style' 'text' 'string' 'Hz' 'fontsize' textFontSize},...
        {},...
        ...
        {'style' 'checkbox' 'string' '(optional) Plot Brain Components' 'tag' 'plotCheck' 'callback' cb_plotCheck 'value' 1},... %roehl edited default value 2022-9-6
        {},...
        ...
        {},...
        ...
        {},...
        {'style' 'text' 'string' 'Plot comps with R^2<' 'fontsize' textFontSize},...
        {'style' 'edit' 'string' '0.5' 'tag' 'maxCompsCoeff' 'enable' 'on'},...
        {},...
        ...
        {},...
        ...
        {},...
        {'style' 'text' 'string' 'Plot at most' 'fontsize' textFontSize},...
        {'style' 'edit' 'string' '35' 'tag' 'maxComps' 'enable' 'on'},...
        {'style' 'text' 'string' 'component(s)' 'fontsize' textFontSize},...
        {},...
        ...
        {},...
        ...
        {'style' 'text' 'string' 'Method to compute inverse weights:'},...
        {'style' 'popupmenu' 'string' 'pinv|regression' 'tag' 'pinvBool'},...
        {},...
        };
    
%% instantiate GUI

    [~, ~, ~, structout] = inputgui('geometry',geomhoriz,'geomvert',geomvert,'uilist',uilist,...
                           'helpcom','pophelp(''pop_findBrainComps'')',...
                           'title','Find brain components -- pop_findBrainComps()');
    
    if(isempty(structout))
        return; %return early if GUI is cancelled
    end

%% parse outputs

% save as previous run
    homefolder = deblank(evalc('!echo %USERPROFILE%'));
    save(fullfile(homefolder,'findBrainComps_options.mat'),'structout');

% get input params from user
    chanTypeOfInterest = structout.chans;
    EEG_chans = find(strcmpi({EEG.chanlocs.type},chanTypeOfInterest));

%remove any existing ICA results
EEG.icaact = [];
EEG.icaweights = [];
EEG.icawinv = [];
EEG.icasphere = [];
% try EEG = rmfield(EEG,'dipfit'); catch; end
try EEG.dipfit=[]; catch; end
try EEG.etc.ic_classification = rmfield(EEG.etc.ic_classification,'ICLabel'); catch; end

% duplicate and filt chans
    filtTypes = {'no' 'LP' 'HP' 'BP' 'Notch'};
    filterType = cell2mat(filtTypes(structout.filter));
    if strcmp(filterType,'no')
        %Grab channels of interest 
        tempEEG = pop_select(EEG, 'channel',EEG_chans);
        %update chan names
        for ch_i = 1:tempEEG.nbchan
            tempEEG.chanlocs(ch_i).labels = ['Filt-', tempEEG.chanlocs(ch_i).labels];
            tempEEG.chanlocs(ch_i).type = ['Filt-', tempEEG.chanlocs(ch_i).type];
        end
        %append  filtered data
        EEG.data = [EEG.data; tempEEG.data];
        try
        EEG.chanlocs = [EEG.chanlocs, tempEEG.chanlocs];
        catch
        EEG.chanlocs = [EEG.chanlocs; tempEEG.chanlocs];
        end
    else
        filterFreqNum = sort(str2num(structout.filterFreq));
        switch filterType
            case 'LP'
                locutoff = [];
                hicutoff = filterFreqNum;
            case 'HP'
                locutoff = filterFreqNum;
                hicutoff = [];
            case 'BP'
                locutoff = filterFreqNum(1);
                hicutoff = filterFreqNum(2);
            case 'Notch'
                locutoff = filterFreqNum(2);
                hicutoff = filterFreqNum(1);
        end
        EEG = duplicateAndFilterChans(EEG,EEG_chans,locutoff,hicutoff);
    end
    
    Noise_chans = find(strcmpi(['Filt-', chanTypeOfInterest],{EEG.chanlocs.type}));

  
    if structout.pinvBool==1 %if you want to do pseudo inverse
%         EEG.icawinv = pinv(EEG.icaweights*EEG.icasphere);
            projectionMethod = 0;
    else
%         EEG.icawinv = mrdivide(EEG.data,EEG.icaact);
            projectionMethod = 1;
    end

% make equiv ICA
    EEG = makeEquivICA_iCanClean(EEG,EEG_chans,Noise_chans,projectionMethod); %main thing that happens


    
% overwrite or new dataset
%     [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 ); %new dataset
%     [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET ); %overwrite current set
%     eeglab redraw;

% optionally view components
    if structout.plotCheck==1
		numCompsAvail = size(EEG.icaact,1);
        maxCompsToPlot = str2num(structout.maxComps);
        numChans = length(EEG_chans);
        
        if numCompsAvail > maxCompsToPlot
            numToPlot = maxCompsToPlot;
        else %have less components available than desired
            numToPlot = numCompsAvail;
        end
        compIndToPlot = (numCompsAvail-numToPlot+1):numCompsAvail;
        
        % NOT IMPLEMENTED???? Auto exclude components (from compIndToPlot) with r^2 > maxCompsCoeff 
        try
            pop_viewprops( EEG, 0, compIndToPlot, {'freqrange', [2 80]}, {}, 1, '' ) %optional viewing properties
        catch
            warning('Please make sure you have installed the ViewProps plugin for eeglab');
        end
    end
    
% Output data
%     OUTEEG = EEG;
    
end

