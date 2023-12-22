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

function [EEG] = duplicateAndFilterChans(EEG,chanInd,locutoff, hicutoff)
%Duplicate and filter EEG channels
%   Literally takes your EEG channels, duplicates them, then filters them
%
%	EEG is your EEG data set
%	chanInd are what channels you want to duplicate (subset of all channels)
%	locutoff is the low cut off frequency (hz) for a band pass filter
%	hicutoff is the high cut off frequency for a band pass filter
% 	Note: if you make hicutoff > locutoff, you have a band pass filter
%			and if you make hicutoff < locutoff, you have a notch (band stop) filter
%			and you you input just hicutoff or just lowcutoff, it should be just a LP or HP filter.

%Grab channels of interest 
tempEEG = pop_select(EEG, 'channel',chanInd);

%update chan names
for ch_i = 1:tempEEG.nbchan
    tempEEG.chanlocs(ch_i).labels = ['Filt-', tempEEG.chanlocs(ch_i).labels];
    tempEEG.chanlocs(ch_i).type = ['Filt-', tempEEG.chanlocs(ch_i).type];
end

%Filter data
if hicutoff < locutoff %should not be true if one of locutoff or hicutoff is empty
    %disp('hi ryan');
    tempEEG = pop_eegfiltnew(tempEEG, 'locutoff',locutoff,'hicutoff',hicutoff,'revfilt',true); %notch filter
else
tempEEG = pop_eegfiltnew(tempEEG, 'locutoff',locutoff,'hicutoff',hicutoff); %BP
end

%append  filtered data
EEG.data = [EEG.data; tempEEG.data];

%fill potentially empty fields (dumb bug)
% find all possible fields
f = fieldnames(EEG.chanlocs);
for field_i = 1:length(f)
    FIELD = f{field_i};
    emptyMat = [];
    if ~isfield(tempEEG.chanlocs,FIELD)
        [tempEEG.chanlocs.(FIELD)] = deal(emptyMat);
    end
end

%append chanlocs
try
EEG.chanlocs = [EEG.chanlocs, tempEEG.chanlocs];
catch
EEG.chanlocs = [EEG.chanlocs; tempEEG.chanlocs];
end


end

