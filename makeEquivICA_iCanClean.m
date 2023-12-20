%[EEG_out] = makeEquivICA_iCanClean(EEG_in,xChan,yChan,projectionMethod)
% 
% makeEquivICA_iCanClean runs iCanClean but shows you the components as if
% it were ICA rather than try to clean your channel-level data
%
% EEG_in = your EEG dataset
% xChan = eeg chan
% yChan = ref noise (or pseudo ref noise) chan
% projectionMethod = 0 for pinv, 1 for regression

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

function [EEG_out] = makeEquivICA_iCanClean(EEG_in,xChan,yChan,projectionMethod)

nChanX = length(xChan);

X = double(EEG_in.data(xChan,:))';
Y = double(EEG_in.data(yChan,:))';


[A B R U V] = canoncorr(X,Y);

EEG_out = pop_select( EEG_in,'channel',xChan); %keep only EEG chans (xChan) to make our ICA equiv dataset 



switch projectionMethod
    case 0
        EEG_out.icasphere = eye(nChanX);%row dim matters
        EEG_out.icaweights = A'; %column dim matters
        EEG_out.icawinv = pinv(EEG_out.icaweights*EEG_out.icasphere);   %option 1
        pop_editoptions('option_computeica',1);
        EEG_out = eeg_checkset(EEG_out,'ica');
    case 1
%     disp('overwritting inverse weights EXPERIMENTAL!!!!');
    secretSauce = mrdivide(EEG_out.data,EEG_out.icaact);            %option 2
    EEG_out.icawinv = secretSauce;%alternative topograph calc option
    otherwise
        error('you did not provide a projection method (pinv or regression)')
end


end

