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

function X = iCanClean_reref(X, rankString) % function EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankBool)
%Perform average re-reference cort to cort, noise to noise, and emg to emg
%   Usage: X = iCanClean_reref(X, 'fullrank')
%   X is expected to be row vectors. For example, EEG.data(desiredChannels,:)
%

fullRankBool = strcmpi(rankString, 'fullrank');
loseRankBool = strcmpi(rankString, 'loserank');
validStringBool = fullRankBool || loseRankBool;

if(~validStringBool)
    error('Error: invalid reref string, choose either ''fullrank'' or ''loserank''.')
end

if fullRankBool
    disp('Rereferencing to average (with zero loss of rank)');
else
    disp('Rereferencing to average (lose 1 rank)');
end

nchansin = min(size(X));
if fullRankBool
    refmatrix = eye(nchansin)-ones(nchansin)*1/(nchansin+1);%full rank method
else
    refmatrix = eye(nchansin)-ones(nchansin)*1/nchansin; %normal avg ref (rank losing method)
end
X = refmatrix*X;


end