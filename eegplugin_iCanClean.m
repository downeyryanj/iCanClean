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

% plugin function for iCanClean tool
function vers = eegplugin_iCanClean(fig, try_strings, catch_strings)
vers = '1.0.2'; %2024-06-30

% create a callback function to pop_iCanClean
cmd = '[EEG,LASTCOM] = pop_iCanClean(EEG, 1);'; %typically, uimenu calls pop function
cb_findBrainComps = '[EEG] = pop_findBrainComps(EEG);';


% check for a non-empty EEG dataset, add to history and ALLEEG variable
finalcmd = [ try_strings.no_check cmd ];
finalcmd = [ finalcmd catch_strings.new_and_hist ];
finalcb_findBrainComps = [ try_strings.no_check cb_findBrainComps ];
finalcb_findBrainComps = [ finalcb_findBrainComps catch_strings.new_and_hist ];

% create a menu item under Tools
toolsmenu = findobj(fig, 'tag', 'tools');
iCCmenu = uimenu(toolsmenu,'label','iCanClean');
uimenu(iCCmenu,'label','Clean dataset','callback',finalcmd);
uimenu(iCCmenu,'label','Find brain components (experimental)','callback',finalcb_findBrainComps);


