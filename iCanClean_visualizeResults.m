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

function [outputArg1,outputArg2] = iCanClean_visualizeResults(cleanEEG,EEG,xChanPlot,yChanPlot)
%Visualize spectra and time series results
%   Detailed explanation goes here
tempEEG = cleanEEG;



%% visualize spectra results
vis_spec = false;
if vis_spec
    % timeWindow = find(EEG.times/1000 >= 0 & EEG.times/1000 <= 200); %arbitrarily calculating PSD no smaller window so we don't have to wait forever
    timeWindow = find(EEG.times/1000 >= 0 & EEG.times/1000 <= inf);
    [spectra_in,FREQ,~,~,~] = spectopo(EEG.data(xChanPlot,timeWindow), length(timeWindow), EEG.srate, 'freqfac',1, 'plot','off');
    [spectra_out,FREQ,~,~,~] = spectopo(tempEEG.data(xChanPlot,timeWindow), length(timeWindow), EEG.srate, 'freqfac',1, 'plot','off');
    tempDiffMat =tempEEG.data(xChanPlot,timeWindow)-EEG.data(xChanPlot,timeWindow);
    if any(tempDiffMat(:)~=0) %make sure there is some non-zero difference between the datasets somewhere in any channel or time point
        spectra_diff =  spectopo( tempDiffMat, length(timeWindow), EEG.srate, 'freqfac',1, 'plot','off');
    else
        spectra_diff = nan(size(spectra_in)); disp('No difference between the clean and dirty data sets!')
    end
    % rnk = rank(double(tempEEG.data)); % Calculate rank of output data (component removal reduces rank)
    % fprintf(['CCA EEG rank out: ' int2str(rnk) '\n'])
    
    figure; % Plot spectra of mean channel input and output
    set(gcf,'Position',[500 200 500 400]);
    set(gcf,'Color','w');
    title('Power spectra');
    hold on;
    plot(FREQ, median(spectra_in,1), 'r','LineWidth', .5); plot(FREQ, prctile(spectra_in,25), 'r--','LineWidth', .5); plot(FREQ, prctile(spectra_in,75), 'r--','LineWidth', .5);
    set(gca,'FontSize',12);
    plot(FREQ, median(spectra_out,1), 'b','LineWidth', .5); plot(FREQ, prctile(spectra_out,25), 'b--','LineWidth', .5); plot(FREQ, prctile(spectra_out,75), 'b--','LineWidth', .5);
    plot(FREQ, median(spectra_diff,1), 'k','LineWidth',.5);
    legend('EEG in Q2', 'Q1','Q3', 'EEG out Q2','Q1','Q3','Diff Q2','Location', 'North', 'Orientation', 'horizontal');
    legend boxoff;
    xlabel('Frequency (Hz)');
    ylabel('Power 10*log_{10} (\muV^{2}/Hz)');
    hold off;
    xlim([.5 70]);  set(gca,'Box','off');
    YLIM = get(gca,'YLim');
    set(gca,'XTick', [0.5,4,8,13,30,60]);
    set(gca, 'XTickLabel', {'','4','8','13','30','60'});
end


myVisArt = vis_artifacts(tempEEG,EEG,'ChannelSubset',[xChanPlot yChanPlot]);


end

