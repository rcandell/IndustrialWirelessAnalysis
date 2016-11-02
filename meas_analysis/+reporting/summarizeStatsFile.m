function summarizeStatsFile(filename)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if nargin < 1
    filename = 'stats.dat';
end

% ACTIVE STATS
%diary('summary_stats.txt')
stats = active_stats(filename);
stats.RunType = cell(size(stats,1),1);
stats.Site = cell(size(stats,1),1);
for ii = 1:size(stats,1)
    if strfind(lower(cell2mat(stats.Run(ii))),'aaplant')
        stats.Site(ii) = {'aaplant'};
    elseif strfind(lower(cell2mat(stats.Run(ii))),'gburg')
        stats.Site(ii) = {'shop'};
    elseif strfind(lower(cell2mat(stats.Run(ii))),'oats')
        stats.Site(ii) = {'oats'};        
    else
        stats.Site(ii) = {'steam'};
    end
    if strfind(cell2mat(stats.Run(ii)),'internal')
        stats.RunType(ii) = {'internal'};
    else
        stats.RunType(ii) = {'external'};
    end
end

ds = table2dataset(stats);

disp('MOBILE MEASUREMENTS =============================')
disp('=================================================')
disp('=================================================')
disp(' ')

disp('ALL ====================================')
% matrix2latex(ds(:,[1 2 18 3 5 7 12 15 9]), '../pub/alldata.tex', 'columnLabels', ...
%     {'Run','$f$ (GHz)','Route','RX Ant','TX Ant','Slope (dB/dec)', ...
%         '$\bar{\tau}$ (ns)','$\bar{S}$ (ns)','$\bar{K}$ (dB)'}, ...
%         'format', '%-.1f', 'alignment', 'c');

statarray = grpstats(ds,{'Run'},'mean','DataVars',{'Freq','Path_Gain_Poly_Slope','Path_Gain_Poly_YInt', 'Mean_Delay_ns', 'Mean_Delay_Spread_ns', 'Mean_K','Max_K','Min_K'})
for ii = 1:length(statarray.Run)
    statarray.Run{ii} = strrep(statarray.Run{ii},'_',' ');
end
export(statarray,'XLSfile','all_data.xls')
statarray1 = statarray(1:end,[1 3 4 6 7 8]);
N = length(statarray1);
N2 = round(N/2);
matrix2latex(statarray1(1:N2,:), '../pub/allruns1.tex', 'columnLabels', ...
    {'Run','$f$ (GHz)','Slope (dB/dec)', ...
        '$\bar{\tau}$ (ns)','$\bar{S}$ (ns)','$\bar{K}$ (dB)',...
        'Max K (dB)'}, 'format', '%-.1f', 'alignment', 'c');
matrix2latex(statarray1(N2+1:N,:), '../pub/allruns2.tex', 'columnLabels', ...
    {'Run','$f$ (GHz)','Slope (dB/dec)', ...
        '$\bar{\tau}$ (ns)','$\bar{S}$ (ns)','$\bar{K}$ (dB)',...
        'Max K (dB)'}, 'format', '%-.1f', 'alignment', 'c');
    
statarray = grpstats(ds,{'Site','Freq','RunType'},'mean','DataVars',{'Path_Gain_Poly_Slope','Path_Gain_Poly_YInt', 'Mean_Delay_ns', 'Mean_Delay_Spread_ns', 'Mean_K','Max_K','Min_K'})
export(statarray,'XLSfile','mobile_gain.xls')
matrix2latex(statarray(1:end,[1 2 3 5 6]), '../pub/all.tex', 'columnLabels', ...
    {'Site','$f$ (GHz)','Route','Slope (dB/dec)','Intercept (dB)', ...
        '$\bar{\tau}$ (ns)','$\bar{S}$ (ns)','$\bar{K}$ (dB)',...
        'Max K (dB)','Min K (dB)'}, 'format', '%-.1f', 'alignment', 'c');


disp('CHANNEL GAIN  ====================================')
statarray = grpstats(ds,{'Site','Freq','RunType'},'mean','DataVars',{'Path_Gain_Poly_Slope','Path_Gain_Poly_YInt'})
export(statarray,'XLSfile','mobile_gain.xls')
matrix2latex(statarray(1:end,[1 2 3 5 6]), '../pub/gain.tex', 'columnLabels', {'Site','$f$ (GHz)','Route','Slope (dB/dec)','Intercept (dB)'}, 'format', '%-.1f', 'alignment', 'c')

% disp('POWER DELAY  ====================================')
% statarray = grpstats(ds,{'Site','Freq','RunType'},'mean','DataVars',{'Mean_Delay_ns'})
% export(statarray,'XLSfile','mobile_du.xls')
% matrix2latex(statarray(1:end,[1 2 3 5]), 'delay.tex', 'columnLabels', {'Site','f (GHz)','Route','$\bar{\tau}$ (ns)'}, 'format', '%-.1f', 'alignment', 'c') 

disp('DELAY AND DELAY SPREAD ====================================')
statarray = grpstats(ds,{'Site','Freq','RunType'},'mean','DataVars',{'Mean_Delay_ns', 'Mean_Delay_Spread_ns'})
export(statarray,'XLSfile','mobile_ds.xls')
matrix2latex(statarray(1:end,[1 2 3 5 6]), '../pub/spread.tex', 'columnLabels', {'Site','$f$ (GHz)','Route','$\bar{\tau}$ (ns)','$\bar{S}$ (ns)'}, 'format', '%-.1f', 'alignment', 'c') 

disp('RICIAN K FACTOR==================================')
statarray = grpstats(ds,{'Site','Freq','RunType'},'mean','DataVars',{'Mean_K','Max_K','Min_K'})
export(statarray,'XLSfile','mobile_K.xls')
matrix2latex(statarray(1:end,[1 2 3 5 6 7]), '../pub/kfactor.tex', 'columnLabels', {'Site','$f$ (GHz)','Route','$\bar{K}$ (dB)','Max K (dB)','Min K (dB)'}, 'format', '%-.1f', 'alignment', 'c') 

return
% CLOUD STATS
cloudstats = cloud_stats();
cloudstats.RunType = cell(size(cloudstats,1),1);
cloudstats.Site = cell(size(cloudstats,1),1);
for ii = 1:size(cloudstats,1)
    if strfind(cell2mat(cloudstats.Run(ii)),'internal')
        cloudstats.RunType(ii) = {'internal'};
    else
        cloudstats.RunType(ii) = {'external'};
    end
    cloudstats.Site(ii) = {'steam'};
end

ds_c = table2dataset(cloudstats);

disp('CLOUD MEASUREMENTS ==============================')
disp('=================================================')
disp('=================================================')
disp(' ')

disp('POWER DELAY  ====================================')
statarray = grpstats(ds_c,{'Site','Freq','RunType'},'mean','DataVars',{'Mean_Delay_ns','Max_Delay_ns','Min_Delay_ns'})
export(statarray,'XLSfile','cloud_du.xls')

disp('DELAY SPREAD ====================================')
statarray = grpstats(ds_c,{'Site','Freq','RunType'},'mean','DataVars',{'Mean_Delay_Spread_ns','Max_Delay_Spread_ns','Min_Delay_Spread_ns'})
export(statarray,'XLSfile','cloud_ds.xls')


% diary off

end

function stats = active_stats(filename)
%% Import data from text file.
% Script for importing data from the following text file:
%
%    E:\FACTORY_MEASUREMENTS2\active_meas\stats.dat
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2016/09/29 13:13:07

%% Initialize variables.
delimiter = '\t';
startRow = 2;

%% Format string for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: text (%s)
%	column4: double (%f)
%   column5: text (%s)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%s%f%s%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
stats = table(dataArray{1:end-1}, 'VariableNames', {'Run','Freq','RX_Ant_Type','RX_Ant_Gain','TX_Ant_Type','TX_Ant_Gain','Path_Gain_Poly_Slope','Path_Gain_Poly_YInt','Mean_K','Min_K','Max_K','Mean_Delay_ns','Min_Delay_ns','Max_Delay_ns','Mean_Delay_Spread_ns','Min_Delay_Spread_ns','Max_Delay_Spread_ns'});

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

end

function cloudstats = cloud_stats(filename)
%% Import data from text file.
% Script for importing data from the following text file:
%
%    E:\FACTORY_MEASUREMENTS2\active_meas\Boulder_c\cloud_stats.dat
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2016/09/29 14:06:40

%% Initialize variables.
delimiter = '\t';
startRow = 2;

%% Format string for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: text (%s)
%	column4: double (%f)
%   column5: text (%s)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%s%f%s%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
cloudstats = table(dataArray{1:end-1}, 'VariableNames', {'Run','Freq','RX_Ant_Type','RX_Ant_Gain','TX_Ant_Type','TX_Ant_Gain','Mean_Delay_ns','Min_Delay_ns','Max_Delay_ns','Mean_Delay_Spread_ns','Min_Delay_Spread_ns','Max_Delay_Spread_ns'});

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

end
