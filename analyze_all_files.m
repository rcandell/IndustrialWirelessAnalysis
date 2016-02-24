% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

clear;

doall = false;
figvis = false;

% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')

%addpath('c:\Users\rnc4\git\rf_analysis')
%dirs = { 'AAPlant D1', 'AAPlant D2', 'AAPlant D3' };
%pattern = 'AA*.mat';
dirs = { 'Gaithersburg Day 1', 'Gaithersburg Day 2', 'Gaithersburg Day 3', 'AAPlant D1', 'AAPlant D2', 'AAPlant D3' };
pattern = '*__array.mat';

dbstop error
parfor kk = 1:length(dirs)
    disp(['entering ' dirs{kk}])
    chdir(dirs{kk})
    %estimate_channel(pattern);
    estimate_channel_arr(pattern, doall, figvis);
    chdir('..')
    disp(['now in ' chdir()])
end



