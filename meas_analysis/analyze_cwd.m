function analyze_cwd(OPTS, TEST_DATA)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

set(0,'DefaultFigureWindowStyle','docked')
set(0, 'DefaultFigureVisible', 'on')

TESTING = false;
if nargin == 2
    TESTING = true;
end

if nargin < 1 || isempty(OPTS)
    OPTS = ...
    [ ...   
        0; ...  % compute path gain
        0; ...  % K factor
        0; ...  % delay spread
        1; ...  % compute average CIR from data
        0; ...  % NTAP estimation
        0; ...  % write stats file
    ]; 
    disp('using the following options')
    disp(OPTS)
elseif strcmp(OPTS,'all')
    disp('enabling all options')
    OPTS = ...
    [ ...   
        1; ...  % compute path gain
        1; ...  % K factor
        1; ...  % delay spread
        1; ...  % compute average CIR from data
        0; ...  % ntap approx
        1; ...  % write stats file        
    ];     
end

pattern = '*.mat';
if TESTING
    estimate_channel_cwd(pattern, OPTS, TEST_DATA);
else
    estimate_channel_cwd(pattern, OPTS);
end


set(0, 'DefaultFigureVisible', 'on')

