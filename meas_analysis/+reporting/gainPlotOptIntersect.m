function [p1, p2] = gainPlotOptIntersect(pattern,freq, location_str, loc_logic, ~)
% Analyze complex impulse responses from measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

if nargin < 4
    loc_logic = true;
end
if nargin < 3
    location_str = NaN;
end
if nargin < 2
    error 'frequency specification is required'
end
if nargin <1
    error 'patter is required'
end

p1 = [];
p2 = [];

%
% query the list of stats files
%
files = dir(pattern);
Nfiles = length(files);

% local variables
R = [];
G = [];

%
% Process each file in turn
%
for fk = 1:Nfiles
    
    % use test data or the real thing
    try 
        
        stats_file_path = files(fk).name;      
        if strfind(stats_file_path,'Oats')
            continue;
        end
        stats = load(stats_file_path);
        stats = stats.stats;
        
    catch me
        warning('Problem reading mat file, trying again then skipping.');
        disp(me.message)         
        warning('Skipping file...');
        continue;
    end
    
    try
        meta = stats.meta;
    catch me
        warning('problem reading meta data')
        disp(me.message);
        continue;
    end
    
    % filter files
    if ~reporting.keepFile(meta, freq, location_str, loc_logic)
        continue;
    end
    
    disp(['Processing file: ' meta.MatFile_str])
    
    % path gain
    R = [R; stats.path_gain_range_m(:)]; %#ok<AGROW>
    G = [G; stats.path_gain_dB(:)]; %#ok<AGROW>
    
end

if isempty(R)
    return
end

% least squares fit
[p1, p2, x_intersect] = calLSqFit(R,G);

% plots 
plotThePathGains(meta, R,G, p1, p2, x_intersect);

end
    

function [p1, p2, x_intersect, gof] = calLSqFit(R,G)
    [fr, gof] = reporting.createPwLFit(log10(R), G);
    p1 = [fr.a1, fr.b1];
    p2 = [fr.a2, fr.b2];
    Q = fr.Q;
    %x_intersect = 10^(fr.Q);
    x_intersect = fzero(@(x) polyval(p1-p2,x),10^(fr.Q));
end

function h = plotThePathGains(meta, R,G, p1, p2, x_intersect)

    if 10^(x_intersect) < min(R)
        x_intersect = [];
    end

    if ~isempty(x_intersect)
        d_val1 = log10(logspace(log10(min(R)), x_intersect, 5));
        g_val1 = polyval(p1, d_val1); 
        d_val2 = log10(logspace(x_intersect,log10(max(R)),5));
        g_val2 = polyval(p2, d_val2);   
    else
        d_val2 = log10(logspace(log10(min(R)),log10(max(R)),5));
        g_val2 = polyval(p2, d_val2);          
    end

    % frii as reference
    d_frii = logspace(log10(min(R)), log10(max(R)), 5);
    ff = meta.Frequency_GHz_num;
    c = physconst('LightSpeed');
    frii_fspl_dB = 10*log10(d_frii.^2) + 20*log10(ff) + 20*log10(1e9) + 20*log10(4*pi/c);         

    % plot the gains
    h = figure();
    semilogx(R,G, 'color', [0,0,0]+0.7, 'marker', '.', 'linestyle' , 'none')
    hold on
    if ~isempty(x_intersect)
        semilogx(10.^d_val1,g_val1,'bx-')
    end
    semilogx(10.^d_val2,g_val2,'bo-')
    semilogx(d_frii, -frii_fspl_dB, 'k-') 
    xlim([0.75*min(R) 1.25*max(R)])
    hold off
    if ~isempty(x_intersect)
        legend('data',sprintf('fit1(n=%.1f)',abs(p1(1)/10)),sprintf('fit2(n=%.1f)', abs(p2(1)/10)),'fspl')
    else
        legend('data',sprintf('fit2(n=%.1f)', abs(p2(1)/10)),'fspl')
    end
    legend('Location','Best')

    reporting.setCommonAxisProps()    
    xlabel('distance (m)')
    ylabel('Path Gain (dB)')
    drawnow

end


    