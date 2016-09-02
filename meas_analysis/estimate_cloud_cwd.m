function estimate_cloud_cwd(pattern, OPTS, TEST_DATA)
% Analyze complex impulse responses from cloud measurements
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

TESTING = false;
if nargin == 3
    TESTING = true;
end

OPT_PATH_GAIN = 1;
OPT_DELAY_SPREAD = 3;
OPT_AVGCIR_NTap = 4;
if nargin < 2
    OPTS = ...
        [ ...   
            1; ...  % compute path gain
            1; ...  % delay spread
            1; ...  % compute average CIR from data
        ]; 
end

if nargin < 1
    error 'specify a file search pattern'
end

more off;

%
% query the list of measurement files
%
arr_dir = '.';          % sub directory of stored cir records
files = dir(pattern);

%
% setup data for stats text file 
%
Cstats = {};
Cstats_ii = 1;

%
% Process each file in turn
%
cir_file = [];
if TESTING
    Nfiles = 1;
else
    Nfiles = length(files);
end
for fk = 1:Nfiles
    
    % use test data or the real thing
    try 
        
        mat_fname = files(fk).name; 
        cir_file_path = [arr_dir '\' mat_fname];       
        if TESTING
            disp('opening TEST_DATA');
            cir_file = TEST_DATA;
        else
            disp(['loading file ' mat_fname '  ...']);
            cir_file = load(cir_file_path);
        end
        
    catch me
        warning('Problem reading mat file, trying again then skipping.');
        disp(me.message)         
        warning('Skipping file...');
        continue;
    end
    
    try
        meta = cir_file.Strct_Metadata;
    catch me
        warning('problem reading meta data')
        disp(me.message);
        continue;
    end
    
    % META DATA SECTION
    disp(meta)
    Ts = (1/meta.SampleRate_MHz_num)*1e-6;  % sample rate
    wl =  meta.CodewordLength_num;          % codeword length
    ns = meta.PNOversample_num;             % oversample rate
    t = (0:Ts:Ts*(wl-1));% -Ts*pn_over;     % time array over a burst transmission
    
    %NN = apf*rpa;
    NN = size(cir_file.IQdata,2);

    % Initialize memory for the metrics
    USE = nan(NN,1);  
    K = nan(NN,1);  
    LOS = nan(NN,1);  
    path_gain_dB = nan(NN,1);  
    path_gain_dB_poly = 0;
    rms_delay_spread_sec = nan(NN,1);     
    mean_delay_sec = nan(NN,1); 
    cir_duration = nan(NN,1);    
    
    % Memory for calculation of average cir
    klos = 0; 
    wla = 2*wl+1;   % size of cir avg calculation buffer
    mlos = ceil(wla/2);  %mid-point of cir avg calculation buffer
    cir_sum_los = zeros(wla,1);
    num_los = 0;
    cir_sum_nlos = zeros(wla,1);
    num_nlos = 0;
    
    cir_class = {'los','nlos'};
    for cir_class_ii = 1:length(cir_class)
        cir_avg_st(cir_class_ii).class = cir_class; %#ok<*AGROW>
        cir_avg_st(cir_class_ii).time = [];
        cir_avg_st(cir_class_ii).mag = [];
        cir_avg_st(cir_class_ii).angle = [];
    end

    % setup output directories
    stats_dir = 'stats';
    fig_dir = 'figs';
    png_dir = 'png';
    mkdir(stats_dir);
    mkdir(fig_dir); 
    mkdir(png_dir); 

    %
    % ANTENNA DATA
    % 
    TransmitterAntennaGain_dBi = meta.TransmitterAntennaGain_dBi_num;
    ReceiverAntennaGain_dBi = meta.ReceiverAntennaGain_dBi_num;
    
    
    % plot the x,y positions
    plot(cir_file.IQdata_CloudLocations_m_num.xPositions, ...
        cir_file.IQdata_CloudLocations_m_num.yPositions)
    xlabel('xpos (m)')
    ylabel('ypos (m)')
    
    
    % 
    % Loop through all records within the file
    %
    for kk = 1:NN
        
        % position of the CIR data record
        POS = [ cir_file.IQdata_CloudLocations_m_num.xPositions(kk), ...
                cir_file.IQdata_CloudLocations_m_num.yPositions(kk) ];
        
        % extract the CIR for this record from the data file
        cir = cir_file.IQdata(:,kk);
        
        % compute the magnitude of the CIR samples
        cir_mag2 = abs(cir).^2;
        
        % ignore record if it is empty or the length is less than the
        % expected codeword length.  This indicates that something went
        % wrong with the instrumentation.
        if isempty(cir_mag2)
            continue;
        elseif length(cir_mag2) < wl
            continue;
        end

        % select the sample of the cir that meet threshold criteria
        % also compute the noise floor
        [k_sel, nf] = select_cir_samples(r, cir);
        if isempty(k_sel)
            continue
        else
            if length(k_sel) < 8
                continue
            end
        end
        USE(kk) = 1;
        
        % extract the data at selected indicies
        t_k = t(k_sel);
        cir_k = cir(k_sel);

        % Compute the path loss in the cir
        % note that the CIR contains antenna gains.  We must remove the
        % bulk antenna gains using the assumption that the gain is
        % applied equally to all multi-path components.  We know that
        % this is not the true case, but without ray-tracing it is the
        % only option.
        if OPTS(OPT_PATH_GAIN)
            path_gain_dB(kk) = compute_path_gain(cir_k, ...
                TransmitterAntennaGain_dBi, ...
                ReceiverAntennaGain_dBi);           
        end
        
        % Aggregate the sums for later computation of avg CIR
        % LOS and NLOS are considered as separate classes of CIR's
        if OPTS(OPT_AVGCIR_NTap)
            [K(kk), LOS(kk), klos] = compute_k_factor(t_k, cir_k, ns);
            inds = mlos-klos+1:wla-klos;
            num_los = num_los + 1;
            cir_sum_los(inds) = cir_sum_los(inds) + cir;
        end        

        % compute delay spread parameters of the CIR 
        % because of the wrapping of energy in the FFT-based
        % correlation we must remove the trailing edge.
        if OPTS(OPT_DELAY_SPREAD)
            [mean_delay_sec(kk), rms_delay_spread_sec(kk), cir_duration(kk)] = ...
                compute_delay_spread(t_k, cir_k, nf);
        end

    end
    
    
    %
    % determine if the run produced enough data to form estimates.  The
    % selection of threshold was chosen arbitrarily
    %
    if sum(~isnan(USE)) < 100
        warning 'not enough CIRs passed selection to form metrics'
        disp(mat_fname)
        continue
    end
    
    %
    % Extract range data for off-line analysis
    %
    r = cir_file.IQdata_Range_m(:,3);
    path_gain_dB = path_gain_dB(r~=r(1));
    r = r(r~=r(1));
    
    if OPTS(OPT_PATH_GAIN)
    %
    % Analyze path loss versus distance
    %
    r_p = r;  pl_p = path_gain_dB;
    r_p = r_p(~isnan(pl_p));
    pl_p = pl_p(~isnan(pl_p));
    r_min_fit = 10;
    r_max_fit = 0.9*max(r_p);    
    k_lfit = find(r_p>r_min_fit & r_p<r_max_fit);
    r_p_fit = r_p(k_lfit);
    pl_p_fit = pl_p(k_lfit);
    path_gain_dB_poly = polyfit(r_p_fit,pl_p_fit,1);
    
    h = figure();
    stdPathGain = std(path_gain_dB(~isnan(path_gain_dB)));
    r_p_plot = linspace(min(r_p_fit), max(r_p_fit), 10);
    pl_poly_vals = polyval(path_gain_dB_poly, r_p_plot);
    plot(r_p, pl_p, 'color', [0,0,0]+0.7, 'marker', '.', 'linestyle' , 'none'); 
    hold on
    plot(r_p_plot, pl_poly_vals, 'k-', ...
       r_p_plot, repmat(pl_poly_vals(:),1,2)+stdPathGain*[ones(10,1) -ones(10,1)], ...
        'k--', 'LineWidth', 1.0);
    hold off
    legend({'measured', ...
        sprintf('%0.2fx + %0.1f',path_gain_dB_poly), ...
        '+/- \sigma'}, 'Location', 'best');
    setCommonAxisProps()    
    xlabel('distance (m)')
    ylabel('Path Gain (dB)')
    %title({'Path Loss', strrep(mat_fname,'_','-')})
    drawnow
    savefig(h, [fig_dir '\' mat_fname(1:end-4) '__pl.fig']);
    setFigureForPrinting();
    print(h,[png_dir '\' mat_fname(1:end-4) '__pl.png'],'-dpng','-r300')
    close(gcf)         
    
    end % if OPTS(OPT_PATH_GAIN)
    

    if OPTS(OPT_DELAY_SPREAD) && sum(~isnan(mean_delay_sec)) > 100
        %
        % Analyze the average delay of the CIR's 
        %
        h = figure();      
        [du_counts,du_centers] = hist(1e9*mean_delay_sec,50000);
        du_probs = cumsum(du_counts/sum(du_counts));
        plot(du_centers, du_probs, 'k');
        ylim([0 1]);
        xlim([0 max(du_centers(du_probs<0.995))]);
        str = 'average delay, $$\tau_D$$ (ns)';xlabel(str,'Interpreter','Latex')
        str = 'Pr. $$ \hat{\tau_D} < {\tau_D} $$'; ylabel(str,'Interpreter','Latex'); 
        setCommonAxisProps();
        %title({'CDF of Average Delay', strrep(mat_fname,'_','-')})
        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__du.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__du.png'],'-dpng','-r300')    
        close(gcf)

        %
        % Analyze the delay spread of the CIR's 
        %
        h = figure();      
        [ds_counts,ds_centers] = hist(1e9*rms_delay_spread_sec,5000);
        ds_probs = cumsum(ds_counts/sum(ds_counts));
        plot(ds_centers, ds_probs, 'k');
        ylim([0 1]);
        %xlabel('rms delay spread, S (ns)')
        str = 'rms delay spread, $$S$$ (ns)';xlabel(str,'Interpreter','Latex')
        str = 'Pr. $$\hat{S} < S$$'; ylabel(str,'Interpreter','Latex');
        setCommonAxisProps();
        %title({'CDF of Delay Spread', strrep(mat_fname,'_','-')})
        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__ds.fig']);
        h.PaperPositionMode = 'auto';
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__ds.png'],'-dpng','-r300')    
        close(gcf)
        
    end %if OPTS(OPT_DELAY_SPREAD)

    if OPTS(OPT_KFACTOR)
    end % OPTS(OPT_KFACTOR)
    
    % approximate an N-tap CIR from the measured CIR's
    if OPTS(OPT_AVGCIR_NTap)
        
        NtapApprox_N = 13;
        
        cir_class_name = cir_class(cir_class_ii);
        h = figure(); 
        cir_avg = cir_sum_los/num_los;

        % now remove leading zeros
        cir_avg(1:find(cir_avg, 1,'first')-1) = [];
        cir_avg = cir_avg(1:wl);
        cir_avg = select_for_avg_cir(cir_avg);
        cir_avg = cir_avg/max(abs(cir_avg));
        Ncir_avg = length(cir_avg);
        t_ciravg = t(1:Ncir_avg);
        cir_avg = cir_avg(1:Ncir_avg);
        [r_t,r_h,r_ph] = reduce_taps(cir_avg,NtapApprox_N);
        r_h = r_h/max(r_h);  % normalize the approximated cir

        hold off
        subplot(4,1,1:2)
        stem(1E9*t_ciravg, abs(cir_avg), '+-'); 
        str = '$$\mid{h(t)}\mid$$';ylabel(str, 'Interpreter', 'Latex')
        hold on; 
        stem(1E9*t_ciravg(r_t+1), abs(r_h),'d-');
        set(gca,'XTickLabel','')
        xlim([0 1000]); 
        hold off
        setCommonAxisProps();
        set(gca,'OuterPosition',get(gca,'OuterPosition').*[1 1 0.95 0.95]+[0.05 0.05 0 0])

        subplot(4,1,3:4)
        stem(1E9*t_ciravg, angle(cir_avg), '+-');
        str = '$$\angle{h(t)}$$';ylabel(str, 'Interpreter', 'Latex')
        xlabel('time (ns)')
        set(gca,...
             'ylim',[-2*pi() 2*pi()],...
             'ytick',[-2*pi() 0 2*pi()],...
             'yticklabel',{'-2\pi' '0' '2\pi'})
        xlim([0 1000]); 
        hold on
        stem(1E9*t_ciravg(r_t+1), r_ph, 'r')
        hold off
        setCommonAxisProps();
        set(gca,'OuterPosition',get(gca,'OuterPosition').*[1 1 0.95 0.95]+[0.05 0.05 0 0])

        cir_avg_st(cir_class_ii).time = t_ciravg;
        cir_avg_st(cir_class_ii).mag = r_h;
        cir_avg_st(cir_class_ii).angle = r_ph;     

        drawnow
        savefig(h,[fig_dir '\' mat_fname(1:end-4) '__avgcir_' cell2mat(cir_class_name) '.fig']);
        setFigureForPrinting();
        print(h,[png_dir '\' mat_fname(1:end-4) '__avgcir_' cell2mat(cir_class_name) '.png'],'-dpng','-r300')    
        close(gcf)   
        
    end % OPTS(OPT_AVGCIR_NTap)
    
    % save the metrics
    stats = struct(...
        'meta',meta,...
        'path_gain_dB',path_gain_dB,...
        'rms_delay_spread_sec',rms_delay_spread_sec, ...
        'mean_delay_sec',mean_delay_sec, ...
        'cir_duration_sec',cir_duration, ...
        'avg_cir_st', cir_avg_st);
    
    save([stats_dir '\' mat_fname(1:end-4) '__channel_stats.mat'], 'stats')

    % save the stats for this measurement run
    RxPol = 'U';
    if strfind(meta.ReceiverAntenna_str,'V Pol')
        RxPol = 'V';
    elseif ~isempty(strfind(meta.ReceiverAntenna_str,'Cross Pol')) || ...
            ~isempty(strfind(meta.ReceiverAntenna_str,'X Pol')) || ...
            ~isempty(strfind(meta.ReceiverAntenna_str,'Long Pol'))
        RxPol = 'X';
    end
    TxPol = 'U';
    if ~isempty(strfind(meta.TransmitterAntenna_str,'V Pol'))
        TxPol = 'V';
    elseif ~isempty(strfind(meta.TransmitterAntenna_str,'Cross Pol')) || ...
            ~isempty(strfind(meta.TransmitterAntenna_str,'X Pol')) || ...
            ~isempty(strfind(meta.TransmitterAntenna_str,'Long Pol'))
        TxPol = 'X';
    end    
    Cstats(Cstats_ii,:) = {   
        meta.MatFile_str, meta.Frequency_GHz_num, ...
        RxPol, meta.ReceiverAntennaGain_dBi_num, ...
        TxPol, meta.TransmitterAntennaGain_dBi_num, ...
        stats.path_gain_dB_poly(1), stats.path_gain_dB_poly(2), ...
        nanmean(stats.K), nanmin(stats.K), nanmax(stats.K) ...
        1e9*nanmean(stats.mean_delay_sec), 1e9*nanmin(stats.mean_delay_sec), 1e9*nanmax(stats.mean_delay_sec)...
        1e9*nanmean(stats.rms_delay_spread_sec), 1e9*nanmin(stats.rms_delay_spread_sec), 1e9*nanmax(stats.rms_delay_spread_sec) ...
    };
    Cstats_ii = Cstats_ii + 1;

    % explicit clear of large memory
    cir_file = [];  %#ok<NASGU>

end

% add entry to the stats text file
writeStatsToFile(Cstats);

%
% Create summary data
%

% create the aggregate polynomial for path loss
cmp_pl_poly( '*_stats.mat', '.\stats', '.\figs', '.\png' )

end % function

function setCommonAxisProps()

    alw = 0.75;    % AxesLineWidth
    fsz = 10;      % Fontsize
    lw = 1.5;      % LineWidth
    msz = 3.5;       % MarkerSize
    
%    grid on
    set(gca,'XGrid','on')
    set(gca,'XMinorGrid','off')
    set(gca,'YGrid','on')
    set(gca,'YMinorGrid','off')
    set(gca,'GridAlpha',0.5)
    set(gca,'MinorGridAlpha',0.4)
    set(gca,'Fontsize',fsz)
    set(gca,'LineWidth',alw);
    set(gca,'FontName','TimesRoman')
    
    % set the line properties
    hline = get(gca,'Children');
    for h = hline(:)'
        h.LineWidth = lw;
        h.MarkerSize = msz;
    end
end

function setFigureForPrinting()
    width=3; height=3;
    %set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
end

function writeStatsToFile(X)
    if isempty(X)
        return
    end
    M = cell2table(X);
    file_path = '../stats.dat';
    VariableNames = ...
        {'Run', 'Freq', 'RX_Ant_Type', 'RX_Ant_Gain', 'TX_Ant_Type', 'TX_Ant_Gain', ...
        'Path_Gain_Poly_Slope', 'Path_Gain_Poly_YInt', 'Mean_K', 'Min_K', 'Max_K', ...
        'Mean_Delay_ns', 'Min_Delay_ns', 'Max_Delay_ns', ...
        'Mean_Delay_Spread_ns', 'Min_Delay_Spread_ns', 'Max_Delay_Spread_ns'};
    M.Properties.VariableNames = VariableNames;
    
    M0 = [];
    if exist(file_path,'file')
        M0 = readtable(file_path);
    end
    if ~isempty(M0)
        M = [M0;M];
    end
    M.Properties.VariableNames = VariableNames;
    writetable(M, file_path, 'Delimiter', '\t');
    
end


