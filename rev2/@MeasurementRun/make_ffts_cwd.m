function make_ffts_cwd(obj, pattern, TEST_DATA)
% Make fft records from complex impulse responses 
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

TESTING = false;
if nargin == 3
    TESTING = true;
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

% setup output directories
ffts_dir = 'ffts';
if ~exist(ffts_dir,'dir'), [status, msg, msgID] = mkdir(ffts_dir), end

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
        
        if TESTING
            disp('opening TEST_DATA');
            cir_file = TEST_DATA;
            mat_fname = [TEST_DATA.Strct_Metadata.MatFile_str '_pp.mat'];
        else
            mat_fname = files(fk).name;             
            cir_file_path = [arr_dir '\' mat_fname]; 
            if testForStatsFile(mat_fname(1:end-4)) % && obj.OPTS(obj.OPT_WRITE_STATS)
                disp(['skipping ' mat_fname])
                continue;
            end
            disp(['loading file ' mat_fname '  ...']);
            cir_file = load(cir_file_path);
        end
        
    catch me
        disp('Problem reading mat file, skipping...');
        disp(me.message)         
        continue;
    end
    
    try
        meta = cir_file.Strct_Metadata;
    catch me
        warning('problem reading meta data')
        disp(me.message);
        continue;
    end
    
    % enforce mat file str is correct, some are incorrect in the files
    meta.MatFile_str = mat_fname;
    
    % determine the site
    SITE_IND_OATS    = 0;
    SITE_IND_AAPLANT    = 1;
    SITE_IND_GBURG      = 2;
    SITE_IND_STEAM      = 3;
    if ~isempty(strfind(lower(mat_fname),'aaplant'))
        SITE_IND = SITE_IND_AAPLANT;
    elseif ~isempty(strfind(lower(mat_fname),'gburg'))
        SITE_IND = SITE_IND_GBURG;
    elseif ~isempty(strfind(lower(mat_fname),'oats'))
        SITE_IND = SITE_IND_OATS;
    else
        SITE_IND = SITE_IND_STEAM;
    end
    
    % META DATA SECTION
    %disp(meta)
    Fs = meta.SampleRate_MHz_num;
    Ts = (1/Fs);  % sample period
    wl =  meta.CodewordLength_num;          % codeword length
    ns = meta.PNOversample_num;             % oversample rate
    t = (0:Ts:Ts*(wl-1));% -Ts*pn_over;     % time array over a burst transmission
    
    %NN = apf*rpa;
    NN = size(cir_file.IQdata,2);
    
    % truncate the last section so not to emphasize the end of the run
    NN = round(NN*0.9); 

    % Initialize memory for the metrics
    USE = nan(NN,1);  
    % K = nan(NN,1);  
    % LOS = nan(NN,1);  
    % path_gain_dB = nan(NN,1);
    % path_gain_range_m = nan(NN,1);
    % path_gain_dB_poly = {};
    % rms_delay_spread_sec = nan(NN,1);     
    % mean_delay_sec = nan(NN,1); 
    % cir_duration = nan(NN,1);    
    
    % Memory for calculation of average cir
    % klos = 0; 
    % wla = 2*wl+1;   % size of cir avg calculation buffer
    % mla = floor(wla/2)+1;  %mid-point of cir avg calculation buffer
    % cir_sum_los = zeros(wla,1);
    % num_los = 0;
    % cir_sum_nlos = zeros(wla,1);
    % num_nlos = 0;

    %
    % ANTENNA DATA
    % 
    % TransmitterAntennaGain_dBi = meta.TransmitterAntennaGain_dBi_num;
    % ReceiverAntennaGain_dBi = meta.ReceiverAntennaGain_dBi_num;

    % metrics table 
    metrics_arr = [];
    ffts_cir_arr = [];
    m_kk = 1;

    % 
    % Loop through all records within the file
    %
    for kk = 1:50:NN

        % get the meta data from the manifest
        u = obj.manifest_tbl(strip(obj.manifest_tbl.filename) == strip(mat_fname), :);
        if isempty(u)
            error('data not found in manifest table')
        else
            metrics_arr{m_kk, MeasurementRunMetric.ID} = u.ID;
            metrics_arr{m_kk, MeasurementRunMetric.Site} = u.Location;
            metrics_arr{m_kk, MeasurementRunMetric.TxPos} = u.TxPos;
            metrics_arr{m_kk, MeasurementRunMetric.TxHeight} = u.TxHeight;
            metrics_arr{m_kk, MeasurementRunMetric.Freq} = u.Fc;
            metrics_arr{m_kk, MeasurementRunMetric.Polarization} = u.Pol;
            metrics_arr{m_kk, MeasurementRunMetric.Obstructed} = u.Obstructed;
            metrics_arr{m_kk, MeasurementRunMetric.Waveguided} = u.Waveguided;
        end        
        
        % range of the CIR data record
        r = cir_file.IQdata_Range_m(kk,3);
        
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
        [k_sel, ~, cir, pk_pwr] = obj.select_cir_samples(r, cir);
        if isempty(k_sel)
            continue
        end
        USE(kk) = 1;

        % record frequency
        metrics_arr{m_kk, MeasurementRunMetric.Freq} = meta.Frequency_GHz_num;        

        % record the range
        metrics_arr{m_kk, MeasurementRunMetric.Range_m} = r;   
        metrics_arr{m_kk, MeasurementRunMetric.Range_l} = physconst('LightSpeed')/r;      

        % record the X and Y coordinates
        coord_x = cir_file.IQdata_Range_m(kk,4);
        coord_y = cir_file.IQdata_Range_m(kk,5);
        metrics_arr{m_kk,MeasurementRunMetric.CoordX} = coord_x;
        metrics_arr{m_kk,MeasurementRunMetric.CoordY} = coord_y;    

        % capture the fft data for this cir
        fft_cir = fft(cir);
        % take only the central region
        D = 6;  % reduction each side of spectrum
        L = length(cir);
        Lcut = round(L/D);
        ii_central = [1:Lcut, L-Lcut+1:L];
        fft_central = fft_cir(ii_central);
        fft_central = fftshift(fft_central);
        
        % now resample the central region
        Ldes = 25;
        Ndes = round(Lcut*2/Ldes);
        fft_central = resample(fft_central, 1, Ndes);

        cir_fft_r = real(fft_central);
        cir_fft_i = imag(fft_central);
        ffts_cir_arr(m_kk,:) = [cir_fft_r(:).', cir_fft_i(:).'];

        if 0
            subplot(2,1,1), plot(abs(fftshift(fft_cir)))
            subplot(2,1,2), plot(abs(fft_central))
        end

        % write the fft real and imaginary to file
        var_names = MeasurementRunMetric.getNames();
        metrics_tbl = cell2table(metrics_arr, 'VariableNames', var_names(1:MeasurementRunMetric.getNumBaseColumns())); 
        red_fft_with_labels_arr = [metrics_tbl(:,1:MeasurementRunMetric.getNumBaseColumns()), num2cell(ffts_cir_arr)];
        ai_redfft_fname = [ffts_dir '/' obj.AI_FFT_FNAME_OUT];
        writetable(red_fft_with_labels_arr, ai_redfft_fname, 'WriteMode','append');          
               
        m_kk = m_kk + 1;
    end

    % explicit clear of large memory
    cir_file = [];  %#ok<NASGU>

end

end % function

function b = testForStatsFile(mat_fname)
    b = false;
    if exist(['stats/' mat_fname '__channel_stats.mat'],'file')
        b = true;
    end
end






