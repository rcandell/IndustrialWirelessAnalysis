% MeasurementRun.m
% Purpose: Define a MATLAB class to manage RF measurement run data and analysis.
classdef MeasurementRun < handle
    
    properties(Constant)
        OPT_PATH_GAIN = 1;
        OPT_KFACTOR = 2;
        OPT_DELAY_SPREAD = 3;
        OPT_AVGCIR = 4;
        OPT_NTAP_APPROX = 5;
        OPT_WRITE_STATS = 6;
        OPT_DO_PLOTS = 7;     

        NtapApprox_N = 24;  

    end
    
    properties

        OPTS = zeros(7,1);

        AI_METRICS_FNAME_OUT = 'aimetrics.csv';
        AI_REDTAP_FNAME_OUT = 'airedtap.csv';
        AI_FFT_FNAME_OUT = 'airedffts.csv';

        manifest_tbl = [];

    end
    
    methods
        
        function obj = MeasurementRun(OPTS, manifest_tbl)
            obj.manifest_tbl = manifest_tbl;
            if nargin < 1
                obj.OPTS = [ 1; 1; 1; 1; 1; 0; 0];
            else
                obj.OPTS = OPTS;
            end
        end

    end
end

