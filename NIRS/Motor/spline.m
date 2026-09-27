% Legacy launcher for spline processing.
% The main pipeline now lives in CW_motor_processing_unified.m

motionCorrectionMethod = "spline";
run(fullfile(fileparts(mfilename('fullpath')), 'CW_motor_processing_unified.m'));
