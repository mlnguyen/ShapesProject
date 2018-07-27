
function audenv = calc_audenv(audfile, params)
% function audenv = calc_audenv(audfile, params)
% Calculate audio envelop of audfile


%% Set default params

% Set resampling factor (>1 to downsample)
if ~isfield(params, 'resamp_factor')
    params.resamp_factor = 1; end

% Set TR length
if ~isfield(params, 'tr')
    params.tr = 1.5; end

% Set hrf lag
if ~isfield(params, 'tshift')
    params.tshift = 0; end

% Set fft filter
if ~isfield(params, 'freq')
    params.freq = []; end

% Convolve with hrf?
if ~isfield(params, 'doConvolve')
    params.doConvolve = 0; end;

% Downsample to bold resolution?
if ~isfield(params, 'mriResolution')
    params.mriResolution = 0; end;


%% Prep audio files

% Load audio file
[w, audfs] = audioread(audfile);

% Make mono
w = w(:,1);

% Crop, if specified
if params.crop ~= 0
    music_crop = params.crop*audfs;
    w = w(music_crop:end,1);
end

% Resample, if specified
if params.resamp_factor ~= 1
    w = w(:,1);
    w = resample(w, 1, params.resamp_factor);
    audfs = audfs/params.resamp_factor;
end

% Filter, if specified
if ~isempty(params.freq)
    w = filterfft(w, audfs, params.freq);
end


%% Get audio envelope

% Audio envelope
stim_mod = abs(hilbert(w));

% Stim length
stim_time = 0 + ( 0 : (length(w)-1) ) / audfs;

if params.doConvolve
    fs = (audfs/params.resamp_factor);
    dt = 1/fs;
    t = 0:dt:dt*(length(w) - 1); %time base
    hrf_Glover= @(t) t.^6 .* exp( -t/.9)/61.46 - 0.35*t.^12 .* exp( -t/.9)/1.54e7;
    hrf_timebase = 0:dt:30;
    hrf_t = hrf_Glover(hrf_timebase);
    
    stim_conv = conv2(hrf_t, stim_mod', 'full');
    audpow = stim_conv(1:length(t))'; %crop the padded end due to the convolution
end

% Crop to stimulus length
boldcrop = and(stim_time > params.stimLength(1), stim_time < params.stimLength(2)); % 304
audpow = stim_mod(boldcrop);

% Resample to tr length
audpow_slow = resample(audpow,1,(audfs/params.resamp_factor)*params.tr); 

% Shift by hrf lag
audpow_slow = round(1000+100*zscore(audpow_slow));
audpow_slow_shifted_100 = round(1000+100*zscore([audpow_slow(1:params.tshift); audpow_slow(1:end-params.tshift)]));
audenv = audpow_slow_shifted_100;


%% save
[stimdir, name] = fileparts(audfile);
savename = fullfile(stimdir, [name '_audEnv.mat']);
save(savename, 'audenv', 'audfs', 'params');


