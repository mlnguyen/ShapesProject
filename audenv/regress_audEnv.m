% regress_audEnv.m
% Script for regressing out audio envelop from fMRI timecourse

%% setup
if ispc
    projdir = 'Z:\mai\projects\shapesStory';
else
    projdir = '/jukebox/hasson/mai/projects/shapesStory';
end

cd(projdir)
addpath(genpath(fullfile(projdir, 'code')));

% Which group?
group = 'fmri_group3';
params = get_analysisParams(group);

% crop
params.crop = [1 305];

condN = 3;

% aud envs
auddir = '';
audfile = 'phys_audio_wIntro2_audio_audEnv';

% Load kept vox files
keptMov = load('');
keptPhys = load('');

% get intersection
kept_sum = keptMov.keptVox + keptPhys.keptVox;
keptAll = zeros(size(kept_sum));
keptAll(kept_sum==2) = 1;


%% Do regression

% load audio file & crop
aud = load(fullfile(auddir, audfile));

% convolve aud env
h = hrf('twogamma', 1.5);
audenv = conv(aud.audenv, h);
audenv = zscore(audenv(params.crop(1):params.crop(2)));

fprintf(['\naudfile: ' audfile '...\n']);
params.scans = {'physAudio_smooth6mm'};

for i = 1:length(params.subs)
    
    fprintf(['processing ' params.subs{i} '...\n']);
    
   % load data
   datafile = fullfile(params.datadir, params.scans{1}, [params.subs{i} '_' params.scans{1} '.mat']);
   data = load(datafile);
   tc = data.tc(:,params.crop(1):params.crop(2))';
   
   % regress
   b = NaN(1,length(tc));
   for j = 1:length(tc)
       b(j) = regress(tc(:,j), audenv);
   end
   
   % fit
   tc_fit =  audenv * b;
   resid = tc - tc_fit;
   
   %mask
   brain = nan(length(keptAll), size(resid,1));
   if condN == 1
       brain(keptMov.keptVox==1,:) = resid';
   elseif condN == 2
       brain(keptAud.keptVox==1,:) = resid';
   elseif condN ==3
       brain(keptPhys.keptVox==1,:) = resid';
   end
   % mask with new brain mask
   data_masked = brain(keptAll==1,:);
   
   % save
   data.tc = data_masked;
   data.keptVox = keptAll;
   savename = fullfile(params.datadir, params.scans{1}, [params.subs{i} '_' params.scans{1} '_resid.mat']);
   save(savename, '-struct', 'data');
   
end

