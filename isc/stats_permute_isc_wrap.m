% stats_bootstrap_isc_wrap

function stats_permute_isc_wrap(taskid)
% Wrapper function for calculating significance of ISC. This function will
% run x iterations and takes a taskid indicating job number for
% parallelizing on a cmoputer cluster. This fx can also run locally.

%% Set paths

if ispc
    projdir     = 'localpath';
    fastscratch = 'localpath';
    taskid      = 6;
else
    projdir     = 'clusterpath';
    fastscratch = 'clusterpath';
end


% Set paths to hasson dir
codedir     = fullfile(projdir, 'code');
addpath(genpath(codedir));


%% Set analysis params

% Get analysis params
which_group1            = 'shapesMovie_smooth4mm2';
which_group2            = 'tomAudio_smooth4mm2';
params                  = get_analysisParams('fmri_group1');
params2                 = get_analysisParams('fmri_group2');

params.subs             = params.subs(keep);
params.isc.condname     = 'mov2tom_half2';
params.isc.type         = 'between';
params.isc.conds        = [1 2];
params.isc.permute      = 0;
params.iterations       = 4;
params.crop             = [60 305];
params.scans            = {'shapesMovie_smooth4mm2_motionAudioResid'  ...
                            'tomAudio_smooth4mm2_motionAudioResid'};

params.datadir1         = fullfile(fastscratch, 'data_in', which_group1);
params.datadir2         = fullfile(fastscratch, 'data_in', which_group2);
params.savedir          = fullfile(fastscratch, 'data_out', params.isc.condname);
params.savename         = fullfile(params.savedir, ...
                            [params.isc.condname '_' params.isc.type '_motionAudResid_' num2str(taskid) '-']);
params.maskfile1        = fullfile(fastscratch, 'data_in', which_group1, ...
                            ['keptVox_all_' which_group1 '_n36' ]);
params.maskfile2        = fullfile(fastscratch, 'data_in', which_group2, ...
                            ['keptVox_all_' which_group2 '_n18' ]);                        
params.sumfile          = fullfile(fastscratch, 'data_in', which_group1, ...
                            ['sum_all_' params.scans{1} '_half2.mat']);

% Seed random number generator
new_seed = sum(100*clock)+taskid;
fprintf(['rng seed: ' num2str(new_seed) '\n']);
rng(new_seed);

% Make savedir if does not exist
if ~exist(params.savedir, 'dir')
    mkdir(params.savedir);
end

%% Go!
tic;


fprintf(['taskid: ' num2str(taskid) '\n']);
stats_permute_isc(params, params.subs, params.isc.conds, params.isc.permute);

fprintf('\n *** job completed ***');
fprintf(['\n t = ' num2str(toc)]);


