function rsa_permute_wrap(task_id)


%% Set paths

if ispc
    projdir     = '';
    fastscratch = '';
    task_id      = 1;
else
    projdir     = '';
    fastscratch = '';
end

% Add code dir
codedir     = fullfile(projdir, 'code');
addpath(genpath(codedir));

%% Set params

params.isctype = 'within';
params.rsatype = 'pairwise';
params.iters = 1000;
params.rsa_dir = 'gap1_rad0_gray';
params.name = 'shapesMovie_smooth4mm2';

params.filename = ['searchlt_rsa_' params.name '_' params.isctype '_' params.rsatype '_ss36_'];
params.outputdir = fullfile(fastscratch, 'data_out', 'searchlt', params.rsa_dir); 

for i = 1:119
    params.datafiles{i} = fullfile(fastscratch, 'data_in', 'searchlt', ...
        [params.rsa_dir], [params.filename num2str(i) '.mat']);
end

% Seed rng
new_seed = sum(100*clock)+task_id;
fprintf(['rng seed: ' num2str(new_seed) '\n']);
rng(new_seed);

%% Go!
tic;
fprintf(['taskid: ' num2str(task_id) '\n']);
sprintf('data file: %s', params.datafiles{task_id});

params

rsa_searchlt_permute(params, task_id);

fprintf('\n *** job completed ***');
fprintf(['\n t = ' num2str(toc)]);
