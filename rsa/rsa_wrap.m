
function rsa_wrap(taskid)
% Wrapper function for conducting intersubject RSA. Taskid refers to array
% job number. Number of jobs = nVoxels/1000.


%% Set paths
if ispc
    projdir     = '';
    fastscratch = '';
    taskid      = 71;
else
    projdir     = '';
    fastscratch = '';
end

% Set paths to hasson dir
codedir     = fullfile(projdir, 'code');
addpath(genpath(codedir));

%% Set params
rsatype = 'pairwise';
isctype = 'within';

group1 = 'fmri_group1';
group2 = 'fmri_group1';
cubefile = 'cubecount_gap_1_rad_0_gray';
savedir = 'gap1_rad0_gray';
suffix = '_motionResid';

if strcmp(isctype, 'within')
    % Which groups?
    params = get_analysisParams(group1);
    params2 = get_analysisParams(group2);
    
    %Keep which subs?
    params.keepsubs1 = 1:36;
    params.keepsubs2 = 1:36;
    params.subs = params.subs(params.keepsubs1);
    params2.subs = params2.subs(params.keepsubs2);
    
    % Which condition?
    params.cond = 'shapesMovie_smooth4mm2';
    params2.cond = 'shapesMovie_smooth4mm2'; 
    
    % Which LSA matrix?
    params.lsa_grp = 'grp1';
    
    % save dir?
    params.outputdir = fullfile(fastscratch, 'data_out', 'searchlt', savedir);
    
    % which isctype
    params.isctype = isctype;
    
else
    % Which groups?
    params = get_analysisParams(group1);    
    params2 = get_analysisParams(group2);
    
    %Keep which subs?
    params.keepsubs1 = 1:36;
    params.keepsubs2 = 1:36;
    params.subs = params.subs(params.keepsubs1);
    params2.subs = params2.subs(params.keepsubs2);
    
    
    % Which condition?
    params.cond = 'shapesMovie_smooth4mm2';
    params2.cond = 'shapesMovie_smooth4mm2';

    % Which LSA matrix?
    params.lsa_grp = 'grp1';
    
    % savedir
    params.outputdir = fullfile(fastscratch, 'data_out', 'searchlt', savedir);
    
    % which isctype
    params.isctype = isctype;
    
end

% what kind of RSA
params.rsatype = rsatype;

% Which directories?
params.datadir1 = fullfile(fastscratch, 'data_in', params.cond);
params.datadir2 = fullfile(fastscratch, 'data_in', params2.cond);

% Which cube file?
params.cubefile = fullfile(fastscratch, 'data_in', 'searchlt', 'cubefiles', cubefile);

% Which LSA file?
params.lsafile = fullfile(fastscratch, 'data_in',  'lsa', 'lsa_data2.mat');

% what crop?
params.crop = [60 305];

% number cubes/group
params.nCubes = 1000;

% which datafiles?
for i = 1:length(params.subs)
    params.datafiles1{i} = fullfile(params.datadir1, [params.subs{i} '_' params.cond suffix]);
end
for i = 1:length(params2.subs)
    params.datafiles2{i} = fullfile(params.datadir2, [params2.subs{i} '_' params2.cond suffix]);
end


fprintf(['**** LSA searchlight-' params.isctype ' ****\n']);
fprintf(['task ID = ' num2str(taskid) '\n']);
sprintf('cube file: %s', params.cubefile);
sprintf('lsa file: %s', params.lsafile);

fprintf('\n');
params


%% Go
% check if file has already been made
if exist(fullfile(params.outputdir, ['searchlt_rsa_shapesMovie_smooth4mm2_within_avgOthers_ss36_' num2str(taskid) '.mat']), 'file')
    fprintf('file exists already, exiting');
    return;
end


rsa(params, taskid);










