%% Set paths
projdir     = 'Z:\mai\projects\shapesStory';
cd(projdir)
addpath(genpath(fullfile(projdir, 'code')));

% Get analysis params
which_group = 'mov2tom';
iscType     = 'within';
params      = get_analysisParams('fmri_group1');

statsdir = 'Z:\mai\projects\shapesStory\between\analysis\iscWholeBrain\1x2\half2_smooth4mm_motionAudioResid';

%% Max stat procedure
niters = 1000;
maxstat = NaN(1, niters);

fprintf(['calculating max stat for ' which_group '\n']);

% get null files
nullfiles = dir(fullfile(statsdir,  [which_group '_' '*']));
if length(nullfiles) ~= niters
    error(['num nullfiles = ' num2str(length(nullfiles)) ' ~= niters'])
end

% Get max
for j = 1:length(nullfiles)
    if mod(j,100)==0
        fprintf(['null file ' num2str(j) '\n']);
    end
    
    % load bootstrap
    boot = load(fullfile(statsdir,  nullfiles(j).name));

    maxstat(j) = max(boot.bootstrap);
end




%% stats


% fit normal curve
mu = mean(maxstat);
sigma = std(maxstat);

% percentile
crit1 = prctile(maxstat, 95);
crit2 = prctile(maxstat, 99);
crit3 = prctile(maxstat, 99.9);

fprintf('\n\n**** Max stat procedure ****\n');
fprintf(['alpha = .05, crit = ' num2str(crit1) '\n']);
fprintf(['alpha = .01, crit = ' num2str(crit2) '\n']);
fprintf(['alpha = .001, crit = ' num2str(crit3) '\n']);







