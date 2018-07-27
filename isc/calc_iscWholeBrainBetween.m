function calc_iscWholeBrainBetween(params, iscParams, datafiles, avg_others, keptVox)
% Calculate ISC between two groups of subjects: each subject in group1 to
% average of all subjects in group2
%
% ARGUMENTS:
%  - params: struct specifying experiment params with fields - 
%        projdir: ''
%        datadir: ''
%         roidir: ''
%           subs: {1x18 cell}
%        exclude: {[]  []  []}
%          group: 'fmri_group1'
%          scans: {'shapesMovie_smooth4mm2'}
%          crop: [tr1 tr2]
% 
% - iscParams: struct specifying params for running ISC
%      savedir: ''
%         corr: [1 1] %which scans to compare
%         name: 'shapesMovie_smooth4mm2_half1' %name of resulting nifti
%         type: 'within' %type of ISC (within/between)
%         subs: {1x18 cell} %list of subjects
%     savename: ''
%
% - datafiles: list of .mat datafiles to run isc over. Each .mat file
%      contains a struct with the following fields: 
%           keptVox: [902629x1 logical] %binary mask with 1s for kept vox
%           datasize: [91 109 91] %original dims
%           tc: [177632x305 single] %nVoxKept x nTRs
%
% - avgOthers: .mat containing average of subjects in group2. Has field
%       avg_all
%

fprintf(['\n *** Calculating ISC (' iscParams.type '-group): ' iscParams.name '***\n']);

% load average of others
others = load(avg_others);
keptOthers = others.keptVox;

for i = 1:length(datafiles)
    
    fprintf([num2str(i) '-']);
    
    % load subdata
    dataStruct = load(datafiles{i});
    keptData = dataStruct.keptVox;
    
    % crop
    data = dataStruct.tc(:, params.crop1(1):params.crop1(2));
    othersdata = others.avg_all(:, params.crop2(1):params.crop2(2));
    
    % zscore
    data = zscore(data, [], 2);
    othersdata = zscore(othersdata, [], 2);
    
    % put in brain order
    data_map = NaN(length(keptData), size(data,2));
    data_map((keptData==1),:) = data;

    others_map = NaN(length(keptOthers), size(othersdata,2));
    others_map((keptOthers==1),:) = othersdata;
    
    % get shared voxels
    data_kept = data_map(keptVox.all,:);
    others_kept = others_map(keptVox.all,:);
    
    % correlate
    corr_data(i,:) = nansum(data_kept'.*others_kept')/(size(data_kept,2)-1);

end
fprintf('done! \n Making files ...');

% Get average correlation
isc = nanmean(corr_data,1)';

%% Save
% Make nifti map
map = NaN(length(keptVox.grp1),1);
map(keptVox.all) = isc;
map = reshape(map, 91, 109, 91);
standard_map = fullfile(params.projdir, params.group, 'data', 'MNI152_T1_2mm_brain.nii');
nii = make_niftiMap(map, standard_map, fullfile(iscParams.savedir, ...
    ['isc_' iscParams.name '_' iscParams.type '_ss' num2str(length(datafiles)) '.nii']));

% Save ISC data
save(fullfile(iscParams.savedir,['isc_' iscParams.name '_' iscParams.type '_ss' num2str(length(datafiles)) '_2.mat']), ...
    'corr_data', 'keptVox', 'iscParams');
fprintf('done! \n');

