
function calc_iscWholeBrainWithin(params, iscParams, datafiles, keptVox)
% Calculate ISC within a group of subjects (each subject to average of
% others). Saves results both as nifti and .mat file
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


fprintf(['\n *** Calculating ISC (' iscParams.type '-group): ' iscParams.name '***\n']);

%% Calc sum of all
fprintf('calculating sum of all...\n')
for i = 1:length(datafiles)
    fprintf([num2str(i) '-']);

    % load subdata
    data = load(datafiles{i});
    
    if i == 1
        sum_all = data.tc;
    else
        sum_all = sum_all + data.tc;
    end
end    
savename = fullfile(iscParams.savedir, ['sum_all_' iscParams.name]);
save(savename, 'sum_all');

%% Do ISC
fprintf('\ndoing isc...\n');
for i = 1:length(datafiles)

    fprintf([num2str(i) '-']);
    
    % load subdata
    data = load(datafiles{i});
    
   % get avg others
    others = (sum_all - zscore(data.tc, [], 2)) / (length(params.subs)-1);
    
    % crop and zscore
    data = zscore(data.tc(:, params.crop(1):params.crop(2)),[],2);
    others = zscore(others(:, params.crop(1):params.crop(2)), [],2);
        
    % correlate
    corr_data(i,:) = nansum(data'.*others')/(size(data,2)-1);

end

fprintf('done! \n Making files ...');

% Get average correlation
isc = nanmean(corr_data,1)';

%% Save
% Make nifti map
map = NaN(length(keptVox),1);
map(keptVox == 1) = isc;
map = reshape(map, 91, 109, 91);
standard_map = fullfile(params.projdir, 'data', 'MNI152_T1_2mm_brain.nii');
nii = make_niftiMap(map, standard_map, [iscParams.savename '.nii']);

% Save ISC data
save(iscParams.savename, 'corr_data', 'keptVox', 'iscParams');
fprintf('done! \n');
