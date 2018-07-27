
% Load data
datadir = 'Z:\mai\projects\shapesStory\between\analysis\iscWholeBrain\1x2';
data1 = load(fullfile(datadir, 'isc_mov2tom_smooth4mm2_motionAudResid_half1_between_ss18_2'));
data2 = load(fullfile(datadir, 'isc_mov2tom_smooth4mm2_motionAudResid_half2_between_ss18_2'));
standard_map = fullfile('Z:\mai\projects\shapesStory\fmri_group1\data', 'MNI152_T1_2mm_brain.nii');


% calculate difference
[h,p,ci,stats] = ttest2(data1.corr_data, data2.corr_data);

% unmask ttest results and diff_isc
results_unmask = NaN(91*91*109, 1);
results_unmask(data1.keptVox.all) = stats.tstat*-1;
results_unmask_p = NaN(91*91*109, 1);
results_unmask_p(data1.keptVox.all) = p;

% load mask
%mask = load_nii('Z:\mai\projects\shapesStory\fmri_group1\analysis\iscWholeBrain\smooth4mm\mask_mov2mov_n36.nii');
mask = load_nii(fullfile(datadir, 'mask_mov2tom_n36.nii'));
mask_flat = reshape(mask.img, [91*109*91,1]);
results_unmask(isnan(mask_flat)) = NaN;

results_unmask_p(isnan(mask_flat)) = NaN;
nii = make_niftiMap(reshape(results_unmask, [91, 109,91,1]), standard_map, ...
    fullfile(datadir, 'half_ttest_mask_mov2tom_motionAudResid_tstat.nii'));