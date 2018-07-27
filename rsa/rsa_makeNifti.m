% rsa_makeNifti.m
% Following rsa searchlight and permute, create niftis and get significance

datadir = '';
savedir = fullfile(datadir,'niftis');
mnidir = '';
codedir = '';

addpath(genpath(codedir));

filename = 'shapesMovie_smooth4mm2_motionResid';
%% Get data

% load base nifti
nifti = load_nii(fullfile(mnidir, 'MNI152_T1_2mm_brain.nii'));

% which datafiles?
datafiles = dir(fullfile(datadir, 'searchlt*'));

% initialize matrices
r_spearman = NaN(numel(nifti.img),1);
p_spearman = NaN(numel(nifti.img),1);

% loop through all datafiles
for i = 1:length(datafiles)
    
    % load data
    fprintf(['datafile ' num2str(i) '...\n']);
    data = load(fullfile(datadir, datafiles(i).name));
    
    % get center of cubes
    cube_inds = sub2ind([91 109 91], data.cube_centers(:,1), data.cube_centers(:,2), data.cube_centers(:,3));
    
    % get correlation vals
    r_spearman(cube_inds) = data.r_spear; 
 
    % get stats
    mu = mean(data.stats.null_spearman,2); 
    sigma = std(data.stats.null_spearman,[],2); 
    p_spearman(cube_inds) = normcdf(data.r_spear, mu, sigma, 'upper'); 
    p_boot(cube_inds) = data.stats.p_spearman;
    
end

% FDR correct
[pcrit4, sigval, siginds4] = fdr_BH(p_spearman(~isnan(p_spearman)), .05);
fprintf(['spearman critical val (norm): ' num2str(pcrit4) '\n']);

%% Make niftis

if ~exist(savedir, 'dir'), mkdir(savedir); end;

% set nifti params
nifti.hdr.dime.datatype = 64;
nifti.hdr.dime.bitpix = 32;
nifti.hdr.dime.cal_min=0.15;
nifti.hdr.dime.cal_max=0.3;

% spearman r_values
nifti.img = reshape(r_spearman, size(nifti.img,1), size(nifti.img,2), size(nifti.img,3));
save_nii(nifti, fullfile(savedir, [filename '_r_spearman4.nii']));

% spearman p values
nifti.img = reshape(p_spearman, size(nifti.img,1), size(nifti.img,2), size(nifti.img,3));
save_nii(nifti, fullfile(savedir, [filename '_p_spearman4_all.nii']));

% spearman -log(p) values
nifti.img = reshape(-1*log(p_spearman), size(nifti.img,1), size(nifti.img,2), size(nifti.img,3));
save_nii(nifti, fullfile(savedir, [filename '_logp_spearman4_all.nii']));



