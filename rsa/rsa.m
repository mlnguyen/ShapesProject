function rsa(params, task_id)
% function rsa(params, task_id)
% runs intersubject RSA. Params are set in wrapper function

%% Load cubes
startTime = tic;

fprintf('\n ***** Loading cubes ***** \n');
[~,name,~] = fileparts(params.cubefile);
fprintf(['Cube file: ' name '\n']);
cubes = load(params.cubefile);

toc(startTime);

%% Load LSA
tic;
fprintf('\n ***** Loading LSA data ***** \n');
lsa = load(params.lsafile);

% get comparison matrix
lsa_mat = lsa.(params.lsa_grp)(params.keepsubs1, params.keepsubs2);

% get average comparison lsa similarity if doing avgOthers comparison
if strcmp(params.rsatype, 'avgOthers')
    fprintf('Averaging lsa similarity...\n')
    if strcmp(params.isctype, 'within')
        lsa_new = NaN(length(params.keepsubs1),1);
        for i = 1:length(params.keepsubs1)
            others = setdiff([1:length(params.keepsubs1)], i);
            lsa_new(i) = mean(lsa_mat(i,others));
        end
        lsa_mat = lsa_new;
    elseif strcmp(params.isctype, 'between')
        lsa_mat = mean(lsa_mat,2);
    else
        fprintf('invalid rsa type')
        return
    end
end

toc;

%% Load data

tic;
fprintf('\n ***** Loading neuro data ***** \n');

% Inds of good cubes
cube_inds = find(cubes.keptCubes);

% Which cubes to load?
which_cubes = ((task_id-1)*params.nCubes + 1) : (task_id-1)*params.nCubes + params.nCubes;
if which_cubes(end) > length(cube_inds)
    which_cubes = ((task_id-1)*params.nCubes + 1) : length(cube_inds);
end

% What are the inds of the voxels in these cubes?
vox = cell2mat(cubes.cubes(cube_inds(which_cubes)));
vox_inds = unique(sub2ind([91 109 91], vox(:,1), vox(:,2), vox(:,3)));

% Get the subdata in these voxels
vox_data1 = load_subdata(params.datafiles1, vox_inds);
if strcmp(params.isctype, 'between')
    vox_data2 = load_subdata(params.datafiles2, vox_inds);
else
    vox_data2 = vox_data1;
end

% get cube numbers 
cubeNums = cube_inds(which_cubes);

toc

%% Do searchlight
tic;
fprintf('\n ***** LSA searchlight ***** \n');

curTime = 0;
r_pear = NaN(length(which_cubes),1);
r_spear = NaN(length(which_cubes),1);
kept_cubes = 0;
if strcmp(params.rsatype, 'pairwise')
    isc_mat_all = NaN(length(which_cubes), length(params.datafiles1), length(params.datafiles2));
elseif strcmp(params.rsatype, 'avgOthers')
     isc_mat_all = NaN(length(which_cubes), length(params.datafiles1));
end

for j = 1:length(cubeNums)
    
    if mod(j, 100) == 0
        fprintf(['cube ' num2str(j) ' (t= ' num2str(toc-curTime) ')\n']);
        curTime = toc;
    end;
    
    % get a cube
    this_cube = cubes.cubes{cubeNums(j)};
    cube_inds = sub2ind([91 109 91], this_cube(:,1), this_cube(:,2), this_cube(:,3));
    data_inds = find(ismember(vox_inds, cube_inds));
    
    % get data in cube
    cubedata1 = vox_data1(:, data_inds,:);
    cubedata2 = vox_data2(:, data_inds,:);
    
    % If more than half of voxels are nans, then skip
    nans1 = sum(isnan(cubedata1(1,:,1)));
    nans2 = sum(isnan(cubedata2(1,:,1)));
    if nans1 > size(cubedata1,2)/2 || nans2 > size(cubedata2,2)/2
        continue;
    end
    kept_cubes = kept_cubes+1;
    
    % average across voxels
    mean_cube1 = squeeze(nanmean(cubedata1,2));
    mean_cube2 = squeeze(nanmean(cubedata2,2));
 
    % crop & zscore
    mean_cube1 = mean_cube1(:,params.crop(1):params.crop(2));
    mean_cube2 = mean_cube2(:,params.crop(1):params.crop(2));
    
    % Get neural similarity: pairwise
    if strcmp(params.rsatype, 'pairwise')
        % Find isc cormat
        isc_mat = corr(mean_cube1(params.keepsubs1,:)', mean_cube2(params.keepsubs2,:)');
        isc_mat_all(j,:,:) = isc_mat;
        
        % Which matrix cells to compare?
        if strcmp(params.isctype, 'within')
            keep = tril(ones(length(lsa_mat), length(lsa_mat)), -1);
        else
            keep = ones(numel(lsa_mat),1);
        end
    
    % Get neural similarity: avg others
    elseif strcmp(params.rsatype, 'avgOthers')
        isc_mat = NaN(length(params.keepsubs1),1);
        
        % if within, compare to the average of others within group
        if strcmp(params.isctype, 'within')
            for sub=1:length(params.keepsubs1)
                subdata = mean_cube1(sub,:);
                avg_others = mean(mean_cube2(setdiff(params.keepsubs2,sub),:));
                isc_mat(sub) = corr(subdata', avg_others');
            end
        
        % if between, compare to the average of second group   
        elseif strcmp(params.isctype, 'between')
            avg_others = mean(mean_cube2,1);
            isc_mat = corr(mean_cube1', avg_others');
        end
        
        isc_mat_all(j,:) = isc_mat;
        keep = ones(size(isc_mat));
        
    else
        fprintf('error: invalid rsa type\n')
        return
    end
    
    % Correlate with behavior
    r_pear(j) = corr(isc_mat(keep==1), lsa_mat(keep==1), 'rows', 'complete', 'type', 'pearson');
    r_spear(j) = corr(isc_mat(keep==1), lsa_mat(keep==1), 'rows', 'complete', 'type', 'spearman');
end

cube_centers = cubes.cube_centers(cubeNums,:);
fprintf(['\nCalculated RSA for ' num2str(kept_cubes) ' of ' num2str(length(cubeNums)) ' cubes\n'])
toc

%% Save
tic; fprintf('\n ***** Saving searchlight ***** \n');

savename = fullfile(params.outputdir, ['searchlt_rsa_' params.cond '_' params.isctype ...
   '_' params.rsatype '_ss' num2str(length(params.datafiles1)) '_' num2str(task_id)]);

save(savename, 'r_pear', 'r_spear', 'cube_centers', 'isc_mat_all', 'lsa_mat', 'keep');

fprintf('done! \n');
toc
end

%% Helper functions

function vox_data = load_subdata(datafiles, vox_inds)
for i = 1:length(datafiles)
    
    str = sprintf('loading %s ...', datafiles{i});
    disp(str);
    
    % load datafile
    data = load(datafiles{i});
    
    % move to brain map
    data_map = NaN(size(data.keptVox,1), size(data.tc,2));
    data_map(data.keptVox>0,:) = data.tc;
    
    % get vox data
    vox_data(i,:,:) = data_map(vox_inds,:);
    
end
end




