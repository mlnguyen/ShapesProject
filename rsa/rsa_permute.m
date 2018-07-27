
function rsa_searchlt_permute(params, task_id)

% Load RSA data
data = load(params.datafiles{task_id});


% Allocate space for permutations
r_pearson = NaN(length(data.cube_centers), params.iters);
r_spearman = NaN(length(data.cube_centers), params.iters);
p_pearson = NaN(length(data.cube_centers),1);
p_spearman = NaN(length(data.cube_centers),1);

% Run permutations
tic;
for cube = 1:length(data.cube_centers)
    
    if mod(cube,10) == 0
        fprintf(['cube ' num2str(cube) '...']);
        toc
    end
    
    if isnan(data.r_spear(cube))
        continue;
    end
    
    % Get isc matrix
    isc_mat = squeeze(data.isc_mat_all(cube,:,:));
    
    for iter = 1:params.iters
        
        % Shuffle rows and cols of isc matrix
        if strcmp(params.rsatype, 'pairwise')
            if strcmp(params.isctype, 'within')
                shuff = randperm(length(isc_mat));
                shuffled_isc = isc_mat(shuff, shuff);
            else
                shuff1 = randperm(size(isc_mat,1));
                shuff2 = randperm(size(isc_mat,2));
                shuffled_isc = isc_mat(shuff1, shuff2);
            end
        elseif strcmp(params.rsatype, 'avgOthers')
            shuff = randperm(length(isc_mat));
            shuffled_isc = isc_mat(shuff)';
        else
            fprintf('error: invalid rsa type')
            return;
        end
            
        % get cells to compare
        isc_keep = shuffled_isc(data.keep==1);
        lsa_keep = data.lsa_mat(data.keep==1);
        
        % Correlate
        r_pearson(cube, iter) = corr(isc_keep, lsa_keep, 'type', 'pearson');
        r_spearman(cube,iter) = corr(isc_keep, lsa_keep, 'type', 'spearman');
    end
    
    % calculate p-values
    p_pearson(cube) = 1-length(r_pearson(cube, r_pearson(cube,:) < data.r_pear(cube)))/params.iters;
    p_spearman(cube) = 1-length(r_pearson(cube, r_spearman(cube,:) < data.r_spear(cube)))/params.iters;

end

% save
savename = fullfile(params.outputdir, ['searchlt_rsa_' params.isctype '_'...
    params.rsatype '_ss36_' num2str(task_id)]);
data.stats.null_pearson = r_pearson;
data.stats.null_spearman = r_spearman;
data.stats.p_pearson = p_pearson;
data.stats.p_spearman = p_spearman;

save(savename, '-struct', 'data');


