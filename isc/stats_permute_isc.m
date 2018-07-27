function corr_data = stats_permute_isc(params, subs, conds, permute)
% Run permutations for assessing significance of ISC

subCond = params.scans{conds(1)};
othersCond = params.scans{conds(2)};

fprintf(['\n *** Calc phase-scrambled ISC (' params.isc.type '-group): ' subCond ' v ' othersCond '***\n']);


%% If between, load average of others or sum of others as appropriate
fprintf('Loading average of others \n');

if strcmp(params.isc.type, 'between')
    %Load average of others
    avg_others_file = fullfile(params.datadir2, ['avg_all_' othersCond '_n18']);
    others_data = load(avg_others_file);
    avg_others = others_data.avg_all;
    
else % load sum of all
    if ~exist(params.sumfile, 'file')
        fprintf('sum all file does not exist, creating now ...\n')
        for i = 1:length(subs)
            fprintf([subs{i} '\n']);
            sub_data_file = fullfile(params.datadir1, [subs{i} '_' subCond '.mat']);
            sub_data = load(sub_data_file);
            if i == 1
                sum_all = zscore(sub_data.tc, 0, 2);
            else
                sum_all = sum_all + zscore(sub_data.tc, 0, 2);
            end
        end
        save(params.sumfile, 'sum_all');
    end
    sum_all_data = load(params.sumfile);
end

%% If between, load masks
if strcmp(params.isc.type, 'between')
    mask_group1 = load(params.maskfile1);
    mask_group2 = load(params.maskfile2);
    
    % new mask is intersection of masks
    keptVox1 = find(mask_group1.keptVox==1);
    keptVox2 = find(mask_group1.keptVox==1);
    mask_both = intersect(keptVox1, keptVox2);
end

%% Do stats
for iter = 1:params.iterations
    tic
    fprintf(['\n*** Iter = ' num2str(iter) ' ***\n']);
    savename = [params.savename num2str(iter) '.mat'];
    if exist(savename, 'file')
        fprintf('file already exists, skipping\n')
        continue;
    end
    
    if strcmp(params.isc.type, 'within')
        corr_data_sum = zeros(1,length(sum_all_data.sum_all));
    else
        corr_data_sum = zeros(1,length(mask_both));
    end
    
    for i = 1:length(subs)
        
        %Load data -----------------
        fprintf(['loading ' subs{i} '...']);
        
        % Load subdata
        sub_data_file = fullfile(params.datadir1, [subs{i} '_' subCond '.mat']);
        sub_data = load(sub_data_file);
        
        % If within subject, get average of others
        if strcmp(params.isc.type, 'within')
            avg_others = (sum_all_data.sum_all - zscore(sub_data.tc,[], 2)) ...
                / (length(subs)-1);           
        end
        
        % Other preprocessing -------------
        % Crop
        fprintf('processing...')
        if ~isempty(params.crop)
            sub_data = sub_data.tc(:, params.crop(1):params.crop(2));
            others = avg_others(:, params.crop(1):params.crop(2));
        else
            sub_data = sub_data.tc;
            others = avg_others;
        end
        
        % Zscore
        others = zscore(others,0, 2);
        sub_data = zscore(sub_data,0,2);
        
        % Mask data
        if strcmp(params.isc.type, 'between')
            sub_data = maskData(sub_data, mask_group1, mask_both);
            others = maskData(others, mask_group2, mask_both);
        end
        
        % Phase scram & corr -------------
        fprintf('phase scram...corr...');
        
        others_scram = NaN(size(others));
        chunks = round(length(others)/5000);
        for j = 1:chunks
            % which vox?
            vox_start = (j-1)*5000+1;
            vox_end = j*5000;
            if vox_end > length(others_scram)
                vox_end = length(others_scram);
            end
            
            % phase scramble others
            others_scram(vox_start:vox_end,:) = phase_rand(others(vox_start:vox_end,:)', permute)';
        end
        
        % calc correlation subject to phase scrambled avg others
        corr_data = sum(sub_data'.*others_scram')/(size(sub_data,2)-1);
        corr_data_sum = corr_data_sum+corr_data;
        
        fprintf('done!\n')
    end
    
    % Calc boots ----------------------
    bootstrap = corr_data_sum / length(subs);
    
    save(savename, 'bootstrap', 'params');
    sprintf('\tsaved %s%d', params.savename, iter);
    toc
end
end

%% Helper functions

function data_masked = maskData(data, oldMask, newMask)
    totalVox = length(oldMask.keptVox);
    nTRs = size(data,2);
    
    % unmask data
    brain = NaN(totalVox, nTRs);
    brain(oldMask.keptVox==1,:) = data;

    % remask data
    data_masked = brain(newMask,:);
end









