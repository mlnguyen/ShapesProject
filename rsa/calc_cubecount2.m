
function calc_cubecount2(mask_img,gap,rad,savename)
% Create cubes for conducting searchlight analysis

% Reshape and binarize brain mask (nVox x 1)
mask = reshape(mask_img,[],1);
mask(mask>0) = 1; 

% Define a grid with within the brain map and the given gap. Each vertex
% (x,y,z coordinate) is the center of a cube to be created. nCubes x 3
% (847,547 x)
cube_centers = find_grid_vertices(size(mask_img,1),size(mask_img,2),size(mask_img,3),gap);
cube_inds = sub2ind([91 109 91], cube_centers(:,1), cube_centers(:,2), cube_centers(:,3));

% Make cubes
cubes = cell(length(cube_centers),1);
keptCubes = zeros(length(cube_centers),1);

fprintf('**** Making cubes ... \n'); tic
for i = 1:size(cube_centers,1)
    
    if mod(i,10000) == 0
        fprintf(['cube ' num2str(i) '\n'])
    end
    
    % if the cube center is inside the brian, then get cube
    if mask(cube_inds(i)) == 1
        this_cube = grow_cube(cube_centers(i,1),cube_centers(i,2),cube_centers(i,3),rad);
        
        % Once have cube, check that entirety of cube is inside image
        try
            inds = sub2ind([91 109 91], this_cube(:,1), this_cube(:,2), this_cube(:,3));
        catch
            continue;
        end
        
        % check that at least half the cube is in the brain
        nVoxInBrain = sum(mask(inds));
        if nVoxInBrain > length(this_cube)/2
            cubes{i} = this_cube;
            keptCubes(i) = 1;
        end
    end
end

fprintf('*** done! \n');
toc;

% save
save(savename, 'keptCubes','cubes','cube_centers', '-v7.3');