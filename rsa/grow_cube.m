
function cube = grow_cube(x,y,z,r)
% accepts xyz coordinates and r, and returns indices for a cube of
% radius r centered on xyz.

xcoords = [x-r:x+r]';
ycoords = [y-r:y+r]';
zcoords = [z-r:z+r]';

sidel = length(xcoords);

cube = []; cubex = []; cubey = []; cubez = [];
for xi = xcoords'
    cubex = [cubex; repmat(xi,sidel^2,1)];
end
for yi = ycoords'
    cubey = [cubey; repmat(yi,sidel,1)];
end
cubey = repmat(cubey,sidel,1);
cubez = repmat(zcoords,sidel^2,1);
    
cube = [cubex cubey cubez];

end