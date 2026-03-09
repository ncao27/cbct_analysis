function [mask] = sphere(nx, ny, nz, dx, dy, dz, lesion_radius)
    lesion_cx = (nx / 2);   % lesion center in x-coordinates
    lesion_cy = (ny / 2);   % lesion center in y-coordinates
    lesion_cz = (nz / 2);   % lesion center in z-coordinates
    
    % convert lesion radius from mm to voxels
    lesion_rx = lesion_radius / dx;     % number of lesion voxels in the x-direction
    lesion_ry = lesion_radius / dy;     % number of lesion voxels in the y-direction
    lesion_rz = lesion_radius / dz;     % number of lesion voxels in the z-direction 
    
    % generate a mask for the lesion:
    [X, Y, Z] = ndgrid(1:nx, 1:ny, 1:nz);   % generate the mask grid first
    
    % define the mask for the sphere (notice it's just sphere equation)
    % also this makes mask into a logical mask of 0's and 1's
    mask = ((X - lesion_cx).^2 / lesion_rx^2 + ...
            (Y - lesion_cy).^2 / lesion_ry^2 + ...
            (Z - lesion_cz).^2 / lesion_rz^2) <= 1;
    mask = double(mask);
end