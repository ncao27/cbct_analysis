function [projections] = siddoncpp_two(nx, ny, dx, dy, sid, sdd, du, nu, n_view, theta_array, vol)
    mex siddon_kernel_two.cpp
    
    x_plane = double(((0:nx) - nx/2) * dx); % note that we go from -nx/2 to +nx/2
    y_plane = double(((0:ny) - ny/2) * dy);
    projections = zeros(nu, n_view, 'double');
    
    for view = 1:n_view
    
        theta = theta_array(view);
        fprintf('Calculating view #%d/%d (angle = %.0f deg)\n',view, n_view, theta/pi*180);
        % -------------------------
        % Source position
        % -------------------------
        xs = sid * cos(theta);
        ys = sid * sin(theta);
    
        % -------------------------
        % Detector center
        % -------------------------
        xd0 = xs - (sdd/sid) * xs;
        yd0 = ys - (sdd/sid) * ys;
    
        ux = cos(theta - pi/2); %  unit vector
        uy = sin(theta - pi/2);
    
        u = 1:nu; % u direction on the detector plane
    
        % we preallocate the x,y,z coordinates of each detector pixel, and then
        % loop through them using the cpp kernel
        xd = xd0 + (u - (nu/2+0.5)) * du * ux; 
        yd = yd0 + (u - (nu/2+0.5)) * du * uy;
        
        % C++ kernel doing siddon reconstruction
        % Inputs:
        %   - vol: the CT volume that we read in
        %   - x_plane: the x-dir planes that we define for voxels
        %   - y_plane: the y-dir planes that we define for voxels
        %   - z_plane: the z-dir planes that we define for voxels
        %   - xs: the x-coord of point source
        %   - ys: the y-coord of point source
        %   - zs: the z-coord of point source
        %   - xd: array of x-coordinates for detector plane
        %   - yd: array of y-coordinates for detector plane
        %   - zd: array of z-coordinates for detector plane
        %   - nu: the number of pixels in u direction on detector
        %   - nv: the number of pixels in v direction on detector
    
        % find projections of the volume
        prj_vol = siddon_kernel_two(vol, x_plane, y_plane, ...
                 xs, ys, xd, yd, nu, dx, dy);
        
        % store projections
        projections(:,view) = prj_vol;
    end
end