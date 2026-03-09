function [projections] = siddoncu(nx, ny, nz, dx, dy, dz, sid, sdd, du, dv, nu, nv, n_view, theta_array, vol)
    
    gpuDevice(1);
    mexcuda siddon_kernel.cu
    
    % todo: keep turning everything into a gpuarray on matlab, and then make
    % the necessary changes in cuda. also, this loop where i declare all the
    % variables, maybe i can pre compute all the variables and store into gpu
    % array and then every loop i can just call the gpuarray instead of having
    % to compute and create a separate instance of it, which prolly saves a lot
    % of memory
    
    x_plane = gpuArray(single(((0:nx) - nx/2) * dx)); % note that we go from -nx/2 to +nx/2
    y_plane = gpuArray(single(((0:ny) - ny/2) * dy));
    z_plane = gpuArray(single(((0:nz) - nz/2) * dz));
    projections_vol = zeros(nu, nv, n_view, 'single');
    
    % because CUDA converts it to float and vol used to be double
    vol = single(vol);
    gpu_vol = gpuArray(vol);
    nu = single(nu);
    nv = single(nv);
    
    % preallocate gpu arrays for the detector pixels
    u = gpuArray(single(1:nu)); % u direction on the detector plane
    v = gpuArray(single(1:nv)); % v direction on the detector plan
    
    for view = 1:n_view
    
        theta = theta_array(view);
        fprintf('Calculating view #%d/%d (angle = %.0f deg)\n',view, n_view, theta/pi*180);
        % -------------------------
        % Source position
        % -------------------------
        xs = single(sid * cos(theta));
        ys = single(sid * sin(theta));
        zs = single(0);
    
        % -------------------------
        % Detector center
        % -------------------------
        xd0 = xs - (sdd/sid) * xs;
        yd0 = ys - (sdd/sid) * ys;
        zd0 = 0;
    
        ux = cos(theta - pi/2); %  unit vector
        uy = sin(theta - pi/2);
    
        % we preallocate the x,y,z coordinates of each detector pixel, and then
        % loop through them using the cpp kernel
        xd = single(xd0 + (u - (nu/2+0.5)) * du * ux); 
        yd = single(yd0 + (u - (nu/2+0.5)) * du * uy);
        zd = single(zd0 + (v - (nv/2+0.5)) * dv);
        
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
        prj_vol = gpuArray.zeros(nu, nv, "single");
        
        % call cuda siddon kernel 
        siddon_kernel(gpu_vol, prj_vol, x_plane, y_plane, z_plane, ...
                 xs, ys, zs, xd, yd, zd, nu, nv, dx, dy, dz);
    
        % gather
        prj_vol = gather(prj_vol);
    
        % store projections
        projections_vol(:,:,view) = prj_vol;
    
        if mod(view,10)==0
            fprintf('Finished view %d/%d\n', view, n_view);
        end
    end
    projections = projections_vol;
end