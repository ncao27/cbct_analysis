function [recon] = fdk(nx, ny, nz, dx, dy, dz, sid, sdd, du, dv, nu, nv, n_view, theta_array, projections)
    % call the fdk cuda kernel
    mexcuda fdk_kernel.cu
    
    % set the reconstruction array as a gpu array
    recon_gpu = gpuArray.zeros(nx,ny,nz,'single');
    
    freq = (-nu/2:nu/2-1)'/(nu*du);
    ramp = abs(2*pi*freq);
    
    projections = single(projections);
    
    % loop through the views as a 
    for view = 1:n_view
        
        % get the angle associated with this view
        theta = theta_array(view);
        
        % weighting but with meshgrid so faster
        [uu,vv] = meshgrid( (1:nu)-nu/2, (1:nv)-nv/2 );
        w = sid ./ sqrt(sid^2 + (uu'*du).^2 + (vv'*dv).^2);
        projections(:,:,view) = projections(:,:,view) .* w;
    
        % do ramp filtering
        for v = 1:nv
            
            % first do fft to frequency domain, then center projections on zero
            P = fftshift(fft(projections(:,v,view)));
            P = P .* ramp;
    
            % ifft the projects back
            projections(:,v,view) = real(ifft(ifftshift(P)));
        end
    
        % get the projection from a specific view
        proj_gpu  = gpuArray(projections(:,:,view));
    
        % run fdk on that view, we don't define plhs in cpp/cuda, so don't assign output
        % proj_gpu: the specific 2d gpu array corresponding to the view
        % recon_gpu: the 3d gpu array of the reconstructed volume
        % nx: the number of pixels in x direction of reconstruction
        % ny: the number of pixels in y direction of reconstruction
        % nz: the number of pixels in z direction of reconstruction
        % nu: the number of pixels in u direction of detector panel (local coordinates)
        % nv: the number of pixels in v direction of detector panel (local coordinates)
        % dx: the spacing of pixels in x direction of reconstruction
        % dy: the spacing of pixels in y direction of reconstruction
        % dz: the spacing of pixels in z direction of reconstruction
        % du: the spacing of pixels in u direction of detector panel
        % dv: the spacing of pixels in v direction of detector panel
        % sid: source to middle of volume distance
        % sdd: source to detector distance
        % theta: angle corresponding to the specific view we need to rotate source and detector by
        fdk_kernel(proj_gpu, recon_gpu, nx, ny, nz, nu, nv, dx, dy, dz, du, dv, sid, sdd, theta);
        fprintf('Finished view %d/%d\n', view, n_view);
    end
    
    % bring reconstructed array back to the gpu
    recon = gather(recon_gpu);
end