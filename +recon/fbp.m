function [recon] = fbp(nx, ny, dx, dy, sid, sdd, du, nu, n_view, theta_array, projections)
    disp('Starting Filtered Back Projection...');
    
    % we create a (x,y) 2D grid for the whole image with ndgrid
    xc = ((1:nx) - nx/2 - 0.5) * dx; 
    yc = ((1:ny) - ny/2 - 0.5) * dy;
    [X, Y] = ndgrid(xc, yc); % Create a 2D grid for the whole image

    % because we are doing fan beam (with one source) reconstruction, weigh the different detector pixels
    u_idx = (1:nu) - (nu/2 + 0.5); 
    u_dist = u_idx * du;
    pre_weight = sdd ./ sqrt(sdd^2 + u_dist.^2)'; 
    weighted_prj = projections .* pre_weight; 

    % ramp filtering in the frequency domain
    f = (-nu/2:nu/2-1)'/(nu*du);
    H = abs(2*pi*f);

    % Apply FFT, multiply by Ramp filter, apply IFFT, and remove padding
    fft_prj = fftshift(fft(weighted_prj, [], 1));
    filtered_fft_prj = fft_prj .* H;
    filtered_prj = real(ifft(ifftshift(filtered_fft_prj), [], 1));

    % now we do fbp but weighted by distance / angle
    recon = zeros(nx, ny, 'double');
    d_theta = (2*pi) / n_view; 

    for view = 1:n_view
        theta = theta_array(view);
        
        % Source position
        xs = sid * cos(theta);
        ys = sid * sin(theta);
        
        % Compute distances for the entire X, Y grid simultaneously
        % U: Distance from source to pixel projected onto the central ray
        U = sid - X.*cos(theta) - Y.*sin(theta);
        
        % V: Distance from central ray to the pixel (parallel to detector)
        V = X.*sin(theta) - Y.*cos(theta);
        
        % Find where this pixel projects onto the detector (u_det)
        u_det = V .* (sdd ./ U);
        
        % Convert physical detector distance back to a bin index
        u_bin = (u_det / du) + (nu/2 + 0.5);
        
        % Create a mask so we don't interpolate outside the detector
        valid_mask = (u_bin >= 1) & (u_bin <= nu);
        
        % Interpolate the filtered projection at the exact bin locations
        % (We use 'linear' interpolation for smooth edges)
        proj_interpolated = interp1(1:nu, filtered_prj(:, view), u_bin(valid_mask), 'linear');
        
        % Add to the reconstruction grid, applying the 1/U^2 distance weight
        recon(valid_mask) = recon(valid_mask) + (proj_interpolated ./ (U(valid_mask).^2)) * d_theta;
    end

    % Scale the final image based on the geometry to recover proper attenuation values
    recon = recon * (sid^2 / sdd);

end

%filtered_prj_pad = real(ifft(fft(weighted_prj, N_pad, 1) .* H, [], 1));
%filtered_prj = filtered_prj_pad(1:nu, :);

%%
% ramp filtering in the frequency domain
%N_pad = 2^nextpow2(nu * 2); 
%f = [0:(N_pad/2-1), -(N_pad/2):-1] / N_pad; % Frequency axis
%H = abs(f)'; % Ram-Lak (Ramp) Filter

% Apply FFT, multiply by Ramp filter, apply IFFT, and remove padding
%fft_prj = fftshift(fft(weighted_prj, N_pad, 1));
%filtered_fft_prj = fft_prj .* H;
%filtered_prj = real(ifft(ifftshift(filtered_fft_prj), [], 1));
%filtered_prj = filtered_prj(1:nu, :);
