%% siddon algorithm
fid = fopen(out_filename, 'w');

% -------------------------
% Voxel boundary planes
% -------------------------
x_plane = ((0:nx) - nx/2) * dx; % note that we go from -nx/2 to +nx/2
y_plane = ((0:ny) - ny/2) * dy;
z_plane = ((0:nz) - nz/2) * dz;

x_plane_gpu = gpuArray(x_plane);
y_plane_gpu = gpuArray(y_plane);
z_plane_gpu = gpuArray(z_plane);

projections = gpuArray(zeros(nu, nv, n_view, 'single'));

h = figure; 
parfor view = 1:n_view
    theta = theta_array(view);
    fprintf('Calculating view #%d/%d (angle = %.0f deg)\n',view, n_view, theta/pi*180);
    % -------------------------
    % Source position
    % -------------------------
    xs = sid * cos(theta);
    ys = sid * sin(theta);
    zs = 0;

    % -------------------------
    % Detector center
    % -------------------------
    xd0 = xs - (sdd/sid) * xs;
    yd0 = ys - (sdd/sid) * ys;
    zd0 = 0;

    ux = cos(theta - pi/2); %  unit vector
    uy = sin(theta - pi/2);

    prj = zeros(nu, nv);
    
    % can try parallelizing this part, can use a meshgrid 
    for u = 1:nu
        xd = xd0 + (u - (nu/2+0.5)) * du * ux;
        yd = yd0 + (u - (nu/2+0.5)) * du * uy;

        for v = 1:nv
            zd = zd0 + (v - (nv/2+0.5)) * dv;

            dxr = xd - xs;
            dyr = yd - ys;  
            dzr = zd - zs;

            L = sqrt(dxr^2 + dyr^2 + dzr^2);

            % Siddon parameters
            ax = (x_plane - xs) / dxr;
            ay = (y_plane - ys) / dyr;
            az = (z_plane - zs) / dzr;

            a_min = max([min(ax), min(ay), min(az)]);
            a_max = min([max(ax), max(ay), max(az)]);

            if a_min >= a_max
                continue;
            end

            % Collect intersections
            a = [ax, ay, az];
            a = a(a >= a_min & a <= a_max);
            a = sort(a);

            % Integrate voxel contributions
            for n = 1:length(a)-1
                amid = 0.5*(a(n)+a(n+1));

                xmid = xs + amid*dxr;
                ymid = ys + amid*dyr;
                zmid = zs + amid*dzr;

                i = floor((xmid - x_plane(1))/dx) + 1;
                j = floor((ymid - y_plane(1))/dy) + 1;
                k = floor((zmid - z_plane(1))/dz) + 1;

                if i>=1 && i<=nx && j>=1 && j<=ny && k>=1 && k<=nz
                    prj(u,v) = prj(u,v) + ...
                        vol(i,j,k) * (a(n+1)-a(n)) * L;
                end
            end
        end
    end
    projections(:,:,view) = prj;
end






    



