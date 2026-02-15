clear;

%% data importation and declaring necessary variables
clear;
indir = 'C:\Users\Nathan Cao\OneDrive\Desktop\ct images analysis\cbct_head_phantom\DCT_HEAD_CLEAR_NAT_FILL_FULL_HU_NORMAL_[AX3D]_0009\';
out_filename = '.\projections.raw';
files = dir([indir '*.IMA']);
header = dicominfo([indir files(1).name]); % getting the dicom scan information
nx = single(header.Rows);       % 3D image matrix size along x
ny = single(header.Columns);    % 3D image matrix size along y
nz = length(files);             % 3D image matrix size along z
dx = header.PixelSpacing(1);    % 3D image pixel size. unit: mm
dy = header.PixelSpacing(2);    % 3D image pixel size. unit: mm
dz = header.SliceThickness;     % 3D image pixel size. unit: mm
sid = 750;                      % Source-to-ISO Distance. unit: mm
sdd = 1200;                     % Source-to-Detector Distance. unit: mm
du = 0.154*4;                   % 2D detector pixel size. unit: mm 
dv = 0.154*4;                   % 2D detector pixel size. unit: mm
nu = 2480/4;                    % number of detector pixels along u (native: 1920 × 2480)
nv = 1920/4;                    % number of detector pixels along v (native: 1920 × 2480)
mu_water = 0.02;                % unit: 1/mm

% we want 20 views, so we get 20 angles from 0 to 360
n_view = 360; 
theta_array = linspace(0, 2*pi, n_view);

% define the volume that we store the CT image in 
vol = zeros(nx, ny, nz, 'double');
for z = 1:nz
    image = dicomread([indir files(z).name]);
    vol(:,:,z) = double(image);
end
vol = vol + header.RescaleIntercept;
vol = (vol + 1000) / 1000 * mu_water;        % this step is to convert the HU to mu. 

%% siddon algorithm
fid = fopen(out_filename, 'w');

% -------------------------
% Voxel boundary planes
% -------------------------
x_plane = ((0:nx) - nx/2) * dx; % note that we go from -nx/2 to +nx/2
y_plane = ((0:ny) - ny/2) * dy;
z_plane = ((0:nz) - nz/2) * dz;
projections = gpuArray(zeros(nu, nv, n_view, 'double'));



h = figure; 
for view = 1:n_view
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
    %imshow(prj',[0 5]);
    %title_str = sprintf('Calculating view #%d/%d (angle = %.0f deg)\n',view, n_view, theta/pi*180);
    %title(title_str); 
    %set(h,'position',[200 200 700 600]); pause(0.5);
    %fwrite(fid, prj, 'double');
    disp("one done")
end
fclose(fid);

%% visualize slices of the forward projection
imagesc(projections(:, :, 120)')
axis image

%% fdk reconstruction of the ct volume

% define the reconstructed volume
recon = zeros(nx, ny, nz, 'single');

% define the ramp filter
% idea here: we use a ramp filter because fbp intrinsically emphasizes low
% frequency components while blurring out higher frequencies, which makes
% something like a ramp filter ideal since it gives edges more definition
freq = (-nu/2:nu/2-1)'/(nu*du);
ramp = abs(freq);

% go through the different projects and do fdk
for view = 1:n_view
    theta = theta_array(view);

    % weigh the reconstruction by distance, go through pixel by pixel since
    % they all have different distance weights
    for u = 1:nu
        for v = 1:nv

            % define the weight
            w = sid / sqrt(sid^2 + ((u-nu/2)*du)^2 + ((v-nv/2)*dv)^2);

            % take the projects and weight it by their distances
            projections(u,v,view) = projections(u,v,view)*w;
        end
    end

    % do ramp filtering
    for v = 1:nv
        
        % first do fft to frequency domain, then center projections on zero
        P = fftshift(fft(projections(:,v,view)));
        P = P .* ramp;

        % ifft the projects back
        projections(:,v,view) = real(ifft(ifftshift(P)));
    end

    % do fdk
    for ix = 1:nx
        x = (ix-nx/2)*dx;
        for iy = 1:ny
            y = (iy-ny/2)*dy;
            for iz = 1:nz
                z = (iz-nz/2)*dz;

                xs = sid*cos(theta);
                ys = sid*sin(theta);

                denom = sid - x*cos(theta) - y*sin(theta);
                if denom <= 0, continue; end

                u = (sdd/denom)*( -x*sin(theta) + y*cos(theta) )/du + nu/2;
                v = (sdd/denom)*z/dv + nv/2;

                iu = round(u);      
                iv = round(v);

                if (iu >= 1) && (iu <= nu) && (iv >= 1) && (iv <= nv)
                    recon(ix,iy,iz) = recon(ix,iy,iz) + projections(iu,iv,view)*(sid^2/denom^2);
                end
            end
        end
    end
    disp(view)
end
recon = recon * (2*pi/n_view);

%% for algorithmic testing purposes: do fbp reconstruction of one slice
% so currently the projects matrix is of size 620 by 480 by 360, idk 
% why exactly we did this, but basically 620 is the number of columns
% and 480 is the number of rows, nu is horizontal along the detector 
% axis and nv is verticle along the detector axis
% so first we can take one particular row (ie one slice) of the detector
% projections and all 360 projections, so that would give us a matrix 
% of 620 by 360

slice = projections(:, ceil(nv / 2), :);    % get middle slice
slice = squeeze(slice);                     % get rid of third dim

recon = zeros(nx, ny, 'single');
freq = (-nu/2:nu/2-1)'/(nu*du);
ramp = abs(freq);

for view = 1:n_view

    % get the angle that detector panel is at
    theta = theta_array(view);

    for u = 1:nu
        % define weight
        w = sid / sqrt(sid^2 + ((u-nu/2)*du)^2);
        
        % weigh voxels
        slice(u, view) = slice(u, view)*w;
    end

    % do ramp filtering
    P = fftshift(fft(slice(:,view)));
    P = P .* ramp;
    slice(:,view) = real(ifft(ifftshift(P)));

    % do fbp reconstruction
    for ix = 1:nx
        x = (ix - nx/2)*dx;
        for iy = 1:ny
            y = (iy - ny/2)*dy;

            xs = sid*cos(theta);
            xy = sid*cos(theta);

            denom = sid - x*cos(theta) - y*sin(theta);
            if denom <= 0, continue; end

            u = (sdd/denom)*( -x*sin(theta) + y*cos(theta) )/du + nu/2;

            iu = round(u);

            if (iu >= 1) && (iu <= nu)
                recon(ix, iy) = recon(ix, iy) + slice(iu, view)*(sid^2/denom^2);
            end
        end
    end

    disp(view)
end
recon = recon * (2*pi / n_view);
%%
imagesc(recon)
axis image