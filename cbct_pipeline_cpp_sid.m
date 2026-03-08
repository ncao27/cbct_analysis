%% data importation and declaring necessary varibles
clear;
indir = 'C:\Users\Nathan Cao\OneDrive\Desktop\ct images analysis\cbct_head_phantom\DCT_HEAD_CLEAR_NAT_FILL_FULL_HU_NORMAL_[AX3D]_0009\';
out_filename = '.\projections.raw';
files = dir([indir '*.IMA']);
header = dicominfo([indir files(1).name]); % getting the dicom scan information
nx = single(header.Rows);       % 3D image matrix size along x
ny = single(header.Columns);    % 3D image matrix size along y
nz = length(files);             %y 3D image matrix size along z
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

%% Define lesion
% Notes to self: so lesion insertion above is actually wrong
% Create 3D lesion volume -> project the lesion -> 
% Do some kind of filtering: Gaussian blurring (basic, but it helps soften the edges around the 
% manually created lesion), Poisson editing (preserves the intensity changes of the inserted lesion)
% MTF filtering (good because it takes the CT scanner's measured MTF and probably convolves it with
% the image), NPS filtering (using the measured noise of the CT scanner, since a lot of noise is not 
% random white noise but structured noise, you do some transfer function / filtering using the
% measured noise)

% define parameters
lesion_radius = 10;     % (mm)
lesion_mu = 0.03;       % higher than soft tissue which is typically 0.02
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

% <=1 makes it logical, convert it to double
mask = double(mask);

%% Siddon + Projection Domain Lesion Insertion with C++ Kernel

mex siddon_kernel.cpp

x_plane = double(((0:nx) - nx/2) * dx); % note that we go from -nx/2 to +nx/2
y_plane = double(((0:ny) - ny/2) * dy);
z_plane = double(((0:nz) - nz/2) * dz);
projections_vol = zeros(nu, nv, n_view, 'double');
projections_les = zeros(nu, nv, n_view, 'double');

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

    u = 1:nu; % u direction on the detector plane
    v = 1:nv; % v direction on the detector plan

    % we preallocate the x,y,z coordinates of each detector pixel, and then
    % loop through them using the cpp kernel
    xd = xd0 + (u - (nu/2+0.5)) * du * ux; 
    yd = yd0 + (u - (nu/2+0.5)) * du * uy;
    zd = zd0 + (v - (nv/2+0.5)) * dv;
    
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
    prj_vol = siddon_kernel(vol, x_plane, y_plane, z_plane, ...
             xs, ys, zs, xd, yd, zd, nu, nv, dx, dy, dz);
    
    % find projections of the lesion
    prj_les = siddon_kernel(mask, x_plane, y_plane, z_plane, ...
             xs, ys, zs, xd, yd, zd, nu, nv, dx, dy, dz);
    
    % store projections
    projections_vol(:,:,view) = prj_vol;
    projections_les(:,:,view) = prj_les;

end

%% Projections with lesion
prj_w_lesion = projections_vol + projections_les;

%% fdk reconstruction with CUDA kernel

% call the fdk cuda kernel
mexcuda fdk_kernel.cu

% set the reconstruction array as a gpu array
recon_gpu = gpuArray.zeros(nx,ny,nz,'single');

freq = (-nu/2:nu/2-1)'/(nu*du);
ramp = abs(2*pi*freq);

projections = single(prj_w_lesion);

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

%% Bring reconstructed volume back to CPU
recon = gather(recon_gpu);   % bring back to CPU

%% Scroll through visualization for reconstructed volume
recon = squeeze(recon);

figure;
sliceViewer(recon);
colormap gray;

%% 3D visualization for reconstructed volume
figure;
volshow(recon);