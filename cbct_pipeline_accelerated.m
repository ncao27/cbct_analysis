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

% simulated lesion parameters (spherical)
lesion_radius = 10;     % (mm)
lesion_mu = 0.03;       % higher than soft tissue which is typically 0.02

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

%% Lesion creation

mask = lesion.sphere(nx, ny, nz, dx, dy, dz, lesion_radius);

%% CUDA Siddon

projections_vol = siddon.siddoncu(nx, ny, nz, dx, dy, dz, sid, sdd, du, dv, nu, nv, n_view, theta_array, vol);

%% C++ Siddon

projections_vol = siddon.siddoncpp(nx, ny, nz, dx, dy, dz, sid, sdd, du, dv, nu, nv, n_view, theta_array, vol);

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
%%
one = squeeze(projections_vol(:, :, 1));
imshow(one)