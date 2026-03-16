%% data importation and declaring necessary varibles
clear;
reset(gpuDevice(1))             % clear lingering threads and kill extra processes
indir = 'C:\Users\Nathan Cao\OneDrive\Desktop\ct images analysis\cbct_head_phantom\2d_cbct_imgs\';
imgs = dir([indir '*.raw']);
lesion = dir([indir '*.jpg']);
nx = 512;                       % number of pixels in x-dir (given)
ny = 512;                       % number of pixels in y-dir (given)
nz = 1;
dx = 400/512;           % x-dir pixel size (fov given)
dy = 400/512;           % y-dir pixel size (fov given)
dz = dx;
sid = 750;                      % Source-to-ISO Distance. unit: mm
sdd = 1200;                     % Source-to-Detector Distance. unit: mm
du = 0.154*4;                   % 2D detector pixel size. unit: mm 
dv = dz;                        % 2D detector pixel size. unit: mm
nu = 2480/4;                    % number of detector pixels along u (native: 1920 × 2480) 
nv = 1;                         % number of detector pixels along v
mu_water = 0.02;                % unit: 1/mm

% we want 20 views, so we get 20 angles from 0 to 360
n_view = 360; 
theta_array = linspace(0, 2*pi, n_view);

% store the with lesion image
fid = fopen([indir imgs(1).name],'r','ieee-be');                % get fid of data (raw), 'ieee-be' reads big-endian style data
img_with_lesion = fread(fid, [nx ny],'single');                 % read data, 32-bit representation, so single precision
fclose(fid);                                                    % close file
img_with_lesion = img_with_lesion';                             % transpose
%img_with_lesion = (img_with_lesion + 1000) / 1000 * mu_water;   % no need for conversion, it's already good 

% store the without lesion image
fid = fopen([indir imgs(2).name],'r','ieee-be');                    % get fid of data (raw), 'ieee-be' reads big-endian style data
img_without_lesion = fread(fid, [nx ny],'single');                  % read data, 32-bit representation, so single precision
fclose(fid);                                                        % close file
img_without_lesion = img_without_lesion';                           % transpose
%img_without_lesion = (img_without_lesion + 1000) / 1000 * mu_water; % no need for conversion, it's already good 

% store the segmented lesion image
lesion = imread([indir lesion(1).name]);

%% Image Without Lesion: Forward Projection
projections_img = squeeze(siddon.siddoncpp_two(nx, ny, dx, dy, sid, sdd, du, nu, n_view, theta_array, img_without_lesion));

%% Segmented Lesion: Forward Projection
% cast the lesion to type double
lesion = double(lesion);

% drop the white background 
lesion(lesion == 150) = 0;

projections_lesion = squeeze(siddon.siddoncpp_two(nx, ny, dx, dy, sid, sdd, du, nu, n_view, theta_array, lesion));

%% 

%% FBP reconstruction
reconstructed_img = recon.fbp(nx, ny, dx, dy, sid, sdd, du, nu, n_view, theta_array, projections_img);
%%
viz.double(img_without_lesion)

