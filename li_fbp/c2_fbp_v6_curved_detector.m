%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPU-Accelerated Fan-Beam FBP (MATHEMATICALLY IDENTICAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear; clc;

gpuDevice(1); % initialize GPU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Parameters and Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detector_nums = 2480;
views = 2000;

SID = 750;
SDD = 1200;
du = 0.35;

d_theta = 2*pi/views;
theta = 0:d_theta:2*pi-d_theta;

u = ((1:detector_nums) - (detector_nums/2+0.5)) * du;
d_gamma = du/sdd;    % detector angular spacing
gamma = ((1:detector_nums) - (detector_nums/2+0.5)) * d_gamma; % radians

fov = 400;
image_size = 512;
dx = fov/image_size;

% --- Read projection data ---
file = fopen('projections.raw', 'r', 'l');
projection_data = fread(file, [detector_nums views], 'single');
fclose(file);

projection_data = projection_data';
% projection_data = fliplr(projection_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move to GPU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

projection_data = gpuArray(single(projection_data));
gamma = gpuArray(single(gamma));
theta = gpuArray(single(theta));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Ramp Filtering (EXACT SAME AS CPU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- build ramp kernel with sin(gamma) as argument ---
ramp_kernel = zeros(1, 2*detector_nums-1, 'single', 'gpuArray');

for k = 1:2:detector_nums-1
    val = -1/pi^2/(sin(k*d_gamma))^2; %-1/(k*pi*du)^2;
    ramp_kernel(detector_nums-k) = val;
    ramp_kernel(detector_nums+k) = val;
end
ramp_kernel(detector_nums) = 1/(2*d_gamma)^2;

% --- weighting ---
projection_weighted = cos(gamma) .* projection_data;

% --- convolution ---
filtered_data = zeros(size(projection_data), 'single', 'gpuArray');

ramp_matrix = repmat(ramp_kernel, views, 1);

for col = 1:detector_nums
    filtered_data(:, col) = sum( ...
        projection_weighted .* ...
        ramp_matrix(:, detector_nums-col+1 : 2*detector_nums-col), ...
        2) * d_gamma;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Backprojection (GPU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = image_size;

x1d = ((1:M)-M/2-0.5)*dx;
x1d = gpuArray(single(x1d));

[x, y] = meshgrid(x1d);

img = gpuArray.zeros(M, M, 'single');

for idx_angle = 1:views
    
    if mod(idx_angle,50)==0
        fprintf('View %d / %d\n', idx_angle, views);
    end
    
    ct = cos(theta(idx_angle));
    st = sin(theta(idx_angle));
    
    L = sqrt((x-SID*ct).^2 + (y-SID*st).^2);    %U = (SID - x*ct - y*st);
    gamma0 = asin((x*st-y*ct)./L);     %u0 = SDD*(-x*st + y*ct) ./ U;
    
    % detector index
    k = floor((gamma0 - gamma(1))/d_gamma) + 1;
    
    % valid mask (GPU-friendly)
    valid = (k >= 1) & (k < detector_nums);
    
    k(~valid) = 1; % safe indexing
    
    w = gamma0 - gamma(k);
    w(~valid) = 0;
    
    % interpolation
    val = (1-w).*reshape(filtered_data(idx_angle, k), M, M) + ...
          w.*reshape(filtered_data(idx_angle, k+1), M, M);
    
    img = img + SID ./ (L.^2) .* val * d_theta;
end

img = img / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Move back to CPU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

img = gather(img);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = fopen('img.raw','w');
fwrite(file, img, 'float');
fclose(file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
imagesc(img', [0.015, 0.035]);
colormap gray; axis image; colorbar;
title('FBP Reconstruction (GPU, Exact)');

