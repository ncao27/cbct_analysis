% ==========================================================
% Fan-Beam Forward Projector (GPU)
% ==========================================================

% clear; clc;

% ---------------- Input / Output ----------------
in_filename = '.\without_lesion.raw';
out_filename = '.\projections.raw';


% ---------------- Volume geometry ----------------
fov = 400;
nx_s = single(512);
ny_s = single(512);
dx = single(fov/nx_s);
dy = dx;
% ---------------- CT geometry ----------------
sid = single(750);      % source to isocenter
sdd = single(1200);     % source to detector
nu  = 2480;            % detector elements
du  = single(0.35);  % detector spacing

% ---------------- Projection angles ----------------
n_view = 2000;
d_theta = 2*pi/n_view;
theta_array = 0:d_theta:2*pi-d_theta;

% ---------------- Load image ----------------
fid = fopen(in_filename, 'r', 'b');
vol = fread(fid, [nx_s, ny_s], 'single');
fclose(fid);
vol_gpu = gpuArray(vol);

figure; imagesc(vol', [0.015, 0.035]); colormap gray; axis image; colorbar;
title('Original Reconstruction');
%%
fid = fopen(out_filename,'w');

% ---------------- Detector fan angles ----------------
u = ((1:nu)' - (nu/2 + 0.5)) * du ;

% ---------------- Ray sampling along each detector ----------------
L_max = sid + sdd;              % approx max ray length
nstep =  ceil(L_max / (dx/2));   % step size <= 0.5*dx
t = gpuArray(((0:nstep-1) + 0.5) / nstep);  % midpoint rule

%% ==========================================================
% Forward projection
%% ==========================================================
for view = 1:n_view
    theta = theta_array(view);
    fprintf('View %d / %d\n', view, n_view);

    % -------- Source position --------
    xs = sid*cos(theta);
    ys = sid*sin(theta);

    % -------- Detector positions along circular arc --------
    xd = xs - sdd*cos(theta) + u*cos(theta + pi/2);
    yd = ys - sdd*sin(theta) + u*sin(theta + pi/2);

    % -------- Ray direction --------
    dxr = xd - xs;   % nu x 1
    dyr = yd - ys;

    % -------- Ray lengths & step size --------
    L = sqrt(dxr.^2 + dyr.^2);   % nu x 1
    ds = L / nstep;

    % -------- Vectorized sampling along each ray --------
    xt = xs + dxr * t;  % nu x nstep
    yt = ys + dyr * t;

    % -------- Image coordinates (indices) --------
    ix = xt ./ dx + nx_s/2 + 0.5; %ix = xt ./ dx + nx_s/2 + 1;
    iy = yt ./ dy + ny_s/2 + 0.5; %iy = yt ./ dy + ny_s/2 + 1;

    % flatten for single interp2 call
    ixv = ix(:);
    iyv = iy(:);

    % -------- Interpolation --------
    val = interp2(vol_gpu, ixv, iyv, 'cubic', 0);
    val = reshape(val, [nu, nstep]);

    % -------- Line integral along ray --------
    prj_gpu = sum(val, 2) .* ds;

    % -------- Write projection --------
    fwrite(fid, gather(prj_gpu'), 'single');
end

fclose(fid);
fprintf('Forward projection complete.\n');


