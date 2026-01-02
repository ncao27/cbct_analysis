indir = 'C:\Users\kli\Dropbox\Personal\R01_PCCT_Dropbox\Data\Head_Phantom\PHANTOM_10020202_10022020_KL\__20201002_095834_000000\DCT_HEAD_CLEAR_NAT_FILL_FULL_HU_NORMAL_[AX3D]_0009\';
out_filename = '.\projections.raw';
files = dir([indir '*.IMA']);
header = dicominfo([indir files(1).name]);
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
%%%%% User-defined projection angular range %%%%%
n_view = 20;
theta_array = linspace(0, 2*pi, n_view);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image = dicomread([indir files(100).name]);
% figure;imshow(image,[]);
vol = zeros(nx, ny, nz, 'single');
for z = 1:nz
    image = dicomread([indir files(z).name]);
    vol(:,:,z) = single(image);
end
vol = vol + header.RescaleIntercept;
vol = (vol + 1000)/1000*mu_water;        % this step is to convert the HU to mu. 

%%
fid = fopen(out_filename, 'w');
% -------------------------
% Voxel boundary planes
% -------------------------
x_plane = ((0:nx) - nx/2) * dx;
y_plane = ((0:ny) - ny/2) * dy;
z_plane = ((0:nz) - nz/2) * dz;

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
    imshow(prj',[0 5]);
    title_str = sprintf('Calculating view #%d/%d (angle = %.0f deg)\n',view, n_view, theta/pi*180);
    title(title_str); 
    set(h,'position',[200 200 700 600]); pause(0.5);
    fwrite(fid, prj, 'single');
end
fclose(fid);