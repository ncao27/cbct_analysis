function [postRegistered, D, info] = demons_iodine_registration(preImage, postImage, varargin)
%DEMONS_IODINE_REGISTRATION Register post-iodine CT onto the pre-iodine grid.
%
%   [POSTREGISTERED,D,INFO] = DEMONS_IODINE_REGISTRATION(PRE,POST) accepts
%   grayscale 2-D slices or 3-D intensity volumes. PRE is fixed; POST moves.
%   This standalone function uses Image Processing Toolbox's IMREGDEMONS
%   (multi-resolution registration in the Thirion demons algorithm family).
%   It requires MATLAB with Image Processing Toolbox; a GPU is not required.
%
%   Where to use it:
%     Reconstruct/load each scan -> put both on the SAME spatial grid ->
%     rigidly align if needed -> call this function -> inspect registration ->
%     subtract the pre-injection image from the registered post-injection image.
%   Inputs must already share size, orientation, spacing, and field of view.
%   Array size alone cannot establish spatial agreement. DICOM metadata is not
%   inferred and neither scanner parameters nor reconstruction code is changed.
%   Large rotations/translations should be handled before local registration.
%
%   Example with reconstructed arrays (HU or attenuation units are retained):
%     [postOnPre,D,info] = demons_iodine_registration(preImage,postImage);
%     enhancement = postOnPre - double(preImage);
%     enhancement(~info.ValidMask) = NaN;
%
%   Run without inputs for a reproducible synthetic-motion/iodine demo:
%     [postOnPre,D,info] = demons_iodine_registration();
%   The demo prints numerical checks and saves one comparison figure under
%   images/demons_registration_demo/ beside this file. Normal calls save no files.
%
%   Name-value options (registration settings, NOT CT geometry parameters):
%     'Iterations'         [100 50 25]: iterations from coarse to full scale.
%                          One entry per pyramid level; a scalar uses one level.
%     'FieldSmoothing'     1.5: accumulated-field Gaussian sigma in voxels.
%     'RegistrationMode'   'local-normalized' (default) or 'intensity'.
%     'NormalizationSigma' 5: local mean/std Gaussian sigma in voxels.
%     'NoiseFloor'         0.03: std floor relative to the fixed feature range.
%     'RegistrationPair'   {} or {fixedFeatures,movingFeatures}, each on its
%                          respective input's grid and with the same array size.
%                          E.g. supply iodine-insensitive structural images:
%       [postOnPre,D,info] = demons_iodine_registration(preImage,postImage,...
%           'RegistrationPair',{preStructural,postStructural},...
%           'RegistrationMode','intensity');
%
%   Iodine changes intensity without necessarily changing anatomy. The default
%   mode estimates motion on locally standardized COPIES to reduce sensitivity
%   to local brightness/contrast differences. It does not make demons invariant
%   to iodine: newly enhanced structures can still bias the field. The optional
%   feature pair lets you drive registration using shared anatomical structure.
%   No histogram matching, normalization, or intensity correction is applied to
%   the returned post-injection image: the ORIGINAL POST is linearly resampled.
%   Interpolation itself can change voxel values and blur small enhanced regions.
%
%   Outputs:
%     POSTREGISTERED: double, same size as PRE, original POST intensity units.
%                     Samples outside the moving image's center grid are NaN.
%     D: backward sampling displacement, in PIXELS/VOXELS, not millimeters.
%        D(:,:,1)=dx (columns), D(:,:,2)=dy (rows) for 2-D;
%        D(:,:,:,1:3) contains [dx,dy,dz] for 3-D. At fixed coordinate (x,y,z),
%        POSTREGISTERED samples POST at (x+dx,y+dy,z+dz). Do not negate D or
%        treat it as a forward transform of moving landmarks.
%     INFO.ValidMask: samples with full linear-interpolation support in POST.
%     INFO.FeatureRMSEBefore/After: feature-space errors on that SAME mask;
%        these describe the optimized features, not iodine signal or truth.
%     INFO.MinJacobian / NonpositiveJacobianFraction: finite-difference checks
%        of the sampling map x+D(x). Inspect nonpositive values before use.
%     INFO.Options / UsedRegistrationPair / ValidFraction: run settings/support.
%
%   All smoothing acts in pixel/voxel units. Anisotropic volumes therefore have
%   anisotropic physical smoothing; choose a common grid deliberately upstream.
%   No patient-image accuracy can be established by the synthetic demo alone.
%
%   References:
%     Thirion, Medical Image Analysis 2(3), 243-260 (1998).
%     https://www.mathworks.com/help/images/ref/imregdemons.html
%     https://www.mathworks.com/help/images/ref/imwarp.html

if nargin == 0
    [postRegistered, D, info] = syntheticDemo();
    return;
end
narginchk(2, inf);
validateImage(preImage, 'preImage');
validateImage(postImage, 'postImage');
if ~isequal(size(preImage), size(postImage))
    error('demons_iodine_registration:GridMismatch', ...
        'PRE and POST must already be resampled onto the same spatial grid and array size.');
end
if exist('imregdemons', 'file') ~= 2
    error('demons_iodine_registration:MissingToolbox', ...
        'This function requires imregdemons from Image Processing Toolbox.');
end

parser = inputParser;
parser.FunctionName = mfilename;
parser.PartialMatching = false;
positiveScalar = @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x > 0;
iterationsOK = @(x) isnumeric(x) && isreal(x) && isvector(x) && ...
    ~isempty(x) && all(isfinite(x)) && all(x >= 1 & x == fix(x));
addParameter(parser, 'Iterations', [100 50 25], iterationsOK);
addParameter(parser, 'FieldSmoothing', 1.5, positiveScalar);
addParameter(parser, 'RegistrationMode', 'local-normalized', ...
    @(x) (ischar(x) && isrow(x)) || (isstring(x) && isscalar(x)));
addParameter(parser, 'NormalizationSigma', 5, positiveScalar);
addParameter(parser, 'NoiseFloor', 0.03, positiveScalar);
addParameter(parser, 'RegistrationPair', {}, @(x) iscell(x) && (isempty(x) || numel(x) == 2));
parse(parser, varargin{:});
options = parser.Results;
options.Iterations = double(options.Iterations(:)');
options.FieldSmoothing = double(options.FieldSmoothing);
options.NormalizationSigma = double(options.NormalizationSigma);
options.NoiseFloor = double(options.NoiseFloor);
options.RegistrationMode = validatestring(options.RegistrationMode, ...
    {'local-normalized', 'intensity'}, mfilename, 'RegistrationMode');
levels = numel(options.Iterations);
if any(size(preImage) < 2^levels)
    error('demons_iodine_registration:TooManyLevels', ...
        'Every image dimension must be at least 2^numel(Iterations). Use fewer pyramid levels.');
end

usedPair = ~isempty(options.RegistrationPair);
if usedPair
    fixed = options.RegistrationPair{1};
    moving = options.RegistrationPair{2};
    validateImage(fixed, 'fixedFeatures');
    validateImage(moving, 'movingFeatures');
    if ~isequal(size(fixed), size(preImage)) || ~isequal(size(moving), size(postImage))
        error('demons_iodine_registration:FeatureSizeMismatch', ...
            'Each registration feature image must match its corresponding input grid.');
    end
else
    fixed = preImage;
    moving = postImage;
end
% Remove large feature arrays from the settings returned in INFO.
options = rmfield(options, 'RegistrationPair');
[fixed, moving] = registrationFeatures(fixed, moving, options);

% Estimate the deformation on registration copies only. The second output of
% imregdemons would be a warped FEATURE image, so do not use it as POSTREGISTERED.
if isequal(fixed, moving)
    D = zeros([size(preImage), ndims(preImage)], 'double');
else
    D = imregdemons(moving, fixed, options.Iterations, ...
        'PyramidLevels', levels, ...
        'AccumulatedFieldSmoothing', options.FieldSmoothing, ...
        'DisplayWaitbar', false);
end
if isempty(D) || any(~isfinite(D(:)))
    error('demons_iodine_registration:InvalidField', ...
        'Demons did not return a finite displacement field. Inspect the input features.');
end

valid = validSamplingMask(D, size(preImage));
% Casting integer CT data to double before interpolation avoids integer rounding.
postRegistered = imwarp(double(postImage), D, 'linear', ...
    'FillValues', 0, 'SmoothEdges', false);
postRegistered(~valid) = NaN;

info = struct;
if nargout > 2
    info.ValidMask = valid;
    info.ValidFraction = nnz(valid) / numel(valid);
    movingAfter = imwarp(moving, D, 'linear', 'FillValues', 0, 'SmoothEdges', false);
    if any(valid(:))
        before = double(moving(valid)) - double(fixed(valid));
        after = double(movingAfter(valid)) - double(fixed(valid));
        info.FeatureRMSEBefore = sqrt(mean(before.^2));
        info.FeatureRMSEAfter = sqrt(mean(after.^2));
    else
        info.FeatureRMSEBefore = NaN;
        info.FeatureRMSEAfter = NaN;
    end
    [info.MinJacobian, info.NonpositiveJacobianFraction] = jacobianSummary(D);
    info.Options = options;
    info.UsedRegistrationPair = usedPair;
end
end

function validateImage(image, name)
validateattributes(image, {'numeric', 'logical'}, ...
    {'real', 'finite', 'nonempty', 'nonsparse'}, mfilename, name);
if isa(image, 'gpuArray') || ndims(image) > 3 || any(size(image) < 2)
    error('demons_iodine_registration:InvalidImage', ...
        '%s must be a CPU grayscale 2-D image or 3-D volume, at least two samples per dimension.', name);
end
end

function [fixed, moving] = registrationFeatures(fixed, moving, options)
% Apply the same affine intensity scale to both registration copies. This
% preserves differences between tissues; it does not histogram-match iodine.
% Local standardization/clipping, when enabled, is a separate feature step.
fixed = double(fixed);
moving = double(moving);
low = min(fixed(:));
span = max(fixed(:)) - low;
if span <= 0 || max(moving(:)) == min(moving(:)) || ~isfinite(span)
    error('demons_iodine_registration:ConstantFeatures', ...
        'Registration images must contain finite spatial intensity variation.');
end
fixed = single((fixed - low) / span);
moving = single((moving - low) / span);
if any(~isfinite(fixed(:))) || any(~isfinite(moving(:)))
    error('demons_iodine_registration:FeatureRange', ...
        'Feature intensity ranges cannot be represented reliably in single precision.');
end
if strcmp(options.RegistrationMode, 'local-normalized')
    fixed = localStandardize(fixed, options.NormalizationSigma, options.NoiseFloor);
    moving = localStandardize(moving, options.NormalizationSigma, options.NoiseFloor);
end
if any(~isfinite(fixed(:))) || any(~isfinite(moving(:))) || ...
        max(fixed(:)) == min(fixed(:)) || max(moving(:)) == min(moving(:))
    error('demons_iodine_registration:FeatureRange', ...
        'Feature preprocessing lost finite spatial contrast. Check image ranges and normalization settings.');
end
end

function normalized = localStandardize(image, sigma, noiseFloor)
if ismatrix(image)
    meanImage = imgaussfilt(image, sigma, 'Padding', 'replicate');
    secondMoment = imgaussfilt(image.^2, sigma, 'Padding', 'replicate');
else
    meanImage = imgaussfilt3(image, sigma, 'Padding', 'replicate');
    secondMoment = imgaussfilt3(image.^2, sigma, 'Padding', 'replicate');
end
variance = max(secondMoment - meanImage.^2, 0);
normalized = (image - meanImage) ./ sqrt(variance + noiseFloor^2);
% Limit the influence of extreme local changes in the registration copy only.
normalized = (min(max(normalized, -3), 3) + 3) / 6;
end

function valid = validSamplingMask(D, imageSize)
% Array dimensions are [row column slice]; field components are [x y z].
order = [2 1 3];
dimensions = numel(imageSize);
valid = true(imageSize);
for component = 1:dimensions
    dimension = order(component);
    coordinateShape = ones(1, dimensions);
    coordinateShape(dimension) = imageSize(dimension);
    coordinate = reshape(1:imageSize(dimension), coordinateShape);
    fieldIndex = repmat({':'}, 1, dimensions + 1);
    fieldIndex{end} = component;
    query = coordinate + D(fieldIndex{:});
    valid = valid & query >= 1 & query <= imageSize(dimension);
end
end

function [minimum, nonpositiveFraction] = jacobianSummary(D)
% Compute det(I + grad(D)) a slice at a time to limit 3-D working memory.
minimum = inf;
nonpositive = 0;
if ndims(D) == 3
    [ux, uy] = gradient(D(:,:,1));
    [vx, vy] = gradient(D(:,:,2));
    determinant = (1 + ux).*(1 + vy) - uy.*vx;
    minimum = min(determinant(:));
    nonpositiveFraction = nnz(determinant <= 0) / numel(determinant);
    return;
end
slices = size(D, 3);
for slice = 1:slices
    [ux, uy] = gradient(D(:,:,slice,1));
    [vx, vy] = gradient(D(:,:,slice,2));
    [wx, wy] = gradient(D(:,:,slice,3));
    lower = max(1, slice - 1);
    upper = min(slices, slice + 1);
    uz = (D(:,:,upper,1) - D(:,:,lower,1)) / (upper - lower);
    vz = (D(:,:,upper,2) - D(:,:,lower,2)) / (upper - lower);
    wz = (D(:,:,upper,3) - D(:,:,lower,3)) / (upper - lower);
    determinant = (1 + ux).*((1 + vy).*(1 + wz) - vz.*wy) ...
        - uy.*(vx.*(1 + wz) - vz.*wx) + uz.*(vx.*wy - (1 + vy).*wx);
    minimum = min(minimum, min(determinant(:)));
    nonpositive = nonpositive + nnz(determinant <= 0);
end
nonpositiveFraction = nonpositive / (size(D,1)*size(D,2)*slices);
end

function [postRegistered, D, info] = syntheticDemo()
% Toy values below describe a synthetic test, not scanner or patient parameters.
[x, y] = meshgrid(1:128, 1:128);
pre = 0.10 + 0.65*exp(-((x-40)/18).^2 - ((y-43)/22).^2) ...
    + 0.45*exp(-((x-87)/16).^2 - ((y-82)/18).^2) ...
    + 0.35*exp(-((x-86)/10).^2 - ((y-35)/12).^2) ...
    + 0.30*exp(-((x-42)/10).^2 - ((y-94)/9).^2);
iodine = 0.35*exp(-((x-87)/10).^2 - ((y-82)/12).^2);
trueDx = 3 + 1.2*sin(2*pi*y/128);
trueDy = -2;
% Fixed-to-moving map is [x+trueDx(y), y+trueDy]. Its inverse is explicit.
sourceY = y - trueDy;
sourceX = x - (3 + 1.2*sin(2*pi*sourceY/128));
post = interp2(pre + iodine, sourceX, sourceY, 'linear', 0.10);
[postRegistered, D, info] = demons_iodine_registration(pre, post);

structuralROI = pre > 0.18 & iodine < 0.005 & info.ValidMask;
errorBefore = post(structuralROI) - pre(structuralROI);
errorAfter = postRegistered(structuralROI) - pre(structuralROI);
info.Demo.StructuralRMSEBefore = sqrt(mean(errorBefore.^2));
info.Demo.StructuralRMSEAfter = sqrt(mean(errorAfter.^2));
endpointError = hypot(D(:,:,1) - trueDx, D(:,:,2) - trueDy);
info.Demo.MedianEndpointErrorPixels = median(endpointError(structuralROI));
info.Demo.ExpectedEnhancement = iodine;
info.Demo.RecoveredEnhancement = postRegistered - pre;
assert(info.Demo.StructuralRMSEAfter < info.Demo.StructuralRMSEBefore, ...
    'demons_iodine_registration:DemoRegistration', 'Synthetic anatomical alignment did not improve.');
assert(info.NonpositiveJacobianFraction == 0, ...
    'demons_iodine_registration:DemoFolding', 'Synthetic deformation contains nonpositive Jacobians.');
fprintf('Synthetic structural RMSE: %.6g before -> %.6g after\n', ...
    info.Demo.StructuralRMSEBefore, info.Demo.StructuralRMSEAfter);
fprintf('Synthetic median displacement error: %.4f pixels\n', info.Demo.MedianEndpointErrorPixels);
fprintf('Minimum Jacobian: %.4f; nonpositive fraction: %.4g\n', ...
    info.MinJacobian, info.NonpositiveJacobianFraction);

fig = figure('Visible', 'off', 'Color', 'white', 'Position', [100 100 1320 760]);
cleanup = onCleanup(@() close(fig));
tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
intensityLimits = [min(pre(:)), max(post(:))];
panels = {pre, post, postRegistered, post-pre, postRegistered-pre, iodine};
labels = {'Pre-injection (fixed)', 'Post-injection + motion (moving)', ...
    'Registered post-injection', 'Post - pre, before registration', ...
    'Registered post - pre', 'Known enhancement (synthetic)'};
for index = 1:6
    ax = nexttile;
    panel = panels{index};
    if index <= 3
        artist = imagesc(ax, panel, intensityLimits);
        colormap(ax, gray(256));
    else
        artist = imagesc(ax, panel, [-0.15 0.4]);
        colormap(ax, parula(256));
    end
    set(artist, 'AlphaData', isfinite(panel));
    axis(ax, 'image'); axis(ax, 'off'); colorbar(ax);
    title(ax, labels{index});
end
sgtitle(sprintf('Synthetic demons test | structural RMSE %.4f -> %.4f | median field error %.2f pixels', ...
    info.Demo.StructuralRMSEBefore, info.Demo.StructuralRMSEAfter, info.Demo.MedianEndpointErrorPixels));
outdir = fullfile(fileparts(mfilename('fullpath')), 'images', 'demons_registration_demo');
if ~isfolder(outdir), mkdir(outdir); end
info.Demo.FigurePath = fullfile(outdir, 'synthetic_registration.png');
exportgraphics(fig, info.Demo.FigurePath, 'Resolution', 150);
fprintf('Synthetic comparison figure: %s\n', info.Demo.FigurePath);
end
