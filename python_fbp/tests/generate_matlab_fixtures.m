% Generate small CPU oracles from untouched li_fbp source scripts.
% Run from the repository root. Outputs go to a temporary directory.
% Only test dimensions/pitch and GPU/display plumbing change in memory.
root = getenv('CBCT_REPOSITORY');
if isempty(root)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
out = getenv('CBCT_FIXTURE_OUTPUT');
if isempty(out), out = fullfile(tempdir, 'cbct_python_validation'); end
if ~isfolder(out), mkdir(out); end
n = 12; detectors = 16; views = 9; pitch = 35;
vol = reshape(0.02 + 0.007*sin(1:n*n), n, n);
fid = fopen(fullfile(out, 'fixture_input.raw'), 'w', 'b');
fwrite(fid, vol, 'single'); fclose(fid);
for geometry = {'linear', 'curved'}
    detector = geometry{1};
    source = fileread(fullfile(root, 'li_fbp', ['c1_forward_projection_v6_' detector '_detector.m']));
    source = strrep(source, '.\without_lesion.raw', fullfile(out, 'fixture_input.raw'));
    source = strrep(source, '.\projections.raw', fullfile(out, [detector '_projections.raw']));
    source = strrep(source, 'single(512)', 'single(12)');
    source = strrep(source, 'nu  = 2480;', 'nu  = 16;');
    source = strrep(source, 'single(0.35)', 'single(35)');
    source = strrep(source, 'n_view = 2000;', 'n_view = 9;');
    source = strrep(source, 'gpuArray(', '(');
    source = strrep(source, 'gather(', '(');
    source = regexprep(source, '(?m)^figure;[^\r\n]*', '');
    eval(source);
    fprintf('Forward interpolation output type: %s\n', class(val));
    source = fileread(fullfile(root, 'li_fbp', ['c2_fbp_v6_' detector '_detector.m']));
    source = strrep(source, 'gpuDevice(1);', '');
    source = strrep(source, 'detector_nums = 2480;', 'detector_nums = 16;');
    source = strrep(source, 'views = 2000;', 'views = 9;');
    source = strrep(source, 'du = 0.35;', 'du = 35;');
    source = strrep(source, 'image_size = 512;', 'image_size = 12;');
    source = strrep(source, '''projections.raw''', ['''' fullfile(out, [detector '_projections.raw']) '''']);
    source = strrep(source, '''img.raw''', ['''' fullfile(out, [detector '_img.raw']) '''']);
    source = strrep(source, 'gpuArray.zeros(', 'zeros(');
    source = strrep(source, ', ''gpuArray''', '');
    source = strrep(source, 'gpuArray(', '(');
    source = strrep(source, 'gather(', '(');
    source = regexprep(source, '(?m)^(figure;|imagesc\([^\r\n]*|colormap gray;[^\r\n]*|title\([^\r\n]*)', '');
    eval(source);
    fid = fopen(fullfile(out, [detector '_filtered.raw']), 'w', 'l');
    fwrite(fid, filtered_data, 'single'); fclose(fid);
    fid = fopen(fullfile(out, [detector '_kernel.raw']), 'w', 'l');
    fwrite(fid, ramp_kernel, 'single'); fclose(fid);
end
v = reshape(sin(1:30), 5, 6);
xq = repmat([0.9 1 1.1 1.7 2 2.4 3.2 4.9 5 5.8 6 6.1], 7, 1);
yq = repmat([0.9;1;1.2;2.6;4.7;5;5.1], 1, 12);
z = interp2(v, xq, yq, 'cubic', 0);
for item = {'v', 'xq', 'yq', 'z'}
    key = item{1};
    fid = fopen(fullfile(out, ['interp_' key '.raw']), 'w', 'l');
    fwrite(fid, eval(key), 'double'); fclose(fid);
end
fprintf('CPU MATLAB fixtures written to %s\n', out);

% Production-size view checks retain every scanner/sampling parameter.
% Only select four entries of the original 2000-view angular grid.
for geometry = {'linear', 'curved'}
    detector = geometry{1};
    source = fileread(fullfile(root, 'li_fbp', ['c1_forward_projection_v6_' detector '_detector.m']));
    source = strrep(source, '.\without_lesion.raw', fullfile(root, 'li_fbp', 'without_lesion.raw'));
    source = strrep(source, '.\projections.raw', fullfile(out, [detector '_full_selected.raw']));
    source = strrep(source, 'for view = 1:n_view', 'for view = [1 251 501 751]');
    source = strrep(source, 'gpuArray(', '(');
    source = strrep(source, 'gather(', '(');
    source = regexprep(source, '(?m)^figure;[^\r\n]*', '');
    eval(source);
end
