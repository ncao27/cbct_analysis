function validate_lesion_insertion(stage, cache)
% Reproduce one paired insertion experiment, preserving the li_fbp sources.
% Run 'compute' first, then one export stage at a time:
% comparison, difference, arrays, pixels, metrics.
% Requires the preceding lesion-free result and MATLAB Image Processing Toolbox.
if nargin < 2, cache=fullfile(tempdir,'cbct_lesion_insertion_validation'); end
out=fileparts(mfilename('fullpath'));
root=fileparts(fileparts(out));
if ~exist(cache,'dir'), mkdir(cache); end
result=fullfile(cache,'computed.mat');
if strcmp(stage,'compute')
    prior=fullfile(root,'images','lesion_free_roundtrip_case_b_del_slice071');
    base=load(fullfile(prior,'03_roundtrip_arrays.mat'));
    priorMetrics=jsondecode(fileread(fullfile(prior,'05_metrics.json')));
    healthy=double(base.ground_truth);
    baseline_reconstruction=base.reconstruction;
    healthy_projections=base.projections;
    assert(isequal(size(healthy),[512 512]));
    assert(isequal(size(healthy_projections),[2000 2480]));
    lesionRelative=fullfile('lesion_insertion_project','lesion modesl','case a','S13_I137_DEL_19x18.raw');
    lesionPath=fullfile(root,lesionRelative);
    info=dir(lesionPath); assert(info.bytes==19*18*2);
    fid=fopen(lesionPath,'r','b'); assert(fid>=0);
    lesion_patch=fread(fid,[19 18],'uint16=>double'); fclose(fid);
    patch_mask=lesion_patch~=0;
    assert(nnz(patch_mask)==246);
    % Reference matrices are transposed for display. No resizing or rotation.
    patch_rows=166:184; patch_columns=251:268;
    lesion_mask=false(512); lesion_mask(patch_rows,patch_columns)=patch_mask;
    replacement=healthy(patch_rows,patch_columns);
    replacement(patch_mask)=lesion_patch(patch_mask);
    image_domain=healthy;
    image_domain(patch_rows,patch_columns)=replacement;
    lesion_delta=image_domain-healthy;
    assert(isequal(image_domain(~lesion_mask),healthy(~lesion_mask)));
    assert(isequal(image_domain(lesion_mask),lesion_patch(patch_mask)));
    assert(isequal(healthy+lesion_delta,image_domain));
    metadata=struct;
    metadata.healthy_source=priorMetrics.run_info.source;
    metadata.healthy_slice_1based=71;
    metadata.healthy_encoding=priorMetrics.run_info.encoding_used;
    metadata.lesion_source=lesionRelative;
    metadata.lesion_encoding='Inferred big-endian uint16; 684 bytes for 19 x 18 pixels. Absolute values; zero is background.';
    metadata.units='Stored pixel units, not HU or calibrated attenuation';
    metadata.mask_definition='Nonzero lesion patch; no mask interpolation';
    metadata.display_patch_rows=[251 268]; metadata.display_patch_columns=[166 184];
    metadata.placement='Fixed illustrative placement; no anatomical targeting or physical lesion-size calibration';
    metadata.lesion_spacing='Unknown; patch pixels used 1:1 on healthy 0.78125 mm grid';
    metadata.geometry=base.geometry;
    metadata.baseline_provenance=priorMetrics.provenance;
    metadata.matlab_version=version;
    metadata.formulas=struct('delta','M .* (lesion_absolute - healthy) on embedded patch; zero outside', ...
        'image_domain','healthy + delta', ...
        'projection_domain','R(A(healthy) + A(delta))', ...
        'matched_reconstruction_control','R(A(healthy + delta))', ...
        'incremental_lesion_error','[projection_domain - R(A(healthy))] - delta');
    metadata.execution='CPU copies of original curved-detector MATLAB scripts; only GPU wrappers, paths and plotting/logging changed';
    metadata.noise='No stochastic noise, Beer-Lambert conversion, HU rescaling, fitted gain/offset, or postprocessing';
    clockStart=tic;
    delta_projections=forwardReference(lesion_delta,root,cache,'delta');
    metadata.delta_forward_seconds=toc(clockStart);
    clockStart=tic;
    image_domain_projections=forwardReference(image_domain,root,cache,'combined');
    metadata.combined_forward_seconds=toc(clockStart);
    % Single precision matches the reference projection file and FBP loader.
    summed_projections=healthy_projections+delta_projections;
    clockStart=tic;
    projection_domain=backReference(summed_projections,root,cache,'sum');
    matched_image_domain=backReference(image_domain_projections,root,cache,'combined');
    reconstructed_delta=backReference(delta_projections,root,cache,'delta');
    metadata.three_backprojections_seconds=toc(clockStart);
    direct_difference=double(projection_domain)-image_domain;
    matched_difference=double(projection_domain)-double(matched_image_domain);
    lesion_increment=double(projection_domain)-double(baseline_reconstruction);
    incremental_difference=lesion_increment-lesion_delta;
    [rr,cc]=ndgrid(1:512,1:512);
    interior_mask=healthy~=0 & min(bwdist(healthy==0),min(min(rr-1,512-rr),min(cc-1,512-cc)))>5;
    ring_mask=conv2(double(lesion_mask),ones(7),'same')>0 & ~lesion_mask;
    masks=struct('all_pixels',true(512),'lesion',lesion_mask, ...
        'lesion_surrounding_ring',ring_mask,'outside_lesion',~lesion_mask,'interior',interior_mask);
    metrics=struct;
    fields=fieldnames(masks);
    for k=1:numel(fields)
        key=fields{k}; mask=masks.(key);
        metrics.direct_comparison.(key)=stats(direct_difference(mask));
        metrics.matched_reconstruction_comparison.(key)=stats(matched_difference(mask));
        metrics.incremental_lesion_error.(key)=stats(incremental_difference(mask));
    end
    metrics.projection_linearity=stats(double(summed_projections)-double(image_domain_projections));
    metrics.projection_linearity.relative_l2=norm(double(summed_projections(:))-double(image_domain_projections(:)))/norm(double(image_domain_projections(:)));
    metrics.reconstruction_additivity=stats(double(projection_domain)-double(baseline_reconstruction)-double(reconstructed_delta));
    metrics.reconstruction_additivity.relative_l2=norm(double(projection_domain(:))-double(baseline_reconstruction(:))-double(reconstructed_delta(:)))/norm(double(projection_domain(:)));
    metrics.matched_relative_l2=norm(matched_difference(:))/norm(double(matched_image_domain(:)));
    metrics.delta=stats(lesion_delta(lesion_mask));
    metrics.lesion_values=struct('minimum',min(lesion_patch(patch_mask)),'maximum',max(lesion_patch(patch_mask)),'mean',mean(lesion_patch(patch_mask)));
    metrics.healthy_values_at_lesion=struct('minimum',min(healthy(lesion_mask)),'maximum',max(healthy(lesion_mask)),'mean',mean(healthy(lesion_mask)));
    metrics.mask_notes='Interior: >5 pixels from image edges and zero padding. Ring: 3-pixel square dilation excluding lesion. No anatomical segmentation.';
    metrics.sign='Projection-domain output minus image-domain output; matched comparison subtracts reconstructed image-domain output.';
    metrics.validation=struct('replacement_exact',true,'outside_mask_unchanged',true, ...
        'projection_relative_l2_tolerance',1e-6,'matched_relative_l2_tolerance',1e-5, ...
        'matched_max_absolute_tolerance',0.25,'tolerance_units','Stored pixel units; numerical equivalence only, not clinical acceptance');
    assert(metrics.projection_linearity.relative_l2<1e-6,'Projection linearity tolerance exceeded.');
    assert(metrics.matched_relative_l2<1e-5,'Matched reconstruction tolerance exceeded.');
    assert(metrics.matched_reconstruction_comparison.all_pixels.max_absolute<0.25,'Matched maximum-error tolerance exceeded.');
    assert(metrics.reconstruction_additivity.relative_l2<1e-5,'Reconstruction additivity tolerance exceeded.');
    save(result,'healthy','baseline_reconstruction','healthy_projections','lesion_patch','lesion_mask', ...
        'lesion_delta','image_domain','delta_projections','image_domain_projections','summed_projections', ...
        'projection_domain','matched_image_domain','reconstructed_delta','direct_difference','matched_difference', ...
        'lesion_increment','incremental_difference','interior_mask','ring_mask','metadata','metrics','-v7.3');
    fprintf('PASS: paired computation, exact replacement, and numerical linearity checks.\n');
    disp(metrics.direct_comparison); disp(metrics.matched_reconstruction_comparison.lesion);
    return
end
s=load(result);
assert(isequal(s.direct_difference,double(s.projection_domain)-s.image_domain));
assert(isequal(s.matched_difference,double(s.projection_domain)-double(s.matched_image_domain)));
switch stage
    case 'comparison'
        f=figure('Visible','off','Color','white','Position',[50 50 1350 900]);
        tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
        showImage(s.healthy,'Healthy slice',false);
        showImage(s.image_domain,'Image-domain replacement',false);
        showImage(s.projection_domain,'Projection-domain insertion + FBP',false);
        showImage(s.image_domain,'Image-domain replacement: zoom',true);
        showImage(s.projection_domain,'Projection-domain insertion: zoom',true);
        showImage(s.matched_image_domain,'Image-domain replacement + FBP: zoom',true);
        sgtitle('Same healthy image, same DEL lesion, same mask and location | shared display window 7000-9000');
        exportgraphics(f,fullfile(out,'01_insertion_comparison.png'),'Resolution',140); close(f);
    case 'difference'
        f=figure('Visible','off','Color','white','Position',[50 50 1450 650]);
        tiledlayout(1,3,'Padding','loose','TileSpacing','loose');
        showDifference(s.direct_difference,200,'Projection - direct image insertion',false);
        showDifference(s.direct_difference,200,'Same subtraction: lesion zoom',true);
        lim=max(abs(s.matched_difference(:))); lim=max(lim,eps);
        showDifference(s.matched_difference,lim,'Projection - reconstructed image insertion',false);
        exportgraphics(f,fullfile(out,'02_method_differences.png'),'Resolution',140); close(f);
    case 'arrays'
        target=fullfile(out,'03_insertion_arrays.mat'); copyfile(result,target);
        check=load(target,'direct_difference','matched_difference','lesion_mask');
        assert(isequal(check.direct_difference,s.direct_difference) && isequal(check.matched_difference,s.matched_difference));
        assert(nnz(check.lesion_mask)==246);
    case 'pixels'
        [r,c]=ndgrid(1:512,1:512);
        names={'healthy','lesion_mask','lesion_delta','image_domain','projection_domain','matched_image_domain','direct_difference','matched_difference','incremental_difference'};
        t=table(r(:),c(:),'VariableNames',{'display_row','display_column'});
        for k=1:numel(names), a=double(s.(names{k}))'; t.(names{k})=a(:); end
        target=fullfile(out,'04_pixel_differences.csv'); writetable(t,target);
        check=readtable(target); assert(height(check)==262144 && nnz(check.lesion_mask)==246);
        assert(max(abs(check.direct_difference-(check.projection_domain-check.image_domain)))<1e-7);
        assert(max(abs(check.matched_difference-(check.projection_domain-check.matched_image_domain)))<1e-7);
    case 'metrics'
        payload=struct('metadata',s.metadata,'metrics',s.metrics);
        target=fullfile(out,'05_metrics.json'); fid=fopen(target,'w'); assert(fid>=0);
        fprintf(fid,'%s\n',jsonencode(payload,'PrettyPrint',true)); fclose(fid);
        check=jsondecode(fileread(target)); assert(check.metrics.direct_comparison.lesion.pixel_count==246);
    otherwise
        error('Unknown export stage: %s',stage);
end
fprintf('PASS: exported and checked %s\n',stage);
end

function p=forwardReference(inputImage,root,cache,tag)
liInput=fullfile(cache,[tag '_input.raw']); liOutput=fullfile(cache,[tag '_forward.raw']);
fid=fopen(liInput,'w','b'); assert(fid>=0); fwrite(fid,inputImage,'single'); fclose(fid);
liCode=fileread(fullfile(root,'li_fbp','c1_forward_projection_v6_curved_detector.m'));
liCode=strrep(liCode,'.\without_lesion.raw',liInput);
liCode=strrep(liCode,'.\projections.raw',liOutput);
liCode=strrep(liCode,'gpuArray(','('); liCode=strrep(liCode,'gather(','(');
liCode=regexprep(liCode,'(?m)^figure;[^\r\n]*','');
liCode=regexprep(liCode,'(?m)^title\([^\r\n]*','');
liCode=strrep(liCode,"fprintf('View %d / %d\n', view, n_view);", ...
    "if mod(view,250)==0, fprintf('Forward %s view %d/%d\n',tag,view,n_view); end");
eval(liCode);
assert(isequal(vol,double(inputImage)),'Float32 input conversion altered values.');
fid=fopen(liOutput,'r','l'); p=fread(fid,[2480 2000],'single=>single')'; fclose(fid);
assert(isequal(size(p),[2000 2480]) && all(isfinite(p(:))));
end

function r=backReference(inputProjections,root,cache,tag)
liInput=fullfile(cache,[tag '_fbp_input.raw']); liOutput=fullfile(cache,[tag '_reconstruction.raw']);
fid=fopen(liInput,'w','l'); assert(fid>=0); fwrite(fid,inputProjections','single'); fclose(fid);
liCode=fileread(fullfile(root,'li_fbp','c2_fbp_v6_curved_detector.m'));
liCode=strrep(liCode,'gpuDevice(1);','');
liCode=strrep(liCode,'''projections.raw''',['''' liInput '''']);
liCode=strrep(liCode,'''img.raw''',['''' liOutput '''']);
liCode=strrep(liCode,'gpuArray.zeros(','zeros(');
liCode=strrep(liCode,', ''gpuArray''','');
liCode=strrep(liCode,'gpuArray(','('); liCode=strrep(liCode,'gather(','(');
liCode=regexprep(liCode,'(?m)^(figure;|imagesc\([^\r\n]*|colormap gray;[^\r\n]*|title\([^\r\n]*)','');
liCode=strrep(liCode,"fprintf('View %d / %d\n', idx_angle, views);",'');
sdd=single(1200); % Original forward-script variable used by c2 for d_gamma.
eval(liCode);
assert(isequal(projection_data,inputProjections));
r=img; assert(isequal(size(r),[512 512]) && all(isfinite(r(:))));
end

function s=stats(values)
d=double(values(:)); assert(~isempty(d) && all(isfinite(d)));
a=sort(abs(d));
s=struct('pixel_count',numel(d),'mae',mean(a),'rmse',sqrt(mean(d.^2)), ...
    'bias',mean(d),'max_absolute',max(a),'minimum_signed',min(d),'maximum_signed',max(d), ...
    'median_absolute',median(a),'p95_absolute',interp1(1:numel(a),a,1+0.95*(numel(a)-1)));
end

function showImage(a,label,zoomed)
ax=nexttile; imagesc(ax,a',[7000 9000]); axis(ax,'image'); colormap(ax,gray(256)); colorbar(ax);
hold(ax,'on'); rectangle(ax,'Position',[165.5 250.5 19 18],'EdgeColor',[1 .3 .1],'LineWidth',1.2);
if zoomed, xlim(ax,[145.5 204.5]); ylim(ax,[230.5 288.5]); else, axis(ax,'off'); end
title(ax,label);
end

function showDifference(a,lim,label,zoomed)
ax=nexttile; imagesc(ax,a',[-lim lim]); axis(ax,'image'); colormap(ax,gray(256)); colorbar(ax);
if zoomed, xlim(ax,[145.5 204.5]); ylim(ax,[230.5 288.5]); else, axis(ax,'off'); end
title(ax,{label,sprintf('Stored pixel units | display +/-%.4g | gray = zero',lim)});
end
