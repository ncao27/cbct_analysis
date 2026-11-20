function validate_direct_reference(stage, cache, destination)
% Direct reference projection -> FBP validation, without baseline correction.
% Preserves li_fbp; corrected FBP is a separately labeled comparison.
% Run compute, then reconstruction, insertion, differences, arrays, metrics.
% Optional destination allows verification before publishing one file at a time.
out=fileparts(mfilename('fullpath')); root=fileparts(fileparts(out));
if nargin<2, cache=fullfile(tempdir,'cbct_direct_reference'); end
if nargin<3, destination=out; end
if ~isfolder(cache), mkdir(cache); end
result=fullfile(cache,'computed.mat');
if strcmp(stage,'compute')
    b=load(fullfile(root,'images/lesion_free_roundtrip_case_b_del_slice071/03_roundtrip_arrays.mat'));
    l=load(fullfile(root,'images/lesion_insertion_case_b_del_slice071/03_insertion_arrays.mat'));
    healthy=double(b.ground_truth); image_domain=double(l.image_domain);
    lesion_mask=l.lesion_mask; lesion_patch=l.lesion_patch;
    interior_mask=b.interior_mask; lesion_delta=l.lesion_delta;
    assert(isequal(healthy,l.healthy) && isequal(healthy+lesion_delta,image_domain));
    assert(nnz(lesion_mask)==246 && isequal(healthy(~lesion_mask),image_domain(~lesion_mask)));
    original_healthy=b.reconstruction; original_insertion=l.projection_domain;
    baseline_check=back(b.projections,root,cache,'baseline_check',false);
    assert(isequal(baseline_check,original_healthy),'Reference baseline no longer matches saved results.');
    % Use the original stored-unit projections directly; no image-derived correction.
    healthy_reconstruction=double(back(b.projections,root,cache,'healthy',true));
    projection_domain=double(back(b.projections+l.delta_projections,root,cache,'insertion',true));
    matched_image_domain=double(back(l.image_domain_projections,root,cache,'matched',true));
    reconstructed_delta=double(back(l.delta_projections,root,cache,'delta',true));
    healthy_error=healthy_reconstruction-healthy;
    insertion_error=projection_domain-image_domain;
    matched_error=projection_domain-matched_image_domain;
    incremental_error=(projection_domain-healthy_reconstruction)-lesion_delta;
    masks=struct('all_pixels',true(512),'interior',interior_mask,'lesion',lesion_mask);
    names=fieldnames(masks); metrics=struct;
    for j=1:numel(names)
        key=names{j}; mask=masks.(key);
        metrics.original_healthy.(key)=stats(double(original_healthy)-healthy,mask);
        metrics.corrected_healthy.(key)=stats(healthy_error,mask);
        metrics.original_insertion.(key)=stats(double(original_insertion)-image_domain,mask);
        metrics.corrected_insertion.(key)=stats(insertion_error,mask);
        metrics.matched_methods.(key)=stats(matched_error,mask);
        metrics.incremental_lesion.(key)=stats(incremental_error,mask);
    end
    additivity=projection_domain-healthy_reconstruction-reconstructed_delta;
    metrics.reconstruction_additivity=stats(additivity,true(512));
    assert(metrics.corrected_healthy.interior.rmse<metrics.original_healthy.interior.rmse);
    assert(metrics.corrected_insertion.interior.rmse<metrics.original_insertion.interior.rmse);
    assert(metrics.matched_methods.all_pixels.max_absolute<0.25,'Matched-method numerical discrepancy exceeds 0.25 stored units.');
    assert(metrics.reconstruction_additivity.max_absolute<0.25,'Numerical additivity check failed.');
    metadata=struct('case','Case B DEL slice 71','lesion','Case A S13_I137_DEL_19x18.raw', ...
        'units','Stored pixel units; inferred raw encoding; not confirmed HU or attenuation', ...
        'method','Reference: R(A(healthy)); comparison: R_corrected(A(healthy))', ...
        'insertion','R_corrected(A(healthy)+A(lesion_delta))', ...
        'matched_control','R_corrected(A(image_domain))', ...
        'corrections','Comparison only: normalized detector interpolation and zero invalid detector contributions', ...
        'limitation','No computational offset, rescaling, baseline restoration, or image-derived artifact correction; raw encoding and calibration remain unverified.', ...
        'smoothing','None; no median, Gaussian, or other image denoising', ...
        'execution','Reuse saved projections from professor curved-detector projector; rerun CPU FBP; geometry and sampling unchanged; GPU not tested', ...
        'interior_roi','Original mask: nonzero pixels more than 5 pixels from image border and zero padding; includes background', ...
        'error_sign','Reconstruction minus target; matched error is projection insertion minus reconstructed image insertion', ...
        'display_offset',8192,'display_window',[-200 300],'geometry',b.geometry,'matlab_version',version);
    save(result,'healthy','image_domain','lesion_mask','lesion_patch','lesion_delta','interior_mask', ...
        'original_healthy','original_insertion','healthy_reconstruction','projection_domain', ...
        'matched_image_domain','reconstructed_delta','healthy_error','insertion_error','matched_error', ...
        'incremental_error','metadata','metrics','-v7.3');
    fprintf('PASS: original reference exact, direct input projections exact, interpolation comparison errors reduced, matched methods and additivity verified.\n');
    disp(metrics.corrected_healthy.interior); disp(metrics.corrected_insertion.interior);
    disp(metrics.corrected_insertion.lesion); disp(metrics.matched_methods.all_pixels);
    return
end
s=load(result);
assert(isequal(s.healthy_error,s.healthy_reconstruction-s.healthy));
assert(isequal(s.insertion_error,s.projection_domain-s.image_domain));
assert(isequal(s.matched_error,s.projection_domain-s.matched_image_domain));
if ~isfolder(destination),mkdir(destination);end
switch stage
    case 'arrays'
        target=fullfile(destination,'04_arrays.mat');copyfile(result,target);
        q=load(target,'healthy_reconstruction','projection_domain','lesion_mask');
        assert(isequal(q.healthy_reconstruction,s.healthy_reconstruction));
        assert(isequal(q.projection_domain,s.projection_domain) && nnz(q.lesion_mask)==246);
    case 'metrics'
        target=fullfile(destination,'05_metrics.json');fid=fopen(target,'w');assert(fid>=0);
        fprintf(fid,'%s\n',jsonencode(struct('metadata',s.metadata,'metrics',s.metrics),'PrettyPrint',true));fclose(fid);
        q=jsondecode(fileread(target));assert(q.metrics.corrected_insertion.lesion.pixel_count==246);
    case {'reconstruction','insertion','differences'}
        f=figure('Visible','off','Color','white','Position',[30 30 1600 1050]);
        cleanup=onCleanup(@()close(f));
        layout=tiledlayout(f,2,3,'Padding','loose','TileSpacing','loose');
        title(layout,'Case B DEL slice 71 | direct projection to FBP | no baseline correction','FontSize',15);
        subtitle(layout,{'Professor forward projector retained. Reference FBP and corrected FBP shown separately. No image smoothing.', ...
            'Image window [-200, 300] after subtracting 8192 for display; subtraction maps use stored units.'},'FontSize',11);
        if strcmp(stage,'reconstruction')
            ct(s.healthy,'Healthy input',false,false);
            ct(s.original_healthy,'Reference FBP (unchanged mathematics)',false,false);
            ct(s.healthy_reconstruction,'FBP with interpolation/index corrections',false,false);
            ct(s.healthy,'Healthy input: detail',true,false);
            ct(s.healthy_reconstruction,'Corrected FBP: detail',true,false);
            difference(s.healthy-s.healthy_reconstruction,50,sprintf('Input - corrected FBP | interior RMSE %.2f',s.metrics.corrected_healthy.interior.rmse),false);
            target=fullfile(destination,'01_reconstruction.png');
        elseif strcmp(stage,'insertion')
            ct(s.image_domain,'Image-domain replacement',false,true);
            ct(s.original_insertion,'Projection insertion: reference FBP',false,true);
            ct(s.projection_domain,'Projection insertion: corrected FBP',false,true);
            ct(s.image_domain,'Image-domain replacement: lesion detail',true,true);
            ct(s.projection_domain,'Projection insertion: corrected FBP detail',true,true);
            ct(s.matched_image_domain,'Image replacement + corrected FBP: detail',true,true);
            target=fullfile(destination,'02_lesion_insertion.png');
        else
            difference(double(s.original_healthy)-s.healthy,200,'Original healthy reconstruction error',false);
            difference(s.healthy_error,200,'Corrected FBP error: same scale',false);
            difference(s.healthy_error,50,'Corrected FBP error: tighter scale',false);
            difference(s.insertion_error,50,'Projection - direct image insertion',true);
            difference(s.incremental_error,50,'Incremental lesion error',true);
            difference(s.matched_error,0.1,'Projection - reconstructed image insertion',true);
            target=fullfile(destination,'03_pixel_differences.png');
        end
        exportgraphics(f,target,'Resolution',130);
    otherwise,error('Unknown stage: %s',stage);
end
fprintf('PASS: exported %s\n',target);
end

function r=back(inputProjections,root,cache,tag,corrected)
liInput=fullfile(cache,[tag '_projections.raw']);liOutput=fullfile(cache,[tag '_reconstruction.raw']);
fid=fopen(liInput,'w','l');assert(fid>=0);fwrite(fid,inputProjections','single');fclose(fid);
code=fileread(fullfile(root,'li_fbp/c2_fbp_v6_curved_detector.m'));
code=strrep(code,'gpuDevice(1);','');
code=strrep(code,'''projections.raw''',['''' liInput '''']);code=strrep(code,'''img.raw''',['''' liOutput '''']);
code=strrep(code,'gpuArray.zeros(','zeros(');code=strrep(code,', ''gpuArray''','');
code=strrep(code,'gpuArray(','(');code=strrep(code,'gather(','(');
code=regexprep(code,'(?m)^(figure;|imagesc\([^\r\n]*|colormap gray;[^\r\n]*|title\([^\r\n]*)','');
code=strrep(code,"fprintf('View %d / %d\n', idx_angle, views);",'');
if corrected
    old='w = gamma0 - gamma(k);';assert(numel(strfind(code,old))==1);
    code=strrep(code,old,'w = (gamma0 - gamma(k)) ./ d_gamma;');
    old='img = img + SID ./ (L.^2) .* val * d_theta;';assert(numel(strfind(code,old))==1);
    code=strrep(code,old,sprintf('val(~valid) = 0;\n    %s',old));
end
sdd=single(1200); % Preserve the original forward-script workspace value.
eval(code);r=img;
assert(isequal(projection_data,inputProjections) && all(isfinite(r(:))));
fprintf('Backprojection complete: %s\n',tag);
end

function m=stats(error,mask)
e=double(error(mask));assert(~isempty(e) && all(isfinite(e)));
m=struct('pixel_count',numel(e),'mae',mean(abs(e)),'rmse',sqrt(mean(e.^2)), ...
    'bias',mean(e),'max_absolute',max(abs(e)));
end

function ct(a,label,zoomed,mark)
ax=nexttile;imagesc(ax,double(a)'-8192,[-200 300]);axis(ax,'image');axis(ax,'off');
colormap(ax,gray(256));c=colorbar(ax);c.Label.String='Stored value - 8192 (assumed)';title(ax,label);
if mark,hold(ax,'on');rectangle(ax,'Position',[165.5 250.5 19 18],'EdgeColor',[1 .3 .1],'LineWidth',1.1);end
if zoomed,xlim(ax,[145.5 204.5]);ylim(ax,[230.5 288.5]);end
end

function difference(e,limit,label,zoomed)
ax=nexttile;imagesc(ax,e',[-limit limit]);axis(ax,'image');axis(ax,'off');colormap(ax,gray(256));
c=colorbar(ax);c.Label.String='Difference in stored units';title(ax,{label,sprintf('Window +/- %g',limit)});
if zoomed,xlim(ax,[145.5 204.5]);ylim(ax,[230.5 288.5]);end
end
