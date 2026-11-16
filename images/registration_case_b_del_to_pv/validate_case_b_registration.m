function validate_case_b_registration(stage,cache)
% Real-case validation of the existing demons function, without changing it.
% Stages: compute, comparison, differences, slab, arrays, metrics.
% All axes use display orientation (transpose of the raw MATLAB slice).
if nargin<2,cache=fullfile(tempdir,'cbct_registration_case_b');end
out=fileparts(mfilename('fullpath'));root=fileparts(fileparts(out));addpath(root);
if ~isfolder(cache),mkdir(cache);end
result=fullfile(cache,'validated.mat');
if strcmp(stage,'compute')
    pvPath=fullfile(root,'lesion_insertion_project/cases without lesion/case b/ABD PEL AX 3 PV_40cm_512x512x171.raw');
    delPath=fullfile(root,'lesion_insertion_project/cases without lesion/case b/ABD PEL AX 3 DEL_40cm_512x512x142.raw');
    pv=readVolume(pvPath,171);del=readVolume(delPath,142);
    moving=double(del(:,:,71));
    % Data-derived slice correspondence; not a DICOM coordinate transform.
    a=sliceFeature(moving);scores=zeros(171,1);
    for k=1:171
        c=normxcorr2(a,sliceFeature(pv(:,:,k)));
        scores(k)=max(c(121:135,121:135),[],'all');
    end
    [~,fixedSlice]=max(scores);fixed=double(pv(:,:,fixedSlice));
    assert(fixedSlice==91,'Case B slice-match result changed; inspect inputs.');
    t=tic;[registered,D,info]=demons_iodine_registration(fixed,moving);seconds=toc(t);
    body=bodyMask(fixed);movingBody=bodyMask(moving);
    roi=imerode(body,strel('disk',5,0)) & info.ValidMask;
    fixedEdges=edges(fixed) & roi;movingEdges=edges(moving) & roi;
    warpedEdges=imwarp(movingEdges,D,'nearest') & roi;
    newEdges=edges(registered) & roi;
    metrics=struct;
    metrics.phase_2d=phaseMetrics(fixed,moving,registered,D,info,body,movingBody,roi);
    metrics.phase_2d.edge_distance_before_pixels=edgeDistance(fixedEdges,movingEdges);
    metrics.phase_2d.propagated_edge_distance_after_pixels=edgeDistance(fixedEdges,warpedEdges);
    metrics.phase_2d.redetected_edge_distance_after_pixels=edgeDistance(fixedEdges,newEdges);
    metrics.phase_2d.seconds=seconds;
    % One bounded feature-setting diagnostic; default remains the main result.
    [trialRegistered,trialD,trialInfo]=demons_iodine_registration(fixed,moving, ...
        'RegistrationPair',{window(fixed),window(moving)},'FieldSmoothing',2);
    common=roi & trialInfo.ValidMask;
    metrics.feature_trial=phaseMetrics(fixed,moving,trialRegistered,trialD,trialInfo,body,movingBody,common);
    metrics.feature_trial.propagated_edge_distance_after_pixels=edgeDistance(fixedEdges & common,imwarp(movingEdges,trialD,'nearest') & common);
    metrics.feature_trial.redetected_edge_distance_after_pixels=edgeDistance(fixedEdges & common,edges(trialRegistered) & common);
    % Known-motion control on the actual PV image, not a patient motion truth.
    translated=imtranslate(fixed,[4 -3],'OutputView','same','FillValues',0);
    [knownRegistered,knownD,knownInfo]=demons_iodine_registration(fixed,translated);
    knownROI=imerode(body,strel('disk',10,0)) & knownInfo.ValidMask;
    endpoint=hypot(knownD(:,:,1)-4,knownD(:,:,2)+3);
    metrics.known_translation=struct('expected_backward_dx',4,'expected_backward_dy',-3, ...
        'median_endpoint_error_pixels',median(endpoint(knownROI)), ...
        'p95_endpoint_error_pixels',prctile(endpoint(knownROI),95), ...
        'rmse_before',rmse(translated(knownROI)-fixed(knownROI)), ...
        'rmse_after',rmse(knownRegistered(knownROI)-fixed(knownROI)), ...
        'min_jacobian',knownInfo.MinJacobian,'fold_fraction',knownInfo.NonpositiveJacobianFraction);
    % Local 3D coverage at half in-plane resolution, preserving 17 slice planes.
    % Unknown slice origin/spacing prevents physical-coordinate claims.
    fixed3=imresize3(single(pv(:,:,83:99)),[256 256 17],'linear');
    moving3=imresize3(single(del(:,:,63:79)),[256 256 17],'linear');
    t=tic;[registered3,D3,info3]=demons_iodine_registration(fixed3,moving3);seconds3=toc(t);
    roi3=fixed3>7700 & info3.ValidMask;roi3(:,:,[1:2 16:17])=false;
    metrics.slab_3d=struct('fixed_slices',[83 99],'moving_slices',[63 79], ...
        'grid',[256 256 17],'seconds',seconds3,'evaluation_voxels',nnz(roi3), ...
        'feature_rmse_before',info3.FeatureRMSEBefore,'feature_rmse_after',info3.FeatureRMSEAfter, ...
        'phase_rmse_before',rmse(double(moving3(roi3))-double(fixed3(roi3))), ...
        'phase_rmse_after',rmse(registered3(roi3)-double(fixed3(roi3))), ...
        'min_jacobian',info3.MinJacobian,'fold_fraction',info3.NonpositiveJacobianFraction, ...
        'valid_fraction',info3.ValidFraction);
    % Identity verifies unchanged values and displacement convention.
    [identity,identityD,identityInfo]=demons_iodine_registration(fixed,fixed);
    assert(isequal(identity,fixed) && all(identityD(:)==0) && all(identityInfo.ValidMask(:)));
    check=imwarp(moving,D,'linear','FillValues',0,'SmoothEdges',false);
    assert(isequal(check(info.ValidMask),registered(info.ValidMask)));
    assert(all(isnan(registered(~info.ValidMask))));
    assert(info.NonpositiveJacobianFraction==0 && info3.NonpositiveJacobianFraction==0);
    assert(metrics.phase_2d.body_dice_after>metrics.phase_2d.body_dice_before);
    assert(metrics.phase_2d.bright_structure_dice_after>metrics.phase_2d.bright_structure_dice_before);
    assert(metrics.phase_2d.propagated_edge_distance_after_pixels<metrics.phase_2d.edge_distance_before_pixels);
    assert(info.FeatureRMSEAfter<info.FeatureRMSEBefore && info3.FeatureRMSEAfter<info3.FeatureRMSEBefore);
    assert(metrics.known_translation.median_endpoint_error_pixels<1 && knownInfo.NonpositiveJacobianFraction==0);
    assert(metrics.known_translation.rmse_after<metrics.known_translation.rmse_before);
    metrics.checks=struct('identity_exact',true,'output_is_original_intensity_warp',true, ...
        'outside_support_is_nan',true,'structural_overlap_improved',true,'no_detected_folds_2d_or_3d',true, ...
        'known_translation_median_error_below_one_pixel',true, ...
        'redetected_edge_metric_improved',metrics.phase_2d.redetected_edge_distance_after_pixels<metrics.phase_2d.edge_distance_before_pixels);
    metadata=struct('fixed_phase','portal venous','moving_phase','delayed', ...
        'fixed_source',pvPath,'moving_source',delPath,'fixed_slice',fixedSlice,'moving_slice',71, ...
        'encoding','Inferred big-endian uint16, consistent with two bytes per voxel. Reported float32 metadata conflicts with file length.', ...
        'units','Stored pixel units, not confirmed HU; displacement in tested-grid pixels/voxels', ...
        'orientation','Raw MATLAB matrices transposed for display before registration', ...
        'inplane_fov_mm',400,'native_inplane_spacing_mm',400/512,'reported_slice_thickness_mm',3, ...
        'spatial_limit','No DICOM origins/orientations or confirmed inter-slice spacing. Content-matched 2D pair and local 3D slab only; not whole-volume physical registration.', ...
        'phase_limit','Both inputs are post-contrast. No pre-injection image is available.', ...
        'function','demons_iodine_registration.m, unchanged, default settings', ...
        'settings',info.Options,'matlab_version',version, ...
        'slice_selection','Maximum bone-focused cross-correlation against DEL 71, scanning all 171 PV slices with up to +/-28 native pixels of translation', ...
        'feature_thresholds','Slice match: clip stored values minus 8450 to [0,1000]. Body proxy >7700, bright-structure proxy >8650. Edges use clipped [7400,8800], Gaussian sigma 1, Canny [0.05,0.15]. These are analysis settings, not HU or scanner parameter changes.', ...
        'edge_limit','Propagated-edge distance measures geometric agreement of source edges. Redetecting edges after interpolation gives a worse result and is also reported. Neither is ground-truth landmark error.', ...
        'feature_trial','Same function, clipped [7400,8800] feature copies and FieldSmoothing=2; retained as diagnostic, not selected as main result');
    save(result,'fixed','moving','registered','D','info','body','movingBody','roi', ...
        'fixedEdges','movingEdges','warpedEdges','newEdges','scores','trialRegistered','trialD','trialInfo', ...
        'translated','knownRegistered','knownD','knownROI','fixed3','moving3','registered3','D3','info3','roi3','metrics','metadata','-v7.3');
    disp(metrics);fprintf('PASS: real-case 2D/3D numerical checks and known-motion control.\n');return
end
s=load(result);
switch stage
    case 'comparison'
        f=figure('Visible','off','Color','white','Position',[50 50 1250 850]);tiledlayout(2,3,'Padding','loose');
        grayPanel(s.fixed,'PV fixed: slice 91');grayPanel(s.moving,'DEL before: slice 71');grayPanel(s.registered,'DEL after demons');
        overlayPanel(s.fixed,s.moving,'Before: PV green / DEL magenta');
        overlayPanel(s.fixed,s.registered,'After: PV green / DEL magenta');
        ax=nexttile;imagesc(ax,hypot(s.D(:,:,1),s.D(:,:,2)));axis(ax,'image');axis(ax,'off');colormap(ax,parula);colorbar(ax);title(ax,'Backward displacement (pixels)');
        exportgraphics(f,fullfile(out,'01_registration_comparison.png'),'Resolution',130);close(f);
    case 'differences'
        f=figure('Visible','off','Color','white','Position',[50 50 1200 560]);tiledlayout(1,2,'Padding','loose');
        a=s.moving-s.fixed;a(~s.info.ValidMask)=NaN;
        b=s.registered-s.fixed;
        for k=1:2
            ax=nexttile;if k==1,p=a;label='Before: DEL - PV';else,p=b;label='After: registered DEL - PV';end
            h=imagesc(ax,p,[-200 200]);h.AlphaData=isfinite(p);set(ax,'Color','black');axis(ax,'image');axis(ax,'off');colormap(ax,gray(256));colorbar(ax);
            title(ax,{label,'Display +/-200 stored units; includes genuine phase differences'});
        end
        exportgraphics(f,fullfile(out,'02_phase_subtractions.png'),'Resolution',130);close(f);
    case 'slab'
        f=figure('Visible','off','Color','white','Position',[50 50 1250 800]);tiledlayout(2,3,'Padding','loose');
        grayPanel(s.fixed3(:,:,9),'3D slab: PV center');grayPanel(s.moving3(:,:,9),'3D slab: DEL before');grayPanel(s.registered3(:,:,9),'3D slab: DEL after');
        grayPanel(s.fixed,'Known-motion control: PV reference');grayPanel(s.translated,'Known shift: +4 x, -3 y pixels');grayPanel(s.knownRegistered,'Recovered from known shift');
        exportgraphics(f,fullfile(out,'03_slab_and_motion_control.png'),'Resolution',130);close(f);
    case 'arrays'
        target=fullfile(out,'04_registration_arrays.mat');copyfile(result,target);
        check=load(target,'registered','D','metrics');assert(isequaln(check.registered,s.registered) && isequal(check.D,s.D));
    case 'metrics'
        target=fullfile(out,'05_metrics.json');fid=fopen(target,'w');assert(fid>=0);
        fprintf(fid,'%s\n',jsonencode(struct('metadata',s.metadata,'metrics',s.metrics),'PrettyPrint',true));fclose(fid);
        check=jsondecode(fileread(target));assert(check.metrics.phase_2d.fold_fraction==0);
    otherwise,error('Unknown stage: %s',stage);
end
fprintf('PASS: exported and checked %s\n',stage);
end

function v=readVolume(path,n)
entry=dir(path);assert(entry.bytes==512*512*n*2);
fid=fopen(path,'r','b');assert(fid>=0);v=reshape(fread(fid,inf,'uint16=>uint16'),512,512,n);fclose(fid);v=permute(v,[2 1 3]);
end
function a=sliceFeature(v)
a=imresize(min(max(single(v)-8450,0),1000),[128 128]);
a(1:10,:)=0;a(110:end,:)=0;a(:,1:10)=0;a(:,119:end)=0;
end
function b=bodyMask(a),b=imfill(bwareafilt(a>7700,1),'holes');end
function a=window(a),a=min(max((double(a)-7400)/1400,0),1);a(~isfinite(a))=0;end
function e=edges(a),e=edge(imgaussfilt(window(a),1),'Canny',[.05 .15]);end
function d=edgeDistance(a,b)
da=bwdist(a);db=bwdist(b);d=mean([da(b);db(a)]);
end
function x=rmse(a),x=sqrt(mean(double(a(:)).^2));end
function m=phaseMetrics(f,moving,r,D,info,body,movingBody,roi)
warpedBody=imwarp(movingBody,D,'nearest');
fixedBright=f>8650 & body;movingBright=moving>8650 & movingBody;
warpedBright=imwarp(movingBright,D,'nearest');
disp=hypot(D(:,:,1),D(:,:,2));
m=struct('evaluation_pixels',nnz(roi),'phase_rmse_before',rmse(moving(roi)-f(roi)), ...
    'phase_rmse_after',rmse(r(roi)-f(roi)), ...
    'phase_mae_before',mean(abs(moving(roi)-f(roi))),'phase_mae_after',mean(abs(r(roi)-f(roi))), ...
    'feature_rmse_before',info.FeatureRMSEBefore,'feature_rmse_after',info.FeatureRMSEAfter, ...
    'body_dice_before',dice(body & info.ValidMask,movingBody & info.ValidMask), ...
    'body_dice_after',dice(body & info.ValidMask,warpedBody & info.ValidMask), ...
    'bright_structure_dice_before',dice(fixedBright & info.ValidMask,movingBright & info.ValidMask), ...
    'bright_structure_dice_after',dice(fixedBright & info.ValidMask,warpedBright & info.ValidMask), ...
    'min_jacobian',info.MinJacobian,'fold_fraction',info.NonpositiveJacobianFraction, ...
    'valid_fraction',info.ValidFraction,'max_displacement_pixels',max(disp(:)), ...
    'body_median_displacement_pixels',median(disp(roi)), ...
    'body_p95_displacement_pixels',prctile(disp(roi),95));
end
function grayPanel(a,label)
ax=nexttile;h=imagesc(ax,a,[7000 9000]);h.AlphaData=isfinite(a);set(ax,'Color','black');axis(ax,'image');axis(ax,'off');colormap(ax,gray(256));colorbar(ax);title(ax,label);
end
function overlayPanel(f,m,label)
ax=nexttile;a=min(max((double(f)-7000)/2000,0),1);b=min(max((double(m)-7000)/2000,0),1);
valid=isfinite(b);a(~valid)=0;b(~valid)=0;image(ax,cat(3,b,a,b));axis(ax,'image');axis(ax,'off');title(ax,label);
end
