function render_display_previews(stage,assumedOffset)
% Display-only previews of saved results. This does NOT establish HU units.
% Default: displayed value = stored value - 8192, window [-200 300].
% The offset is an explicit assumption because raw export calibration is absent.
% Stages: check, reconstruction, insertion, lesion, registration, controls.
% Reads saved MAT files only; no reconstruction/registration or data rewrite.
if nargin<2,assumedOffset=8192;end
validateattributes(assumedOffset,{'numeric'},{'real','finite','scalar'});
out=fileparts(mfilename('fullpath'));imagesRoot=fileparts(out);
roundtrip=fullfile(imagesRoot,'lesion_free_roundtrip_case_b_del_slice071','03_roundtrip_arrays.mat');
insertion=fullfile(imagesRoot,'lesion_insertion_case_b_del_slice071','03_insertion_arrays.mat');
registration=fullfile(imagesRoot,'registration_case_b_del_to_pv','04_registration_arrays.mat');
limits=[-200 300];
if strcmp(stage,'check')
    s=load(roundtrip,'ground_truth','reconstruction','signed_error');
    shiftedDifference=(double(s.reconstruction)-assumedOffset)-(double(s.ground_truth)-assumedOffset);
    assert(max(abs(shiftedDifference(:)-s.signed_error(:)))<1e-9);
    assert(isequal((limits+assumedOffset)-assumedOffset,limits));
    assert(isequal(displayGray(limits+assumedOffset,assumedOffset),[0 1]));
    r=load(registration,'registered');
    assert(isequal(isfinite(r.registered),isfinite(double(r.registered)-assumedOffset)));
    fprintf('PASS: display mapping, endpoint clipping, NaN support, and subtraction invariance.\n');return
end
if any(strcmp(stage,{'insertion','registration','controls'})),shape=[2 3];sz=[1400 950];
elseif strcmp(stage,'lesion'),shape=[1 1];sz=[650 650];
else,shape=[1 3];sz=[1450 650];end
f=figure('Visible','off','Color','white','Position',[50 50 sz]);
cleanup=onCleanup(@()close(f));
layout=tiledlayout(f,shape(1),shape(2),'Padding','loose','TileSpacing','loose');
title(layout,sprintf('Window [-200, 300] | stored value minus %g | assumed offset, NOT confirmed HU',assumedOffset),'FontSize',12);
switch stage
    case 'reconstruction'
        s=load(roundtrip,'ground_truth','reconstruction','signed_error');
        ctPanel(s.ground_truth','Healthy ground truth',assumedOffset,false);
        ctPanel(s.reconstruction','Reconstruction',assumedOffset,false);
        ax=nexttile;imagesc(ax,-s.signed_error',[-200 200]);axis(ax,'image');axis(ax,'off');colormap(ax,gray(256));
        c=colorbar(ax);c.Label.String='Difference in stored units';
        title(ax,{'Ground truth - reconstruction','Difference window remains +/-200'});
        target='01_reconstruction.png';
    case 'insertion'
        s=load(insertion,'healthy','image_domain','projection_domain','matched_image_domain');
        ctPanel(s.healthy','Healthy slice',assumedOffset,true);
        ctPanel(s.image_domain','Image-domain replacement',assumedOffset,true);
        ctPanel(s.projection_domain','Projection-domain insertion + FBP',assumedOffset,true);
        ax=ctPanel(s.image_domain','Image-domain replacement: zoom',assumedOffset,true);zoomPatch(ax);
        ax=ctPanel(s.projection_domain','Projection-domain insertion: zoom',assumedOffset,true);zoomPatch(ax);
        ax=ctPanel(s.matched_image_domain','Image-domain replacement + FBP: zoom',assumedOffset,true);zoomPatch(ax);
        target='02_lesion_insertion.png';
    case 'lesion'
        s=load(insertion,'lesion_patch');
        assert(isequal(size(s.lesion_patch),[19 18]) && nnz(s.lesion_patch)==246);
        ctPanel(s.lesion_patch',{'Case A: S13 I137 DEL','19 x 18 pixels; 246 lesion pixels'},assumedOffset,false);
        target='03_inserted_lesion.png';
    case 'registration'
        s=load(registration,'fixed','moving','registered','D');
        ctPanel(s.fixed,'PV fixed: slice 91',assumedOffset,false);
        ctPanel(s.moving,'DEL before: slice 71',assumedOffset,false);
        ctPanel(s.registered,'DEL after registration',assumedOffset,false);
        overlayPanel(s.fixed,s.moving,'Before: PV green / DEL magenta',assumedOffset);
        overlayPanel(s.fixed,s.registered,'After: PV green / DEL magenta',assumedOffset);
        ax=nexttile;imagesc(ax,hypot(s.D(:,:,1),s.D(:,:,2)));axis(ax,'image');axis(ax,'off');colormap(ax,parula);colorbar(ax);
        title(ax,'Backward displacement (pixels; unchanged)');
        target='04_registration.png';
    case 'controls'
        s=load(registration,'fixed3','moving3','registered3','fixed','translated','knownRegistered');
        ctPanel(s.fixed3(:,:,9),'3D slab: PV center',assumedOffset,false);
        ctPanel(s.moving3(:,:,9),'3D slab: DEL before',assumedOffset,false);
        ctPanel(s.registered3(:,:,9),'3D slab: DEL after',assumedOffset,false);
        ctPanel(s.fixed,'Known-motion control: PV',assumedOffset,false);
        ctPanel(s.translated,'Known shift: +4 x, -3 y pixels',assumedOffset,false);
        ctPanel(s.knownRegistered,'Recovered from known shift',assumedOffset,false);
        target='05_registration_controls.png';
    otherwise,error('Unknown stage: %s',stage);
end
drawnow;
exportgraphics(f,fullfile(out,target),'Resolution',130);
fprintf('Saved %s; display assumption only, original arrays unchanged.\n',target);
end

function ax=ctPanel(a,label,offset,markPatch)
ax=nexttile;shifted=double(a)-offset;
h=imagesc(ax,shifted,[-200 300]);h.AlphaData=isfinite(shifted);
axis(ax,'image');axis(ax,'off');colormap(ax,gray(256));c=colorbar(ax);
c.Label.String=sprintf('Stored value - %g (assumed)',offset);
assert(isequal(ax.CLim,[-200 300]));
assert(isequaln(h.CData,shifted));
title(ax,label);
if markPatch
    hold(ax,'on');rectangle(ax,'Position',[165.5 250.5 19 18],'EdgeColor',[1 .3 .1],'LineWidth',1.2);
end
end
function zoomPatch(ax)
xlim(ax,[145.5 204.5]);ylim(ax,[230.5 288.5]);
end
function rgb=displayGray(a,offset)
rgb=min(max((double(a)-offset+200)/500,0),1);
end
function overlayPanel(f,m,label,offset)
ax=nexttile;a=displayGray(f,offset);b=displayGray(m,offset);
valid=isfinite(f) & isfinite(m);a(~valid)=0;b(~valid)=0;
image(ax,cat(3,b,a,b));axis(ax,'image');axis(ax,'off');title(ax,label);
end
