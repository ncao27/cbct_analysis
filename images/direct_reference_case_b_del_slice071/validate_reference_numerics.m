function validate_reference_numerics(stage,cache,destination)
% Research validation of the professor's cubic ray integration and fan FBP.
% No baseline subtraction, source-derived residual correction, or fitted gain.
% Geometry unchanged. Boundary extension and reconstruction bandwidth explicit.
% Stages: compute, reconstruction, insertion, differences, phantoms, bandwidth,
% arrays, metrics. Each export stage writes one file.
% extend_unwindowed_control adds the paired original-bandwidth control to a completed cache.
out=fileparts(mfilename('fullpath'));root=fileparts(fileparts(out));
if nargin<2,cache=fullfile(tempdir,'cbct_reference_numerics');end
if nargin<3,destination=out;end
if ~isfolder(cache),mkdir(cache);end
file=fullfile(cache,'computed.mat');
if strcmp(stage,'extend_unwindowed_control'),addNativeControl(file,root,cache);return;end
if strcmp(stage,'compute')
 b=load(fullfile(root,'images/lesion_free_roundtrip_case_b_del_slice071/03_roundtrip_arrays.mat'));
 l=load(fullfile(root,'images/lesion_insertion_case_b_del_slice071/03_insertion_arrays.mat'));
 healthy=double(b.ground_truth);image_domain=double(l.image_domain);lesion_delta=l.lesion_delta;
 lesion_mask=l.lesion_mask;lesion_patch=l.lesion_patch;interior_mask=b.interior_mask;
 original_healthy=b.reconstruction;original_insertion=l.projection_domain;
 assert(isequal(healthy,l.healthy)&&isequal(healthy+lesion_delta,image_domain));
 assert(isequal(healthy(~lesion_mask),image_domain(~lesion_mask))&&nnz(lesion_mask)==246);
 assert(b.geometry.fov_mm==400 && b.geometry.views==2000 && b.geometry.detector_elements==2480);
 % Tests use the fixed professor geometry; all intensities remain unchanged.
 dx=400/512;referenceStep=1200/ceil((750+1200)/(dx/2));fineStep=referenceStep/2;
 selected=unique(round(linspace(1,2000,17)));
 [x,y]=meshgrid(((1:512)-256.5)*dx);central_mask=hypot(x,y)<150;
 masks=struct('all_pixels',true(512),'interior',interior_mask,'central150mm',central_mask,'lesion',lesion_mask);
 phantoms=phantomChecks(root,cache,selected,fineStep);
 pCoarse=project(healthy,referenceStep,4,1:2000,'healthy coarse');
 healthy_projections=project(healthy,fineStep,4,1:2000,'healthy fine');
 pQuarter=project(healthy,fineStep/2,4,selected,'quarter-step convergence');
 pMorePadding=project(healthy,fineStep/2,6,selected,'padding invariance');
 tests.projection_coarse_fine_relative_l2=relative(pCoarse,healthy_projections);
 tests.projection_fine_quarter_relative_l2=relative(healthy_projections(selected,:),pQuarter);
 tests.padding4_padding6_relative_l2=relative(pQuarter,pMorePadding);
 pInsufficient=project(healthy,fineStep/2,2,selected,'insufficient-padding control');
 tests.padding2_padding4_relative_l2=relative(pInsufficient,pQuarter);
 save(fullfile(cache,'quadrature_checkpoint.mat'),'pCoarse','healthy_projections','pQuarter','pMorePadding','tests','-v7.3');
 assert(tests.projection_coarse_fine_relative_l2<1e-5,'Quadrature not converged.');
 assert(tests.projection_fine_quarter_relative_l2<1e-5,'Finer selected-view quadrature not converged.');
 assert(tests.padding4_padding6_relative_l2<1e-5,'Padding extension not converged.');
 rCoarse=fbp(pCoarse,root,cache,'healthy_coarse',true);
 unwindowed_healthy=fbp(healthy_projections,root,cache,'healthy_unwindowed',false);
 healthy_reconstruction=fbp(healthy_projections,root,cache,'healthy_grid_bandwidth',true);
 tests.reconstruction_coarse_fine=stats(rCoarse-healthy_reconstruction,true(512));
 assert(tests.reconstruction_coarse_fine.rmse<.1,'Reconstruction quadrature tolerance exceeded.');
 interpolation_only=fbp(double(b.projections),root,cache,'interpolation_only',false);
 delta_projections=project(double(lesion_delta),fineStep,4,1:2000,'lesion delta');
 combined_projections=project(image_domain,fineStep,4,1:2000,'independent image-domain input');
 save(fullfile(cache,'independent_combined.mat'),'combined_projections','-v7.3');
 tests.projection_linearity_relative_l2=relative(healthy_projections+delta_projections,combined_projections);
 assert(tests.projection_linearity_relative_l2<1e-12,'Projection linearity failed.');
 projection_domain=fbp(healthy_projections+delta_projections,root,cache,'projection_insertion',true);
 matched_image_domain=fbp(combined_projections,root,cache,'image_insertion',true);
 reconstructed_delta=fbp(delta_projections,root,cache,'delta',true);
 tests.reconstruction_additivity=stats(projection_domain-healthy_reconstruction-reconstructed_delta,true(512));
 assert(tests.reconstruction_additivity.max_absolute<1e-6,'Backprojection additivity failed.');
 healthy_error=healthy_reconstruction-healthy;insertion_error=projection_domain-image_domain;
 matched_error=projection_domain-matched_image_domain;
 incremental_error=projection_domain-healthy_reconstruction-lesion_delta;
 names=fieldnames(masks);metrics=struct;
 for k=1:numel(names)
  key=names{k};mask=masks.(key);
  metrics.reference_healthy.(key)=stats(double(original_healthy)-healthy,mask);
  metrics.interpolation_only.(key)=stats(interpolation_only-healthy,mask);
  metrics.unwindowed_corrected.(key)=stats(unwindowed_healthy-healthy,mask);
  metrics.grid_bandwidth_healthy.(key)=stats(healthy_error,mask);
  metrics.reference_insertion.(key)=stats(double(original_insertion)-image_domain,mask);
  metrics.grid_bandwidth_insertion.(key)=stats(insertion_error,mask);
  metrics.matched_methods.(key)=stats(matched_error,mask);
  metrics.incremental_lesion.(key)=stats(incremental_error,mask);
 end
 assert(metrics.matched_methods.all_pixels.max_absolute<1e-6,'Paired-method agreement failed.');
 assert(metrics.grid_bandwidth_healthy.interior.rmse<metrics.reference_healthy.interior.rmse);
 metadata=struct('geometry',b.geometry,'case','Case B DEL slice 71','lesion','Case A S13_I137_DEL_19x18.raw', ...
 'units','Stored pixel units; raw decoding inferred; physical attenuation and HU calibration unavailable', ...
 'forward_model','Cubic convolution of the unchanged samples with zero samples outside the rectangular array; midpoint ray integration on its bounded support', ...
 'boundary','Four zero samples added outside each array edge; original coordinates and values unchanged; this is an explicit boundary-model change from hard extrapolation at the last pixel center', ...
 'backprojector','Professor equiangular fan FBP; normalized detector interpolation; complete filtered-convolution domain; double precision', ...
 'filter_domain','Zero extension of raw projections for linear convolution; filtered samples cover all reconstruction angles. Physical 2480 detector samples and their spacing are unchanged.', ...
 'reference_max_step_mm',referenceStep,'selected_max_step_mm',fineStep, ...
 'filter','Band-limited ramp with angular cutoff SID/(2*dx); same fan-beam sin(gamma) kernel correction', ...
 'angular_cutoff_cycles_per_radian',750/(2*dx),'isocenter_cutoff_cycles_per_mm',1/(2*dx), ...
 'filter_limit','Cutoff defined at isocenter; fan-beam spatial response varies with position. Resolution and ringing differ from the original detector-Nyquist ramp.', ...
 'not_applied','No baseline subtraction/restoration, fitted gain or offset, source-derived error-map subtraction, or post-reconstruction denoising', ...
 'display','Uniform stored value minus 8192, limits [-200 300], display only and not confirmed HU', ...
 'execution','MATLAB CPU; GPU path and clinical accuracy not validated','matlab_version',version, ...
 'reference_fbp_source',fileread(fullfile(root,'li_fbp/c2_fbp_v6_curved_detector.m')), ...
 'reference_forward_source',fileread(fullfile(root,'li_fbp/c1_forward_projection_v6_curved_detector.m')));
 save(file,'healthy','image_domain','lesion_delta','lesion_mask','lesion_patch','interior_mask','central_mask', ...
 'original_healthy','original_insertion','interpolation_only','unwindowed_healthy','healthy_reconstruction', ...
 'projection_domain','matched_image_domain','reconstructed_delta','healthy_error','insertion_error','matched_error', ...
 'incremental_error','healthy_projections','delta_projections','metrics','tests','phantoms','metadata','-v7.3');
 addNativeControl(file,root,cache);
 fprintf('PASS: independent phantom, convergence, zero-extension, linearity, and paired reconstruction checks.\n');
 disp(tests);disp(metrics.grid_bandwidth_healthy);disp(metrics.grid_bandwidth_insertion.lesion);
 return
end
s=load(file);assert(isequal(s.healthy_error,s.healthy_reconstruction-s.healthy));
assert(isequal(s.insertion_error,s.projection_domain-s.image_domain));
if ~isfolder(destination),mkdir(destination);end
switch stage
 case 'arrays'
  target=fullfile(destination,'07_arrays.mat');copyfile(file,target);
  q=load(target,'healthy_reconstruction','projection_domain','lesion_mask');
  assert(isequal(q.healthy_reconstruction,s.healthy_reconstruction)&&isequal(q.projection_domain,s.projection_domain)&&nnz(q.lesion_mask)==246);
 case 'metrics'
  metadata=rmfield(s.metadata,{'reference_fbp_source','reference_forward_source'});
  payload=struct('metadata',metadata,'metrics',s.metrics,'tests',s.tests,'phantom_metrics',s.phantoms.metrics);
  target=fullfile(destination,'06_metrics.json');fid=fopen(target,'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(payload,'PrettyPrint',true));fclose(fid);
  q=jsondecode(fileread(target));assert(q.metrics.grid_bandwidth_insertion.lesion.pixel_count==246);
 otherwise
  f=figure('Visible','off','Color','w','Position',[20 20 1600 1050]);cleanup=onCleanup(@()close(f));
  t=tiledlayout(f,2,3,'Padding','loose','TileSpacing','loose');
  title(t,'Direct cubic ray projection + fan-beam FBP | no baseline correction','FontSize',16);
  subtitle(t,{'Numerical boundary/interpolation corrections; ramp bandwidth set by output pixel spacing.', ...
   'Uniform display: stored value - 8192, window [-200, 300]. Assumed display offset only; not calibrated HU.'},'FontSize',11);
  switch stage
   case 'reconstruction'
    ct(s.healthy,'Input',false,false);ct(s.original_healthy,'BEFORE: original reference FBP',false,false);
    ct(s.healthy_reconstruction,'Updated FBP: grid-bandwidth result',false,false);
    ct(s.healthy,'Input: detail',true,false);ct(s.healthy_reconstruction,'Updated result: detail',true,false);
    difference(-s.healthy_error,50,'Input - reconstruction (full image)',false);
    target=fullfile(destination,'01_reconstruction.png');
   case 'insertion'
    ct(s.image_domain,'Direct image-domain replacement',false,true);
    ct(s.unwindowed_insertion,'Projection insertion: original bandwidth',false,true);
    ct(s.projection_domain,'Projection insertion: grid bandwidth',false,true);
    ct(s.image_domain,'Direct image replacement: detail',true,true);
    ct(s.unwindowed_insertion,'Original bandwidth: lesion detail',true,true);
    ct(s.projection_domain,'Grid bandwidth: lesion detail',true,true);
    target=fullfile(destination,'02_lesion_insertion.png');
   case 'differences'
    difference(double(s.original_healthy)-s.healthy,200,'Reference reconstruction error',false);
    difference(s.healthy_error,200,'Updated error: same scale',false);
    difference(s.healthy_error,50,'Updated error: expanded scale',false);
    difference(s.insertion_error,50,'Projection - direct image insertion',true);
    difference(s.incremental_error,50,'Incremental lesion error',true);
    difference(s.matched_error,max(1e-8,max(abs(s.matched_error(:)))),'Projection - reconstructed image insertion',true);
    target=fullfile(destination,'03_pixel_differences.png');
   case 'phantoms'
    title(t,'Independent analytic phantom validation');subtitle(t,'Synthetic dimensionless values; no patient-image calibration or fitted gain.');
    scalarPanel(s.phantoms.gaussian_target,'Analytic Gaussian target',[0 1]);
    scalarPanel(s.phantoms.gaussian_grid,'FBP of exact Gaussian projections',[0 1]);
    scalarPanel(s.phantoms.gaussian_grid-s.phantoms.gaussian_target,'Gaussian reconstruction error',[-.005 .005]);
    scalarPanel(s.phantoms.disk_target,'Analytic unit-disk target',[0 1]);
    scalarPanel(s.phantoms.disk_grid,'FBP of exact disk projections',[0 1]);
    ax=nexttile;plot(ax,s.phantoms.edge_x,s.phantoms.edge_native,'-',s.phantoms.edge_x,s.phantoms.edge_grid,'-',s.phantoms.edge_x,s.phantoms.edge_target,'--');
    xlim(ax,[100 110]);ylim(ax,[-.2 1.2]);grid(ax,'on');legend(ax,{'Detector-Nyquist ramp','Grid-Nyquist ramp','Ideal step'},'Location','southwest');
    xlabel(ax,'x (mm)');ylabel(ax,'Value');title(ax,'Edge resolution and ringing');
    target=fullfile(destination,'04_phantom_validation.png');
   case 'bandwidth'
    ct(s.interpolation_only,'FBP fixes; original projections',false,false);
    ct(s.unwindowed_healthy,'Boundary/integration fixes; original bandwidth',false,false);
    ct(s.healthy_reconstruction,'Same projections; grid-derived bandwidth',false,false);
    difference(s.interpolation_only-s.healthy,100,'FBP-only error',true);
    difference(s.unwindowed_healthy-s.healthy,100,'Original bandwidth error',true);
    difference(s.healthy_error,100,'Grid bandwidth error: same scale',true);
    target=fullfile(destination,'05_bandwidth_comparison.png');
   otherwise,error('Unknown stage: %s',stage);
  end
  exportgraphics(f,target,'Resolution',130);
end
fprintf('PASS: exported %s\n',target);
end

function addNativeControl(file,root,cache)
s=load(file);
independent=load(fullfile(cache,'independent_combined.mat'),'combined_projections');combined=independent.combined_projections;
assert(relative(s.healthy_projections+s.delta_projections,combined)<1e-12);
s.unwindowed_insertion=fbp(s.healthy_projections+s.delta_projections,root,cache,'insertion_native',false);
s.unwindowed_matched_image=fbp(combined,root,cache,'matched_native',false);
s.unwindowed_delta=fbp(s.delta_projections,root,cache,'delta_native',false);
s.unwindowed_insertion_error=s.unwindowed_insertion-s.image_domain;
s.unwindowed_incremental_error=s.unwindowed_insertion-s.unwindowed_healthy-s.lesion_delta;
s.unwindowed_matched_error=s.unwindowed_insertion-s.unwindowed_matched_image;
masks=struct('all_pixels',true(512),'interior',s.interior_mask,'central150mm',s.central_mask,'lesion',s.lesion_mask);
names=fieldnames(masks);
for k=1:numel(names)
 key=names{k};mask=masks.(key);
 s.metrics.unwindowed_insertion.(key)=stats(s.unwindowed_insertion_error,mask);
 s.metrics.unwindowed_incremental_lesion.(key)=stats(s.unwindowed_incremental_error,mask);
 s.metrics.unwindowed_matched_methods.(key)=stats(s.unwindowed_matched_error,mask);
end
s.tests.unwindowed_additivity=stats(s.unwindowed_insertion-s.unwindowed_healthy-s.unwindowed_delta,true(512));
assert(s.tests.unwindowed_additivity.max_absolute<1e-6);
assert(s.metrics.unwindowed_matched_methods.all_pixels.max_absolute<1e-6);
s.metadata.filter_selection='Both bandwidths retained. Grid-derived cutoff reduces background aliasing but can increase lesion-edge error; it is not an unqualified improvement in lesion fidelity.';
save(file,'-struct','s','-v7.3');
fprintf('Original bandwidth lesion RMSE %.8f; incremental lesion RMSE %.8f\n',s.metrics.unwindowed_insertion.lesion.rmse,s.metrics.unwindowed_incremental_lesion.lesion.rmse);
fprintf('Grid bandwidth lesion RMSE %.8f; incremental lesion RMSE %.8f\n',s.metrics.grid_bandwidth_insertion.lesion.rmse,s.metrics.incremental_lesion.lesion.rmse);
end

function p=project(input,step,padding,viewIndices,label)
% Cubic interpolation preserves the original samples; padding is exterior only.
assert(isequal(size(input),[512 512])&&all(isfinite(input(:))));
original=input;h=padarray(double(input),[padding padding],0,'both');
assert(isequal(h(padding+(1:512),padding+(1:512)),original));
M=512;dx=400/M;nu=2480;sid=750;sdd=1200;dg=.35/sdd;
gamma=((1:nu)'-(nu/2+.5))*dg;lo=-(M/2-.5+padding)*dx;hi=-lo;p=zeros(numel(viewIndices),nu);
for j=1:numel(viewIndices)
 v=viewIndices(j);theta=(v-1)*2*pi/2000;xs=sid*cos(theta);ys=sid*sin(theta);
 ux=-cos(theta+gamma);uy=-sin(theta+gamma);
 [entry,span]=rayInterval(xs,ys,ux,uy,lo,hi,sdd);active=span>0;
 if ~any(active),continue;end
 N=ceil(max(span)/step);t=(.5:N-.5)/N;ell=entry(active)+span(active).*t;
 ix=(xs+ux(active).*ell)/dx+M/2+.5+padding;iy=(ys+uy(active).*ell)/dx+M/2+.5+padding;
 values=interp2(h,ix,iy,'cubic',0);p(j,active)=sum(values,2).*span(active)/N;
 if mod(j,500)==0,fprintf('%s: view %d/%d\n',label,j,numel(viewIndices));end
end
assert(all(isfinite(p(:))));
end

function [entry,span]=rayInterval(xs,ys,ux,uy,lo,hi,sdd)
% Slab intersection with the interpolation domain, not a segmented image mask.
a=(lo-xs)./ux;b=(hi-xs)./ux;c=(lo-ys)./uy;d=(hi-ys)./uy;
entry=max(max(min(a,b),min(c,d)),0);leave=min(min(max(a,b),max(c,d)),sdd);
span=max(leave-entry,0);
end

function ph=phantomChecks(root,cache,selected,step)
[x,y]=meshgrid(((1:512)-256.5)*400/512);gamma=((1:2480)-(1240+.5))*(.35/1200);
theta=(0:1999)'*2*pi/2000;xs=750*cos(theta);ys=750*sin(theta);ux=-cos(theta+gamma);uy=-sin(theta+gamma);
cx=32;cy=-21;sx=20;sy=30;ph.gaussian_target=exp(-.5*(((x-cx)/sx).^2+((y-cy)/sy).^2));
a=(ux/sx).^2+(uy/sy).^2;b=(xs-cx).*ux/sx^2+(ys-cy).*uy/sy^2;c=((xs-cx)/sx).^2+((ys-cy)/sy).^2;
analytic=sqrt(2*pi./a).*exp(-.5*(c-b.^2./a));
numerical=project(ph.gaussian_target,step,4,selected,'Gaussian selected rays');
ph.metrics.gaussian_projection_relative_l2=relative(numerical,analytic(selected,:));
assert(ph.metrics.gaussian_projection_relative_l2<1e-4);
ph.gaussian_native=fbp(analytic,root,cache,'gaussian_native',false);
ph.gaussian_grid=fbp(analytic,root,cache,'gaussian_grid',true);
ph.metrics.gaussian_native=stats(ph.gaussian_native-ph.gaussian_target,true(512));
ph.metrics.gaussian_grid=stats(ph.gaussian_grid-ph.gaussian_target,true(512));
ph.metrics.gaussian_grid_relative_l2=relative(ph.gaussian_grid,ph.gaussian_target);
assert(ph.metrics.gaussian_grid_relative_l2<.005);
cx=25;cy=-17;radius=80;distance=abs((cx-xs).*uy-(cy-ys).*ux);analytic=2*sqrt(max(radius^2-distance.^2,0));
ph.disk_target=double(hypot(x-cx,y-cy)<=radius);
ph.disk_native=fbp(analytic,root,cache,'disk_native',false);ph.disk_grid=fbp(analytic,root,cache,'disk_grid',true);
inner=hypot(x-cx,y-cy)<radius-3*(400/512);
ph.metrics.disk_native_interior=stats(ph.disk_native-ph.disk_target,inner);
ph.metrics.disk_grid_interior=stats(ph.disk_grid-ph.disk_target,inner);
ph.metrics.disk_grid_all=stats(ph.disk_grid-ph.disk_target,true(512));
assert(ph.metrics.disk_grid_interior.rmse<.01);
[~,row]=min(abs(y(:,1)-cy));ph.edge_x=x(row,:);ph.edge_native=ph.disk_native(row,:);ph.edge_grid=ph.disk_grid(row,:);ph.edge_target=ph.disk_target(row,:);
ph.metrics.edge_native_10_90_mm=edgeWidth(ph.edge_x,ph.edge_native,cx+radius);
ph.metrics.edge_grid_10_90_mm=edgeWidth(ph.edge_x,ph.edge_grid,cx+radius);
% Constant-valued square: analytic line lengths independently test quadrature.
angles=[1 126 251 387];numerical=project(ones(512),step,0,angles,'constant square');exact=zeros(size(numerical));old=exact;
dx=400/512;lo=-255.5*dx;N=ceil((750+1200)/(dx/2));ds=1200/N;t=(.5:N-.5)*ds;
for j=1:numel(angles)
 v=angles(j);[~,len]=rayInterval(xs(v),ys(v),ux(v,:)',uy(v,:)',lo,-lo,1200);exact(j,:)=len;
 qx=xs(v)+ux(v,:)'.*t;qy=ys(v)+uy(v,:)'.*t;old(j,:)=sum(qx>=lo&qx<=-lo&qy>=lo&qy<=-lo,2)*ds;
end
ph.metrics.square_original=stats(old-exact,true(size(exact)));
ph.metrics.square_bounded=stats(numerical-exact,true(size(exact)));
assert(ph.metrics.square_bounded.max_absolute<1e-9);
end

function width=edgeWidth(x,p,center)
levels=[.9 .1];cross=zeros(1,2);
for k=1:2
 ids=find(p(1:end-1)>=levels(k)&p(2:end)<levels(k)&abs(x(1:end-1)-center)<4);
 assert(~isempty(ids));[~,j]=min(abs(x(ids)-center));i=ids(j);
 cross(k)=x(i)+(levels(k)-p(i))*(x(i+1)-x(i))/(p(i+1)-p(i));
end
width=cross(2)-cross(1);
end

function v=relative(a,b)
v=norm(double(a(:))-double(b(:)))/max(norm(double(b(:))),eps);
end
function m=stats(error,mask)
e=double(error(mask));assert(~isempty(e)&&all(isfinite(e)));
m=struct('pixel_count',numel(e),'mae',mean(abs(e)),'rmse',sqrt(mean(e.^2)),'bias',mean(e),'max_absolute',max(abs(e)));
end
function scalarPanel(a,label,limits)
ax=nexttile;imagesc(ax,a',limits);axis(ax,'image');axis(ax,'off');colormap(ax,gray);colorbar(ax);title(ax,label);
end

function r=fbp(inputProjections,root,cache,tag,gridBandwidth)
physical=inputProjections;assert(isequal(size(physical),[2000 2480]));
dg=.35/1200;extent=sqrt(2)*(512/2-.5)*(400/512);
padding=ceil((asin(extent/750)-(2480/2-.5)*dg)/dg)+2;
inputProjections=padarray(physical,[0 padding],0,'both');
assert(isequal(inputProjections(:,padding+(1:2480)),physical));
physicalGamma=((1:2480)-1240.5)*dg;
filterGamma=((1:size(inputProjections,2))-(size(inputProjections,2)/2+.5))*dg;
assert(isequal(filterGamma(padding+(1:2480)),physicalGamma));
liInput=fullfile(cache,[tag '_projections.raw']);liOutput=fullfile(cache,[tag '_reconstruction.raw']);
fid=fopen(liInput,'w','l');assert(fid>=0);fwrite(fid,inputProjections','double');fclose(fid);
code=fileread(fullfile(root,'li_fbp/c2_fbp_v6_curved_detector.m'));
assert(numel(strfind(code,'detector_nums = 2480;'))==1);
code=strrep(code,'detector_nums = 2480;',sprintf('detector_nums = %d;',size(inputProjections,2)));
code=strrep(code,'gpuDevice(1);','');
code=strrep(code,'''projections.raw''',['''' liInput '''']);code=strrep(code,'''img.raw''',['''' liOutput '''']);
code=strrep(code,'gpuArray.zeros(','zeros(');code=strrep(code,', ''gpuArray''','');
code=strrep(code,'gpuArray(','(');code=strrep(code,'gather(','(');
code=regexprep(code,'(?m)^(figure;|imagesc\([^\r\n]*|colormap gray;[^\r\n]*|title\([^\r\n]*)','');
code=strrep(code,"fprintf('View %d / %d\n', idx_angle, views);",'');
if true
    old='w = gamma0 - gamma(k);';assert(numel(strfind(code,old))==1);
    code=strrep(code,old,'w = (gamma0 - gamma(k)) ./ d_gamma;');
    old='img = img + SID ./ (L.^2) .* val * d_theta;';assert(numel(strfind(code,old))==1);
    code=strrep(code,old,sprintf('val(~valid) = 0;\n    %s',old));
end
code=strrep(code,'''single''','''double''');code=strrep(code,'single(','double(');
if gridBandwidth
code=char(code);start=strfind(code,'% --- build ramp kernel');stop=strfind(code,'% --- weighting ---');
assert(numel(start)==1 && numel(stop)==1);
replacement=sprintf(['n=-(detector_nums-1):(detector_nums-1); g=n*d_gamma; cutoff=SID/(2*dx);\n' ...
'ramp_kernel=zeros(size(g)); nz=n~=0;\n' ...
'ramp_kernel(nz)=(cutoff*g(nz).*sin(2*pi*cutoff*g(nz))/pi+(cos(2*pi*cutoff*g(nz))-1)/(2*pi^2))./sin(g(nz)).^2;\n' ...
'ramp_kernel(~nz)=cutoff^2;\n']);
code=[code(1:start-1) replacement code(stop:end)];
end
sdd=double(1200); % Preserve the original forward-script workspace value.
eval(code);r=img;
assert(isequal(projection_data,inputProjections) && all(isfinite(r(:))));
fprintf('Backprojection complete: %s\n',tag);
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
