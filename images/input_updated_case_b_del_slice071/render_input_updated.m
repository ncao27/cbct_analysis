function render_input_updated(stage,destination)
% Export two-panel input/updated comparisons from existing validated arrays.
% No reconstruction is recomputed and no numerical source file is modified.
% Stages: healthy_reference, healthy_grid, lesion_reference, lesion_grid.
out=fileparts(mfilename('fullpath'));
if nargin<2,destination=out;end
source=fullfile(fileparts(out),'direct_reference_case_b_del_slice071','07_arrays.mat');
switch stage
 case 'healthy_reference'
  inputName='healthy';updatedName='unwindowed_healthy';
  heading='Healthy image | reference ramp bandwidth';
  filename='01_healthy_reference_bandwidth.png';
 case 'healthy_grid'
  inputName='healthy';updatedName='healthy_reconstruction';
  heading='Healthy image | grid-derived ramp bandwidth';
  filename='02_healthy_grid_bandwidth.png';
 case 'lesion_reference'
  inputName='image_domain';updatedName='unwindowed_insertion';
  heading='Lesion insertion | reference ramp bandwidth';
  filename='03_lesion_reference_bandwidth.png';
 case 'lesion_grid'
  inputName='image_domain';updatedName='projection_domain';
  heading='Lesion insertion | grid-derived ramp bandwidth';
  filename='04_lesion_grid_bandwidth.png';
 otherwise
  error('Unknown comparison stage: %s',stage);
end
s=load(source,inputName,updatedName);
input=s.(inputName);updated=s.(updatedName);
assert(isequal(size(input),[512 512])&&isequal(size(updated),[512 512]));
assert(all(isfinite(input(:)))&&all(isfinite(updated(:))));
isLesion=startsWith(stage,'lesion_');
if isLesion
 labels={'Input: image-domain lesion replacement', ...
         'Updated reconstruction: projection-domain insertion'};
else
 labels={'Input image','Updated reconstruction'};
end
f=figure('Visible','off','Color','w','Position',[20 20 1400 760]);
cleanup=onCleanup(@()close(f));
t=tiledlayout(f,1,2,'Padding','loose','TileSpacing','loose');
title(t,heading,'FontSize',17,'FontWeight','bold');
subtitle(t,'Display: stored value - 8192; window [-200, 300]. Assumed offset, not calibrated HU.', ...
 'FontSize',11);
values={input,updated};
for k=1:2
 ax=nexttile(t);
 shown=double(values{k})'-8192;
 im=imagesc(ax,shown,[-200 300]);
 axis(ax,'image');axis(ax,'off');colormap(ax,gray(256));
 title(ax,labels{k},'FontSize',14);
 c=colorbar(ax);c.Label.String='Stored value - 8192 (assumed)';
 if isLesion
  hold(ax,'on');
  rectangle(ax,'Position',[165.5 250.5 19 18],'EdgeColor',[1 .3 .1],'LineWidth',1.1);
 end
 assert(isequal(im.CData,shown)&&isequal(ax.CLim,[-200 300]));
end
assert(numel(findall(f,'Type','image'))==2);
if ~isfolder(destination),mkdir(destination);end
target=fullfile(destination,filename);
exportgraphics(f,target,'Resolution',180);
info=imfinfo(target);assert(info.Width>1000&&info.Height>500);
fprintf('PASS: %s and %s rendered with matching display scales to %s\n',inputName,updatedName,target);
end
