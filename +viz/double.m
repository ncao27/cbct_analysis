function [] = double(img)
    figure;
    % We manually force the color limits (clim) from 0 to roughly the max 
    % value of the soft tissue in the original image.
    imagesc(img, [0, max(img(:)) * 0.6]);
    colormap gray;
    colorbar;
    title('');
end