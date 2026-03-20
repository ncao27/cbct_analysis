function [] = double(img)
    figure;
    % We manually force the color limits (clim) from 0 to roughly the max 
    % value of the soft tissue in the original image.
    imagesc(img, [0, max(img(:)) * 0.6]);
    colormap gray;
    %caxis([0.015 0.04]) this is for better visualization with orig imgs
    colorbar;
    title('Reconstructed Image (With Lesion Inserted)');
end