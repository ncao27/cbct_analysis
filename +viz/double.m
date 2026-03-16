function [] = double(img)
    figure;
    imagesc(img);
    colormap gray;
    colorbar;
    title('Check the background value!');
end