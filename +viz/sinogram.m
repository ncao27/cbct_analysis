function [] = sinogram(projections)
    imshow(projections,[])
    colormap gray
    xlabel('View angle')
    ylabel('Detector pixel')
end