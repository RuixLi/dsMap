function h = im_view(im)
% im_view visualize im
% im can be a vector, imv reshape im into square matrix
% im can be a matrix, imv treat im as gray-scale image
% im can be a X-T-3 tensor, imv treat im as RGB image 
% im can be a X-Y-T tensor, imv average im along T
% im can be a cell (future plan)

if size(im,1) == 1 && sqrt(size(im,2)) == round(sqrt(size(im,2)))
   d = sqrt(size(im,2));
   if d ~= round(d); return; end
   im = reshape(im,d,d);
end

if size(im,2) == 1 && sqrt(size(im,1)) == round(sqrt(size(im,1)))
   d = sqrt(size(im,1));
   im = reshape(im,d,d);
end

if ismatrix(im)
        clim = [quantile(im(:),0.01),quantile(im(:),0.99)];
        if clim(2) < clim(1); t = clim(1); clim(1) = clim(2); clim(2) = t; end
        if clim(2) == clim(1); clim(2) = clim(1) + .1; end
        imH = imagesc(im,clim);
        axis tight; axis equal; axis off;
        colormap gray
        colorbar
    
elseif size(im,3) == 3
    imH = imshow(im);
else
    im = mean(im,3);
    clim = [quantile(im(:),0.01),quantile(im(:),0.99)];
    if clim(2) < clim(1); t = clim(1); clim(1) = clim(2); clim(2) = t; end
    if clim(2) == clim(1); clim(2) = clim(1) + .1; end
    imH = imagesc(im,clim);
    axis tight; axis equal; axis off;
    colormap gray
    colorbar
    
if nargout >= 1; h = imH; end

end