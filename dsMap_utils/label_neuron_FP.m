function imL = label_neuron_FP(FP, FOV, cm)

% lable FP on FOV with specific colormap
% 
FOV = .7*rescale(FOV);
im = cat(3,FOV,FOV,FOV);
[d1,d2,~] = size(im);
im = reshape(im,d1*d2,3);
if ~exist('cm','var'); cm = lines; end
n = size(cm,1);
for i = 1:size(FP,2)
pidx = FP(:,i)>0;
pvalue = rescale_map(FP(pidx,i),'n');
im(pidx,:) = bsxfun(@times, cm(mod(i,n)+1,:),pvalue);
end
im = reshape(im,d1,d2,3);

if nargout == 0
imagesc(im)
title('ROIs location')
axis tight; axis off; axis equal;
else
imL = im;
end

end