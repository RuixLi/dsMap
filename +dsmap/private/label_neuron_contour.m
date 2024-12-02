function map = label_neuron_contour(FP, FOV, indx)
% label FP on FOV with contour

% INPUTS:
% FP[P x K matrix], each column is the spatial footprint of a element
% FOV, an image
% indx, the order of FP, optional

% OUTPUT:
% map: an RGB image of labeled FP, ranged between [0, 1]
if ~exist('FOV','var')
    fovsiz = sqrt(size(FP,1));
    FOV = zeros(fovsiz,fovsiz);
end

FP = full(FP);
FOV = .7*rescale(FOV);
im = cat(3,FOV,FOV,FOV);
[d1,d2,~] = size(im);
im = reshape(im,d1*d2,3);
k = size(FP,2);
if ~exist('indx','var'); indx = 1:k; end
se = strel('disk',1);

for i = 1:k
    x = double(reshape(FP(:,indx(i))>0,d1,d2));
    y = sign(imdilate(x,se));
    pv = find(y-x);
    im(pv,1) = 1;
    im(pv,2) = 0;
    im(pv,3) = 0;
end

im = reshape(im,512,512,3);

if nargout
    map = im;
else
    figure; imshow(im);
end

end