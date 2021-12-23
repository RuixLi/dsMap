function sdMap = sd_val_map(Y)
% project std along temproal axis of xyt image data into a 2d map
% assume Y has [Y,X,T] shape

[d1,d2,T] = size(Y);
d = d1*d2;
Yr = reshape(Y,d,T);

mpV = std(Yr,0,2);
map = reshape(mpV,d1,d2);

if nargout < 1
    figure; imshow(rescale(map)); 
else
    sdMap = map;
end