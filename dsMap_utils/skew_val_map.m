function skMap = skew_val_map(Y)
% project skewness along temproal axis of xyt image data into a 2d map
% assume Y has [Y,X,T] shape

%maxV = rescale(max(Y,[],3));
[d1,d2,T] = size(Y);
d = d1*d2;
Y = reshape(Y,d,T);

mpV = skewness(Y,0,2);
map = reshape(mpV,d1,d2);
%map = maxV.*map;

if nargout < 1
    figure; im_view(map); 
else
    skMap = map;
end
