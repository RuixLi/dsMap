function mpMap = max_val_map(Y)
% project maxium along temproal axis of XYT image data into a 2d map
% assume Y has [Y,X,T] shape

map = max(Y,[],3);

if nargout < 1
    figure; imview(map); 
else
    mpMap = map;
end

end