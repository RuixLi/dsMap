function generate_color_wheel(nDir,saveDir)
% generate 2 color wheels
% indicate directions and orientations in HSV pseudocolor map

rInner = 80;     % inner radius of the colour ring
rOuter = 200;    % outer radius of the colour ring
[x, y] = meshgrid(-rOuter:rOuter);
[theta, rho] = cart2pol(x, y);
dirhue = wrapTo2Pi(-theta) / (2 * pi);
dirhue = floor(dirhue * nDir) / nDir;
dirhue=hue_KO2HSV(dirhue);
dirsatu = ones(size(dirhue)); 
dirval = double(rho >= rInner & rho <= rOuter); 
dirColorWheel = hsv2rgb(cat(3, dirhue, dirsatu, dirval));
imwrite(dirColorWheel,fullfile(saveDir,'_dirColorWheel.tif'));
orihue = wrapTo2Pi(-theta);
orihue(orihue>pi) = orihue(orihue>pi)-pi;
orihue = orihue / pi;
orihue = floor(orihue * (nDir/2)) / (nDir/2);
orihue=hue_KO2HSV(orihue);
orisatu = ones(size(orihue)); 
orival = double(rho >= rInner & rho <= rOuter); 
oriColorWheel = hsv2rgb(cat(3, orihue, orisatu, orival));
imwrite(oriColorWheel,fullfile(saveDir,'_oriColorWheel.tif'));