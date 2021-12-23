function hue2=hue_KO2HSV(hue1)
% convert hue values in KO color space to hue vaules in HSV color space
% adapt from hueKO2HSV by Kenichi Ohki  2008.11.4
hue1(hue1<0)=0;
hue1(hue1>1)=1;
hue1=hue1*2/3;
hue1(hue1>(1/3))=hue1(hue1>(1/3))*2-1/3;
hue2=hue1;
end