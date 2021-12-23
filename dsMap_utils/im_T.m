function [dataout, info] = im_T(datain,mask)

% transform the imaging data array form XYT to TP or from TP to XYT

% INPUT
% datain[XYT/TP], data to transfrom
% mask[struct], mask struct (optional)
% -if mask is specific, only transform pixels where mask == 1

% OUTPUT
% dataout[TP/XYT], transformed data
% info[struct], transform information [d1,d2,(pixN),T]

% updated by RuiX.Li from IMRESH in Jan, 2021

if nargin == 1
% no mask
if length(size(datain)) == 3
    [d1,d2,T] = size(datain);
    dataout = permute(reshape(datain,[d1*d2,T]),[2,1]);
    
    info = [d1,d2,T];
    
elseif length(size(datain)) == 2
    [T,d] = size(datain);
    dd = sqrt(d);
    if floor(dd) == dd
        dataout = reshape(datain',[dd,dd,T]);
        info = [dd,dd,T];
        
    else
        error('need mask for correct transform')
    end
else
    error('wrong data dimentionality')
end

end

if nargin == 2
% get mask
if length(size(datain)) == 3
    [d1,d2,T] = size(datain);
    pixN = length(mask.maskind);
    data = permute(reshape(datain,[d1*d2,T]),[2,1]);
    dataout = data(:,mask.maskind);
    
    info = [d1,d2,pixN,T];
    
elseif length(size(datain)) == 2
    [T,pixN] = size(datain);
    [d1,d2] = size(mask.bitmap);
    dataout = zeros(T,d1*d2);
    dataout(:,mask.maskind) = datain;
    dataout = reshape(dataout',[d1,d2,T]);
    
    info = [d1,d2,pixN,T];
else
    error('wrong data dimentionality')
end

end
end
