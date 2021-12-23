function dataOUT = rescale_map(dataIN, outRange)
% rescale the input map to outRange

% INPUT
% dataIN, input data, rescale will ingnore the structure
% outRange, [min, max] indicate range after rescale, deflut [0,1]
% if outRange = 'z', zscore the data
% if outRange = 's', sigmoid rescale the data
% - sigmoid normalization is defined by Y = 1 / (1 + exp(-X));
% if outRange = 'n', non-linearly rescale the range into [0,1]
% - this is better if the data contains extreme value

% OUTPUT
% dataOUT, data after rescale

% wirtten by Ruix.Li at Apr,2019
% modified in Nov, 2019, upgrade the function
% modified in Feb, 2020, add the zscore function
% modified in Jun, 2020, fix NaN bug
% modified in Jul, 2021, add non-linear rescale

if nargin < 2
    outRange = [0,1]; % default range is [0,1]
end

dataIN = double(dataIN);

if strcmpi(outRange,'z')
    dataIN = dataIN-mean(dataIN(:));
    dataST = std(dataIN(:),0,1);
    if isequal(dataST,0)
        warning('zscore error, std = 0')
        dataOUT = dataIN;
    else
        dataOUT = dataIN / dataST;
    end
    return
end

if strcmpi(outRange,'s')
    dataOUT = sigmoidn(dataIN);
    return
end

if strcmpi(outRange,'n')
    outRange = [0, 1];
    nel = min(1e12,numel(dataIN));
    inMin=quantile(dataIN(1:nel),0.01);
    inMax=quantile(dataIN(1:nel),0.99);
    
else
    inMin=min(dataIN(:));
    inMax=max(dataIN(:));
end

if inMin == inMax
    dataOUT = zeros(size(dataIN));
else
    scaleFactor = (outRange(2)-outRange(1)) / (inMax-inMin);
    dataOUT=(dataIN - inMin).* scaleFactor + outRange(1);
end

dataOUT(dataOUT > outRange(2)) = outRange(2);
dataOUT(dataOUT < outRange(1)) = outRange(1);

end

function Y = sigmoidn(X)
% sigmd returen sigmoid normalized Y from X
% The sigmoid operation is defined by
% Y = 1 / (1 + exp(-X));
Y = zeros(size(X),class(X));
for i = 1:numel(X); Y(i) = 1/(1+exp(-X(i))); end
end