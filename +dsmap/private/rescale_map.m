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

nel = min(1e12,numel(dataIN));
inMin=quantile(dataIN(1:nel),0.01);
inMax=quantile(dataIN(1:nel),0.99);

dataIN(dataIN<inMin) = inMin;
dataIN(dataIN>inMax) = inMax;

scaleFactor = (outRange(2)-outRange(1)) / (inMax-inMin);
dataOUT=(dataIN - inMin).* scaleFactor + outRange(1);

dataOUT(dataOUT > outRange(2)) = outRange(2);
dataOUT(dataOUT < outRange(1)) = outRange(1);

end
