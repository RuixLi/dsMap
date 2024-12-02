function [Avg,trialType,trialIndL,trialIndS] = sort_and_avg_rand_trial(data, trialLength, trialList)

% trialAvg average data across different trial types randomly arranged

% INPUT
% data[2D or 3D array] data with different trials mixed
% - if data is 2D assume 1st dim contains different trials
% - eg. data is a T x K, contains N x C types of trial along time bin T
% - N is repeat time of each trial type, C is number of trial type
% - K is number of ROIs or neurons, T = N x C x S
% - if data is 3D assume 3st dim contains different trials
% - eg. data is a XYT imaging data, contains different trials along T
% trialLength[scalar], length S of a single trial
% trialList[list] sequence of trial
% - the length can have length of N x C
% - eg. [0 1 2 3 2 1 3 0 1 3 2 0]
% - C is 4, N is 3, S is unkonwn
% - the length can also be the same with T
% - eg. [0 0 0 4 4 4 1 1 1 4 4 4 1 1 1 0 0 0]
% - S is 5, C is 3, N is 2

% OUTPUT
% Avg [2+1D or 3+1D array] the last dim is average of each type
% trialType[C x 1 list]
% trialIndL[T x C binary] indicate the location of all trials
% - T is the length of data, C is the type number of trials
% trialIndS[N x C binary] indicate the location of all trials
% - N is the number of each trial type, C is the number of trial type


% written by Ruix.Li in May 2021

%% get trial type category
dataType = class(data);
trialType = str2double(categories(categorical(trialList)));
trialTypeNum = length(trialType);
%fprintf('find %d types of trial from the sequence list\n', trialTypeNum)
%fprintf('check data completeness ...\n')

ndFlg = 0;

if ismatrix(data)
    dataLength = size(data,1);
    %fprintf('assume 1st dim contains different trials\n')
    K = size(data,2);
    ndFlg = 2;
end

if ndims(data) == 3
    dataLength = size(data,3);
    %fprintf('assume 3st dim contains different trials\n')
    ndFlg = 3;
end

if ~ndFlg; error('only support 2D and 3D array data'); end

if dataLength ~= (length(trialList) * trialLength) && dataLength ~= length(trialList)
    error('trial length and sequence do not match data length')
end

%% get longSequence (the same length as data) from sequence or vice versa
if dataLength == length(trialList)
    longSequence = trialList(:);
    s =  reshape(longSequence,trialLength,[]);
    trialList = s(1,:)';
else
    trialList = trialList(:);
    longSequence = reshape(repmat(trialList,1,trialLength)',dataLength,1);
end

%% get logical index
trialIndL = zeros(dataLength,trialTypeNum);
trialIndS = zeros(length(trialList),trialTypeNum);

for i = 1:trialTypeNum
    trialIndL(:,i) = ismember(longSequence,trialType(i));
    trialIndS(:,i) = ismember(trialList,trialType(i));
end

trialIndS = logical(trialIndS);
trialIndL = logical(trialIndL);

%disp('average trials ...')

if ndFlg == 2
    data = reshape(data',K,trialLength,[]);
    data = permute(data,[2,1,3]);
    Avg = zeros(trialLength, K, trialTypeNum);
    for i = 1:trialTypeNum
        Avg(:,:,i) = mean(data(:,:,trialIndS(:,i)),3);
    end
end

if ndFlg == 3
    Avg = zeros(size(data,1), size(data,2), trialLength, trialTypeNum);
    for i = 1:trialTypeNum
        Avg(:,:,:,i) = mean(reshape(data(:,:,trialIndL(:,i)),size(data,1),size(data,2),trialLength,[]),4);
    end
end

Avg = cast(Avg,dataType);

end