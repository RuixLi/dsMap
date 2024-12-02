function stack = load_raw_data(fname,indx,nCH,pixelSize)
% load thorlabs raw file, much faster IO speed than mat file
% specific the index with rdx
% stack is a multi-dimentional array
%	1st Dimension - Image Y axis.
%	2nd Dimension - Image X axis.
%	3rd Dimension - time point
%	4th Dimension - number of channels.

% INPUT
% fname
% indx, [startIndex, endIndex] or startIndex only
% nCH, channel number
% pixelSize, [nRow, nCol]

if nargin < 4 || isempty(pixelSize)
    pixelSize = [512,512];
end

if nargin < 3 || isempty(nCH)
    nCH = 1;
end

if nargin < 2 || isempty(indx)
    indx = [1,-1];
end

if isscalar(indx); indx = [indx -1]; end

ubyte = 2; % 16 bit
imSize = pixelSize(1)*pixelSize(2);
frameChunkSize = pixelSize(1)*pixelSize(2)*ubyte;
fid = fopen(fname);
currentSeek = ftell(fid);
fseek(fid, 0, 1);
fileLength = ftell(fid);
fseek(fid, currentSeek, -1);
nFrame = fileLength/frameChunkSize;

if indx(2) == -1
    indx(2) = nFrame;
    indx(1) = nCH*indx(1);
else
    indx = nCH*indx;
end

if indx(2) > nFrame
    warning('give index more than exist')
end

if nFrame ~= round(nFrame); error('wrong pixelSize'); end

nRead = indx(2) - indx(1) + 1;
fprintf('load %d frames ... ', nRead)

fseek(fid,(indx(1)-1)*frameChunkSize,-1);
tic;

stack = fread(fid,imSize*nRead,'uint16=>uint16');

if nCH == 1
    stack = reshape(stack,[pixelSize(1),pixelSize(2),nRead]);
    stack = permute(stack,[2,1,3]);
else
    stack = reshape(stack,[pixelSize(1),pixelSize(2),nCH,nRead/nCH]);
    stack = permute(stack,[2,1,4,3]);
end

fclose(fid);
fprintf('done in %2.2f\n', toc);

end