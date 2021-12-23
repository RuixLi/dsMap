function data = load_bin_data(fname,idx,fsiz)
% load saved binary file, only support uint16 type
% if the data is not saved as uint16, the result will get wrong

% INPUT
% fname[string] full path to the data
% idx [start, end] indicate the start and end point of read
% fsiz, [d1, d2] indicate the size of imaging data

if nargin < 2; idx = [1,-1]; end
if nargin < 3; fsiz = [512,512]; end

fid = fopen(fname, 'r');

if idx(2) == -1
    ubyte = 2; % 16 bit
    pixelSize = [512,512];
    frameChunkSize = pixelSize(1)*pixelSize(2)*ubyte;
    currentSeek = ftell(fid);
    fseek(fid, 0, 1);
    fileLength = ftell(fid);
    fseek(fid, currentSeek, -1);
    idx(2) = fileLength/frameChunkSize;
end

n = idx(2) - idx(1) + 1;
fseek(fid,(idx(1)-1)*fsiz(1)*fsiz(2)*2,-1);

data = fread(fid,  fsiz(1)*fsiz(2)*n, 'uint16=>uint16');
data = reshape(data,fsiz(1),fsiz(2),[]);
fclose(fid);
end