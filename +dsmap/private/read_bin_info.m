function info = read_bin_info(fname,pixelSize)
% read binary file info

if nargin < 2
    pixelSize = [512,512];
    warning('data size of bin file not provided, assuming 512x512');
end

if numel(pixelSize) < 2
    pixelSize = [512,512];
    warning('data size of bin file not provided, assuming 512x512');
end

ubyte = 2; % 16 bit
frameChunkSize = pixelSize(1)*pixelSize(2)*ubyte;
fid = fopen(fname);
currentSeek = ftell(fid);
fseek(fid, 0, 1);
fileLength = ftell(fid);
fseek(fid, currentSeek, -1);
nFrame = fileLength/frameChunkSize;
info.Loops = nFrame;
fclose(fid);
end