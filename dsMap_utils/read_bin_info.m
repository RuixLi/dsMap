function info = read_bin_info(fname)
% read binary file info
ubyte = 2; % 16 bit
pixelSize = [512,512];
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