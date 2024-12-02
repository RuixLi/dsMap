function [im_ch1, frametime] = load_nd2_data(file)
% A simple functon to read single channel image stack from nd2 file.
% Currently works only for 1 channel (should be easy to modify to multi
% channel data).
% written by T. Matsui (2024/5/25)
%
% Modified and simplified from nd2finfo.m and nd2read.m by Chao-Yuan Yeh
%
% Key modifications:
% 1. Added frametime read. ImageDataSeq starts with frametime in float64.
%    So, image data in each ImageDataSeq starts 8 bytes after DataStartPos.
%
% 2. fread entire image plane at once instead of fread line by line.
%    This is much faster.
%


%% get file info (width, height, nframes, nchannels)
% modified by TM from nd2finfo.m by Chao-Yuan Yeh

fid = fopen(file, 'r');

% Each datasegment begins with these signature bytes
sigbyte = [218, 206, 190, 10];

% Second segment always begins at 4096
fseek(fid, 4096, 'bof');
count = 1;
fs = struct;

flag = 0;
while flag == 0
  signature = fread(fid, 4, '*uint8')';
  if ~isempty(signature) && sum(signature == sigbyte) == 4
    fs(count).nameLength = fread(fid, 1, 'uint32');
    fs(count).dataLength = fread(fid, 1, 'uint64');
    fs(count).nameAttribute = fread(fid, fs(count).nameLength, '*char')';
    fs(count).dataStartPos = ftell(fid);
    flag = fseek(fid, fs(count).dataLength, 'cof');
    count = count + 1;
  else
    break
  end
end

for ii = 1 : length(fs)
  if strfind(fs(ii).nameAttribute, 'ImageAttributesLV!')
    img_attrib_idx = ii;
  elseif strfind(fs(ii).nameAttribute, 'ImageDataSeq|0')
    img_data_idx = ii;
  elseif strfind(fs(ii).nameAttribute, 'ImageMetadataSeq')
    img_metadata_idx = ii;
  end
end

fseek(fid, fs(img_attrib_idx).dataStartPos, 'bof');
img_attrib = fread(fid, fs(img_attrib_idx).dataLength, '*char')';

strloc(fid,fs, img_attrib_idx, img_attrib, 'uiWidth');
width = fread(fid, 1, '*uint32');
strloc(fid, fs, img_attrib_idx, img_attrib, 'uiHeight');
height = fread(fid, 1, '*uint32');
strloc(fid, fs, img_attrib_idx, img_attrib, 'uiSequenceCount');
nframe = fread(fid, 1, '*uint32');

fseek(fid, fs(img_metadata_idx).dataStartPos, 'bof');
img_metadata = fread(fid, fs(img_metadata_idx).dataLength, '*char')';

pos = strfind(img_metadata, insert0('ChannelIsActive'));
nchannel = length(pos);

%% read image data
% modified by TM from nd2read.m by Chao-Yuan Yeh

im_ch1 = zeros(width, height, nframe, 'uint16');
frametime = zeros(nframe,1);

% get start positions for all frames
all_dataStartPos = [fs(strncmp('ImageDataSeq',{fs(:).nameAttribute}, 12)).dataStartPos];

for myframe = 1:nframe
  fseek(fid,all_dataStartPos(myframe),'bof');
  
  % Each ImageDataSeq starts with frame time(ms) in float64.
  frametime(myframe) = fread(fid, 1, 'float64');
  
  % Image extracted from ND2 has image width defined by its first dimension.
  im_ch1(:, :,myframe) = reshape(fread(fid, nchannel * width * height, '*uint16'),...
    [height width]);
end

fclose(fid);

im_ch1 = permute(im_ch1, [2 1 3]);


% ND2 file has text sequence with 0 between characters.
  function out = insert0(in)
    num = in + 0;
    out = char([reshape(cat(1, num, zeros(size(num))), [1, length(num)*2]), 0]);
  end

  function strloc( fid, fs, fsidx, text, str )
    idx = strfind(text, insert0(str)) + length(insert0(str));
    fseek(fid, fs(fsidx).dataStartPos + idx(1), 'bof');
  end

end