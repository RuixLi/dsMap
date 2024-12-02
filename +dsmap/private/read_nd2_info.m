function info = read_nd2_info(file)

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
  end
end

fseek(fid, fs(img_attrib_idx).dataStartPos, 'bof');
img_attrib = fread(fid, fs(img_attrib_idx).dataLength, '*char')';

strloc(fid,fs, img_attrib_idx, img_attrib, 'uiWidth');
info.imWidth = double(fread(fid, 1, '*uint32'));
strloc(fid, fs, img_attrib_idx, img_attrib, 'uiHeight');
info.imHeight = double(fread(fid, 1, '*uint32'));
strloc(fid, fs, img_attrib_idx, img_attrib, 'uiSequenceCount');
info.Loops = double(fread(fid, 1, '*uint32'));

end

function out = insert0(in)
num = in + 0;
out = char([reshape(cat(1, num, zeros(size(num))), [1, length(num)*2]), 0]);
end

function strloc( fid, fs, fsidx, text, str )
idx = strfind(text, insert0(str)) + length(insert0(str));
fseek(fid, fs(fsidx).dataStartPos + idx(1), 'bof');
end