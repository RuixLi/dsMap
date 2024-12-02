function info = read_mat_info(fname)
% read matfile data info
m = matfile(fname);
info.imHeight = size(m.stack,1);
info.imWidth = size(m.stack,2);
info.Loops = size(m.stack,3);
end