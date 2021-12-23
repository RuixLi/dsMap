function ccMap = corr_val_map(Y,sz,normFlg)

% calculate correlation map for cell detection based on neighboing pixels

% INPUT
% Y[XYT], raw 2p imaging data
% sz[int/vector/matrix], define the relative location of neighbours. 
% -int, number of nearest neighbours, either 4(default)
% -vector, [rmin, rmax], the range of neighbour distance
% -matrix: a squared matrix (2k+1)*(2k+1)indicating the location of neighbours
% d1,d2: spatial dimensions (required with Y is matrix)
% normFlg: indicate whether Y need normalization (default 1)

% OUTPUT
% ccMap[matrix], correlatin map
% -if no specified output, show ccMap directly

% adapt from Eftychios A. Pnevmatikakis, Simons Foundation, 2015
% -modified by Ruix. Li, in Jun, 2020

% preprocess the raw data
if nargin<3; normFlg = 1; end
if ~strcmpi(class(Y),'single'); Y = single(Y); end
if (nargin<2) || isempty(sz); sz = [0,1,0; 1,0,1; 0,1,0]; end

[d1,d2,~] = size(Y);    % image dimension

if normFlg
    % centering
    mY = mean(Y,3);
    Y = bsxfun(@minus, Y, mY);
    sY = sqrt(mean(Y.*Y,3));
    Y = bsxfun(@times, Y, 1./(sY+1));
end

% construct a matrix indicating location of the matrix
if  isscalar(sz)
    if sz == 8      % 8 nearest neighbours
        sz = [1,1,1; 1,0,1; 1,1,1];
    elseif sz==4
        sz = [0,1,0; 1,0,1; 0,1,0];
    end
elseif length(sz(:)) == 2
    % the specified neighbours has a distance within the domain [dmin,
    % dmax)
    sz = ceil(sz);
    dmin = min(sz); dmax = max(sz);
    rsub = (-dmax+1):(dmax-1);      % row subscript
    csub = rsub;      % column subscript
    [cind, rind] = meshgrid(csub, rsub);
    R = sqrt(cind.^2+rind.^2);
    sz = (R>=dmin) .* (R<dmax);
end

% compute the correlation
Yconv = imfilter(Y, sz,'conv');        % sum over the neighbouring pixels
MASK = imfilter(ones(d1,d2), sz, 'conv');   % count the number of neighbouring pixels
Cn = mean(Yconv.*Y, 3)./MASK;   % compute correlation and normalize

if nargout < 1
    figure; im_view(Cn);
else
    ccMap = Cn;
end

end