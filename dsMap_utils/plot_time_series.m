function imH = plot_time_series(TSs,labels,ops)
% compactly plot time series

% INPUT
% TSs, a 2d array of time series or a cell of 2d arrays
% fps, frame per sec in sampling
% ops, specific options to control the appearance of figure

% OUTPUT
% imH, figure handel

% written by Ruix.Li in Jul,2020
% modified in Oct, 2020 add ops

fps = [];
defaultOps.normFlg = 0;
defaultOps.cm = lines(16);
defaultOps.spaceBetweenTraces = 3;
defaultOps.lineWidth = 1;
defaultOps.YTickL = [];
defaultOps.fontSize = 11;

if nargin < 4; ops = defaultOps; end
if ~isfield(ops,'norm'); ops.normFlg = defaultOps.normFlg; end
if ~isfield(ops,'cm'); ops.cm = defaultOps.cm; end
if ~isfield(ops,'spaceBetweenTraces'); ops.spaceBetweenTraces = defaultOps.spaceBetweenTraces; end
if ~isfield(ops,'lineWidth'); ops.lineWidth = defaultOps.lineWidth; end
if ~isfield(ops,'YTickL'); ops.YTickL = defaultOps.YTickL; end
if ~isfield(ops,'fontSize'); ops.fontSize = defaultOps.fontSize; end

if nargin < 3 || isempty(fps); fps = []; end
if ~iscell(TSs); TSs = {TSs}; end
nTS = length(TSs);
if nargin < 2 || isempty(labels)  
    labels = cell(1,nTS);
end

maxT = 0;
maxCh = 0;
maxD = [];

if nTS == 1; ops.cm = [0.1,0.1,0.1]; end

for i = 1:nTS
    if ops.normFlg
        base = median(TSs{i},1)+1;
        TSs{i} = bsxfun(@minus, TSs{i}, base);
        TSs{i} = bsxfun(@times, TSs{i}, 1./base);
    end
    maxT = max(maxT, size(TSs{i},1));
    maxCh = max(maxCh, size(TSs{i},2));
    for iCh = 1: size(TSs{i},2)
        maxD = cat(1, maxD, max(TSs{i}(:,iCh)) - min(TSs{i}(:,iCh)));
    end
    if isempty(labels{i})
        labels{i} = ['data #' num2str(i)];
    end
end

figWidth = round(min([1050 maxT*4]));
figHeight = round(max(200, min([600 maxCh*80])));
tcsFigSize = [100,100,figWidth,figHeight]; % [450,50,1050,600] is the max size

ncm = size(ops.cm,1);
if nargout > 0
h = figure('Position',tcsFigSize);
end

if isempty(ops.spaceBetweenTraces)
    intendStep = median(maxD,'omitnan');
else
    intendStep = ops.spaceBetweenTraces;
end
intendV = 0:intendStep:((maxCh-1)*intendStep);
hold on
for i = 1:nTS
    if isempty(fps)
        xVals = 1:size(TSs{i},1);
    else
        xVals = (0:(size(TSs{i},1)-1)) / double(fps);
    end
        
    for nCh = 1:size(TSs{i},2) 
     p(i) = plot(xVals, TSs{i}(:,nCh) + intendV(nCh),...
         'Linewidth',ops.lineWidth,'Color',ops.cm(rem(i-1,ncm)+1,:));
    end
    
end

if nTS > 1
    ld = legend(p,labels);
    ld.Box = 'off';
end

hold off
box off
YTicKs = intendV;
YTicKL = cell(1,maxCh);

if isempty(ops.YTickL) || length(ops.YTickL) ~= maxCh
for nCh = 1:maxCh
    YTicKL{nCh} = num2str(nCh);
end
else
    YTicKL = ops.YTickL;
end

if maxCh >= 30; YTicKL = []; end % dont show label when too many traces

set(gca,'XGrid','off','GridLineStyle','--',...
    'FontSize',ops.fontSize,'lineWidth',1.5,'FontWeight','Bold', ...
    'XAxisLocation','top','YTick',YTicKs,'YTickLabel',YTicKL,...
    'YLim',[-1*intendStep, intendV(end)+intendStep],'XLim',[0,max(xVals)]);

ylabel(['trial averaged response 1 u. = ' num2str(intendStep) 'SD'],'FontWeight','bold');

if isempty(fps)
    xlabel('Frame#','FontWeight','bold');
else
    xlabel('Time(s)','FontWeight','bold');
end 

if nargout > 0; imH = h; end


end
