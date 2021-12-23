function getOriStat(dataDir,fpDir,cnfgDir)
% for given neuron FP
% calculate statistics direction and orientation selectivity

% INPUT
% datadir[string] directory of data
% - data type supported:
% MATLAB .mat file, binary .bin file and THORLABS .raw file

% fpDir[string] directory of matfile where neruon footprint is stored

% config directory of config
% - config must have the following fields:
% config.animalID, animal ID
% config.TTStamp, a string of time stamp
% config.dirList, a list of stimuli direction, usually 12 directions
% config.trialList, a ordered list of stimuli directions in experiment
% config.trialNum, how many trials of experiment
% config.onOffFrames, [ON OFF] frame number per stimulus

onIndx = 2:5; % index of frames treated as stim ON
offIndx = 6:8; % index of frames treated as stim OFF or baseline

if nargin == 0
    [dataName,dataPath] = uigetfile('*.*','select data file ...');
    if dataName == 0; return; end
    dataDir = fullfile(dataPath, dataName);
end

if nargin < 2 || isempty(fpDir)
[fpName,fpPath] = uigetfile('*.mat','select FP file ...');
if fpName == 0; return; end
fpDir = fullfile(fpPath, fpName);
end

if nargin < 3 || isempty(cnfgDir)
[configName,configPath] = uigetfile('*.mat','select config file ...');
if configName == 0; return; end
configDir = fullfile(configPath, configName);
end

savDir = fileparts(configDir);

%% load data, FP and config
d = load(fpDir,'FP');
FP = d.FP;
cnfg = load(configDir);

DTstamp = cnfg.TTStamp;
session = cnfg.sessionName;
% to polar coordinates, where 0 deg is upward and increase clockwise
% the grating deg in psychopy increases cw, so add minus sign
dirList = wrapTo360(-cnfg.dirList);
trialList = wrapTo360(-cnfg.trialList);
trialNum = cnfg.trialNum;
onOffFrames = cnfg.onOffFrames;
%onOffFrames = config.offOnFrames;
trialLength = sum(onOffFrames);

[~, ~, ext] = fileparts(dataDir);

switch ext
    case '.mat'
        info = read_mat_info(dataDir);
    case '.raw'
        info = read_bin_info(dataDir);
    case '.bin'
        info = read_bin_info(dataDir);
end

expectT = trialNum * length(dirList) * trialLength;

if expectT ~= info.Loops; error('data and configration not match'); end

switch ext
    case '.mat'
        d = load(dataDir);
        data = d.stack;
    case '.bin'
        data = load_bin_data(dataDir);
    case '.raw'
        data = load_raw_data(dataDir);
end

data = single(data);
FOV = mean(data,3);
FOV = uint8(rescale(FOV,0,255));
TS = im_T(data)*FP;

%% compute averaged neuron time series
K = size(TS,2);
nDir = length(dirList);
base = std(TS,0,1)+1;
DF = (TS-repmat(median(TS,1),1))./repmat(base,size(TS,1),1);
[DFavg,dirList,~,trialInd] = sort_and_avg_rand_trial(DF,trialLength,trialList);
DFtrialAvg = reshape(permute(DFavg,[1,3,2]),nDir*trialLength,K);

%% compute averaged neuron response in each direction and orientation
nOri = nDir/2;
dirRspAvg = zeros(K,nDir);
oriRspAvg = zeros(K,nOri);

for i = 1:nDir
    dirRspAvg(:,i) = max(DFavg(onIndx,:,i),[],1) - mean(DFavg(offIndx,:,i),1);
end

dirRspAvg(dirRspAvg<0) = 0.00001;

oriList = dir_to_ori(dirList);
oriList = oriList(1:nOri);

for i = 1:nOri
    oriRspAvg(:,i) = min(dirRspAvg(:,i),dirRspAvg(:,i+nOri));
end

%% compute prefered direction/orientation and tune
dirList = dirList(:);
oriList = 2*oriList(:);

dirPolarVector = [cos(deg2rad(dirList)) sin(deg2rad(dirList))];
oriPolarVector = [cos(deg2rad(oriList)) sin(deg2rad(oriList))];

dirRspVector = dirRspAvg*dirPolarVector;
oriRspVector = oriRspAvg*oriPolarVector;

dirPrefAngle = wrapTo360(rad2deg(atan2(dirRspVector(:,2),dirRspVector(:,1))));
oriPrefAngle = wrapTo360(rad2deg(atan2(oriRspVector(:,2),oriRspVector(:,1))))/2;

dirVectorNorm = sqrt(dirRspVector(:,1).^2+dirRspVector(:,2).^2);
oriVectorNorm = sqrt(oriRspVector(:,1).^2+oriRspVector(:,2).^2);

dirTune = dirVectorNorm./sum(dirRspAvg,2);
oriTune = oriVectorNorm./sum(oriRspAvg,2);

oriList = oriList./2;

%% fit tuning curves
trialTS = zeros(nDir,trialLength,trialNum,K);
TSt = reshape(DF,trialLength,nDir*trialNum,K);
for i = 1:nDir
trialTS(i,:,:,:) = TSt(:,trialInd(:,i),:);
end
dirRsp = squeeze(max(trialTS(:,onIndx,:,:),[],2)-mean(trialTS(:,offIndx,:,:),2));
oriRsp = zeros(nOri,trialNum,K);
for i = 1:nOri
oriRsp(i,:,:) = min(dirRsp(i,:,:),dirRsp(i+nOri,:,:));
end

%% fit tuning curve with wraped gaussian
disp('fit tuning curve with wraped gaussian ...')
[dirTunCurv,dirStat] = gauss_fit_dir_tuning(dirRspAvg,dirList);
[oriTunCurv,oriStat] = gauss_fit_ori_tuning(oriRspAvg,oriList);

%%
savDir = fullfile(savDir,['dsStat-' session '-' DTstamp]);
mkdir(savDir);

%% 1) FP
disp('plot results ...')

im = label_neuron_FP(FP,uint8(FOV));
imwrite(im,fullfile(savDir,'1_neruon_FP.tif'))

%% 2) averaged timecourses of all neurons
imh = figure('Position',[1700,10,700,1000]);
stimOnset = 0:trialLength:nDir*trialLength;
plot_time_series(DFtrialAvg)
plot_stim(0:trialLength:nDir*trialLength,onOffFrames(1))
set(gca,'XTick',stimOnset,'XTickLabel',dirList)
xlabel('Deg')
export_fig(imh, fullfile(fullfile(savDir,'2_trial_averaged_time_series.tif')))
close(imh)

%% 3) averaged response in each direction
theta = 2*pi*dirList/360;
n = 1:24:K;
n = [n,K+1];
p = 1;
s = 1;
imH = figure('Position',[100,100,800,550]);
while n(p) <= K
clf(imH)
figure(imH)
for i = n(p):n(p+1)-1
    sub_tight_plot(4,6,s,[0.04,0.04]);
    polarplot([theta;0],[dirRspAvg(i,:) dirRspAvg(i,1)]','lineWidth',1.5)
    set(gca,'fontSize',6)
    s = s+1;
    title(sprintf('neuron #%d, gDSI %.2f',i,dirTune(i)))
end
export_fig(imH, fullfile(fullfile(savDir,sprintf('3_trial_averaged_dir_response_%d.tif',p))))
p = p+1;
s = 1;
end

close(imH)

%% 3) averaged response in each orientation
theta = pi*oriList/180;
n = 1:24:K;
n = [n,K+1];
p = 1;
s = 1;
imH = figure('Position',[100,100,800,550]);
while n(p) <= K
clf(imH)
figure(imH)
for i = n(p):n(p+1)-1
    sub_tight_plot(4,6,s,[0.04,0.04]);
    polarplot(theta,oriRspAvg(i,:)','lineWidth',1.5)
    set(gca,'fontSize',6)
    s = s+1;
    title(sprintf('neuron #%d, gOSI %.2f',i,oriTune(i)))
end
export_fig(imH, fullfile(fullfile(savDir,sprintf('3_trial_averaged_ori_response_%d.tif',p))))
p = p+1;
s = 1;
end

close(imH)

%% 4) plot prefered ori/dir histograms by vector average
imh = figure('Position',[-700,100,600,600]);
figure(imh)
polarhistogram(deg2rad(dirPrefAngle),nDir,'EdgeAlpha',0,'FaceColor','k')
export_fig(imh, fullfile(fullfile(savDir,'4_pref_Dir_Ang_hist.tif')))
clf(imh)

figure(imh)
polarhistogram(deg2rad(oriPrefAngle),nOri,'EdgeAlpha',0,'FaceColor','k')
export_fig(imh, fullfile(fullfile(savDir,'4_pref_Ori_Ang_hist.tif')))
clf(imh)

%% plot prefered ori/dir tuning (gOSI/gDSI) by vector average
figure(imh)
histogram(dirTune,0:0.1:1,'EdgeAlpha',0,'FaceColor','k')
export_fig(imh, fullfile(fullfile(savDir,'4_gDSI_hist.tif')))
clf(imh)

figure(imh)
histogram(oriTune,0:0.1:1,'EdgeAlpha',0,'FaceColor','k')
export_fig(imh, fullfile(fullfile(savDir,'4_gOSI_hist.tif')))
clf(imh)

close(imh)

%% 5) plot dir tuning curves
n = 1:24:K;
n = [n,K+1];
p = 1;
s = 1;
imH = figure('Position',[100,100,800,550]);
while n(p) <= K
clf(imH)
figure(imH)
for i = n(p):n(p+1)-1
    subplot(4,6,s)
    plot(dirList,dirRsp(:,:,i),'lineWidth',0.5,'Color',[0.7,0.7,0.7])
    hold on
    plot(dirList,dirRspAvg(i,:),'lineWidth',1.5,'Color','k')
    plot(0:1:359,dirTunCurv(i,:),'lineWidth',1,'Color','r')
    xlim([0,359])
    set(gca,'fontSize',5,'XTick',2*oriList)
    s = s+1;
    title(sprintf('neuron #%d, DSI %.2f',i,dirStat(i).dsi))
end
export_fig(imH, fullfile(fullfile(savDir,sprintf('5_dir_tuning_curve_%d.tif',p))))
p = p+1;
s = 1;
end

close(imH)

%% plot ori tuning curves
n = 1:24:K;
n = [n,K+1];
p = 1;
s = 1;
imH = figure('Position',[100,100,800,550]);
while n(p) <= K
clf(imH)
figure(imH)
for i = n(p):n(p+1)-1
    subplot(4,6,s)
    plot(oriList,oriRsp(:,:,i),'lineWidth',0.5,'Color',[0.7,0.7,0.7])
    hold on
    plot(oriList,oriRspAvg(i,:),'lineWidth',1.5,'Color','k')
    plot(0:1:179,oriTunCurv(i,:),'lineWidth',1,'Color','r')
    xlim([0,179])
    set(gca,'fontSize',5,'XTick',oriList)
    s = s+1;
    title(sprintf('neuron #%d, OSI %.2f',i,oriStat(i).osi))
end
export_fig(imH, fullfile(fullfile(savDir,sprintf('5_ori_tuning_curve_%d.tif',p))))
p = p+1;
s = 1;
end

close(imH)

%% 6) plot prefered ori/dir histograms by vector average
imh = figure('Position',[-700,100,600,600]);
figure(imh)
polarhistogram(deg2rad([dirStat.prefDir]),nDir,'EdgeAlpha',0,'FaceColor','r')
export_fig(gcf, fullfile(fullfile(savDir,'6_fitPrefDirAng_hist.tif')))
clf(imh)

polarhistogram(deg2rad([oriStat.prefOri]),nOri,'EdgeAlpha',0,'FaceColor','r')
export_fig(gcf, fullfile(fullfile(savDir,'6_fitPrefOriAng_hist.tif')))
clf(imh)

%% plot prefered ori/dir tuning by vector average
histogram([dirStat.dsi],0:0.1:1,'EdgeAlpha',0,'FaceColor','r')
export_fig(gcf, fullfile(fullfile(savDir,'6_fitDSI_hist.tif')))
clf(imh)

histogram([oriStat.osi],0:0.1:1,'EdgeAlpha',0,'FaceColor','r')
export_fig(gcf, fullfile(fullfile(savDir,'6_fitOSI_hist.tif')))
clf(imh)

close(imh)

%% save data
ds.FOV = FOV;
ds.FP = FP;
ds.TS = trialTS;
ds.TStrialAvg = DFtrialAvg;
ds.dirRsp = dirRsp;
ds.dirRspTrialAvg = dirRspAvg;
ds.oriRsp = oriRsp;
ds.oriRspTrialAvg = oriRspAvg;
ds.dirRspVector = dirRspVector;
ds.oriRspVector = oriRspVector;
ds.dirPrefAngle = dirPrefAngle;
ds.oriPrefAngle = oriPrefAngle;
ds.dirTune = dirTune;
ds.oriTune = oriTune;
ds.dirTuningCurve = dirTunCurv;
ds.oriTuningCurve = oriTunCurv;
ds.dirStat = dirStat;
ds.oriStat = oriStat;

save(fullfile(savDir,['dsStat-' session '-' DTstamp '.mat']), '-struct', 'ds')

disp('data saved successfully ...')


end
