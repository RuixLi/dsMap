function DSstat = get_stat(FP,ops)
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

ops.t0 = tic;
dataDir = ops.full_path;
saveDir = ops.save_dir;
dirList = wrapTo360(-ops.dir_list);
trialList = wrapTo360(-ops.trial_list);
trialNum = ops.trial_num;
onOffFrames = ops.ON_OFF_frames;
trialLength = sum(onOffFrames);
expectT = trialNum * length(dirList) * trialLength;

if ischar(FP)
    f = load(FP,'FP');
    FP = f.FP;
end

% read data size if not provided, this could be slow
if any([isempty(ops.d1), isempty(ops.d2), isempty(ops.Loops)]) && ~isempty(ops.file_type)
    fprintf('getting data info ...')
    switch ops.file_type
        case 'mat'
            info = read_mat_info(dataDir);
        case 'raw'
            info = read_bin_info(dataDir, [ops.d1, ops.d2]);
        case 'bin'
            info = read_bin_info(dataDir, [ops.d1, ops.d2]);
        case 'tif'
            info = read_tif_info(dataDir);
        case 'tiff'
            info = read_tif_info(dataDir);
        case 'nd2'
            info = read_nd2_info(dataDir);
    end
    
    if expectT ~= info.Loops; error('data length not match, check parameters are correct'); end
    fprintf('done\n')
end

fprintf('loading data ...')
switch ops.file_type
    case '.mat'
        d = load(dataDir);
        data = d.stack;
    case '.bin'
        data = load_bin_data(dataDir);
    case '.raw'
        data = load_raw_data(dataDir);
    case '.tif'
        data = load_tif_data(dataDir);
    case '.tiff'
        data = load_tif_data(dataDir);
    case '.nd2'
        data = load_nd2_data(dataDir);
end

[d1,d2,Loops] = size(data);
fprintf(' %d frames loaded\n',Loops)
assert(Loops == ops.Loops, 'data length not match, check parameters are correct');
assert(d1*d2 == size(FP,1), 'data size not match, check FP are correct');

data = single(data);
FOV = mean(data,3);
FOV = uint8(rescale(FOV,0,255));
TS = im_T(data)*full(FP);

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
    dirRspAvg(:,i) = max(DFavg(ops.ON_idx,:,i),[],1) - mean(DFavg(ops.OFF_idx,:,i),1);
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
dirRsp = squeeze(max(trialTS(:,ops.ON_idx,:,:),[],2)-mean(trialTS(:,ops.OFF_idx,:,:),2));
oriRsp = zeros(nOri,trialNum,K);
for i = 1:nOri
    oriRsp(i,:,:) = min(dirRsp(i,:,:),dirRsp(i+nOri,:,:));
end

%% fit tuning curve with wraped gaussian
disp('fit tuning curve with wraped gaussian ...')
[dirTunCurv,dirStat] = gauss_fit_dir_tuning(dirRspAvg,dirList);
[oriTunCurv,oriStat] = gauss_fit_ori_tuning(oriRspAvg,oriList);

%% 1) FP
disp('plot results ...')

im = label_neuron_FP(FP,uint8(FOV));
imwrite(im,fullfile(saveDir,'1_neruon_FP.tif'))
im = label_neuron_contour(FP,uint8(FOV));
imwrite(im,fullfile(saveDir,'1_neruon_contour.tif'))
%% 2) averaged timecourses of all neurons
imh = figure('Position',[200,10,500,900],'Visible','off');
stimOnset = 0:trialLength:nDir*trialLength;
n = 1:24:K;
n = [n,K+1];
p = 1;
while n(p) <= K
    clf(imh)
    %figure(imh)
    plot_time_series(DFtrialAvg(:,n(p):n(p+1)-1),[],n(p):n(p+1)-1)
    plot_stim(ops.ON_idx(1):trialLength:Loops,ops.ON_idx(end)-ops.ON_idx(1),[0.95,0.71,0.71])
    plot_stim(ops.OFF_idx(1):trialLength:Loops,ops.OFF_idx(end)-ops.OFF_idx(1),[0.8,0.8,0.8])
    set(gca,'XTick',stimOnset,'XTickLabel',dirList)
    xlabel('Deg')
    export_fig(imh, fullfile(fullfile(saveDir,sprintf('2_trial_averaged_time_series_%d.tif',p))))
    p = p + 1;
end
close(imh)

%% 3) averaged response in each direction
theta = 2*pi*dirList/360;
n = 1:24:K;
n = [n,K+1];
p = 1;
s = 1;
imH = figure('Position',[100,100,800,550],'Visible','off');
while n(p) <= K
    clf(imH)
    %figure(imH)
    for i = n(p):n(p+1)-1
        sub_tight_plot(4,6,s,[0.04,0.04]);
        polarplot([theta;0],[dirRspAvg(i,:) dirRspAvg(i,1)]','lineWidth',1.5)
        set(gca,'fontSize',6)
        s = s+1;
        title(sprintf('neuron #%d, gDSI %.2f',i,dirTune(i)))
    end
    export_fig(imH, fullfile(fullfile(saveDir,sprintf('3_trial_averaged_dir_response_%d.tif',p))))
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
imH = figure('Position',[100,100,800,550],'Visible','off');
while n(p) <= K
    clf(imH)
    %figure(imH)
    for i = n(p):n(p+1)-1
        sub_tight_plot(4,6,s,[0.04,0.04]);
        polarplot(theta,oriRspAvg(i,:)','lineWidth',1.5)
        set(gca,'fontSize',6)
        s = s+1;
        title(sprintf('neuron #%d, gOSI %.2f',i,oriTune(i)))
    end
    export_fig(imH, fullfile(fullfile(saveDir,sprintf('3_trial_averaged_ori_response_%d.tif',p))))
    p = p+1;
    s = 1;
end

close(imH)

%% 4) plot prefered ori/dir histograms by vector average
imh = figure('Position',[100,100,600,600],'Visible','off');
%figure(imh)
polarhistogram(deg2rad(dirPrefAngle),nDir,'EdgeAlpha',0,'FaceColor','k')
export_fig(imh, fullfile(fullfile(saveDir,'4_pref_Dir_Ang_hist.tif')))
clf(imh)

%figure(imh)
polarhistogram(deg2rad(oriPrefAngle),nOri,'EdgeAlpha',0,'FaceColor','k')
export_fig(imh, fullfile(fullfile(saveDir,'4_pref_Ori_Ang_hist.tif')))
clf(imh)

%% plot prefered ori/dir tuning (gOSI/gDSI) by vector average
%figure(imh)
histogram(dirTune,0:0.1:1,'EdgeAlpha',0,'FaceColor','k')
export_fig(imh, fullfile(fullfile(saveDir,'4_gDSI_hist.tif')))
clf(imh)

%figure(imh)
histogram(oriTune,0:0.1:1,'EdgeAlpha',0,'FaceColor','k')
export_fig(imh, fullfile(fullfile(saveDir,'4_gOSI_hist.tif')))
clf(imh)

close(imh)

%% 5) plot dir tuning curves
n = 1:24:K;
n = [n,K+1];
p = 1;
s = 1;
imH = figure('Position',[100,100,800,550],'Visible','off');
while n(p) <= K
    clf(imH)
    %figure(imH)
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
    export_fig(imH, fullfile(fullfile(saveDir,sprintf('5_dir_tuning_curve_%d.tif',p))))
    p = p+1;
    s = 1;
end

close(imH)

%% plot ori tuning curves
n = 1:24:K;
n = [n,K+1];
p = 1;
s = 1;
imH = figure('Position',[100,100,800,550],'Visible','off');
while n(p) <= K
    clf(imH)
    %figure(imH)
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
    export_fig(imH, fullfile(fullfile(saveDir,sprintf('5_ori_tuning_curve_%d.tif',p))))
    p = p+1;
    s = 1;
end

close(imH)

%% 6) plot prefered ori/dir histograms by vector average
imh = figure('Position',[100,100,600,600],'Visible','off');
%figure(imh)
polarhistogram(deg2rad([dirStat.prefDir]),nDir,'EdgeAlpha',0,'FaceColor','r')
export_fig(gcf, fullfile(fullfile(saveDir,'6_fitPrefDirAng_hist.tif')))
clf(imh)

polarhistogram(deg2rad([oriStat.prefOri]),nOri,'EdgeAlpha',0,'FaceColor','r')
export_fig(gcf, fullfile(fullfile(saveDir,'6_fitPrefOriAng_hist.tif')))
clf(imh)

%% plot prefered ori/dir tuning by vector average
histogram([dirStat.dsi],0:0.1:1,'EdgeAlpha',0,'FaceColor','r')
export_fig(gcf, fullfile(fullfile(saveDir,'6_fitDSI_hist.tif')))
clf(imh)

histogram([oriStat.osi],0:0.1:1,'EdgeAlpha',0,'FaceColor','r')
export_fig(gcf, fullfile(fullfile(saveDir,'6_fitOSI_hist.tif')))
clf(imh)

close(imh)

%% save data
DSstat.ops = ops;
DSstat.FOV = FOV;
DSstat.FP = FP;
DSstat.TS = TS;
DSstat.trialAvgDF = DFtrialAvg;
DSstat.dirRsp = dirRsp;
DSstat.dirRspTrialAvg = dirRspAvg;
DSstat.oriRsp = oriRsp;
DSstat.oriRspTrialAvg = oriRspAvg;
DSstat.dirRspVector = dirRspVector;
DSstat.oriRspVector = oriRspVector;
DSstat.dirPrefAngle = dirPrefAngle;
DSstat.oriPrefAngle = oriPrefAngle;
DSstat.dirTune = dirTune;
DSstat.oriTune = oriTune;
DSstat.dirTuningCurve = dirTunCurv;
DSstat.oriTuningCurve = oriTunCurv;
DSstat.dirStat = dirStat;
DSstat.oriStat = oriStat;


fprintf('finished in %.2f s\n',toc(ops.t0))
close all

end
