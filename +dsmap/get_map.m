function DSresult = get_map(ops)
% get HSV pseudocolor code direction and orientation selectivity map

% for every pixel in the data
% this function compute direction and orientation selectivity of every

% using vector-averaging method
% refer to M. Mazurek, et al. 2014 for more details

% INPUT
%
% cnfgDir directory of config
% - ops must have the following fields:
% ops.full_path, full path to the data
% ops.data_dir, directory of data
% ops.save_dir, directory to save results
% ops.file_type, ext of data file
% ops.dir_list, a list of directions of drifting gratings, usually 12 directions
% ops.trial_list, a ordered list of stimuli directions in experiment
% ops.trial_num, how many trials of experiment
% ops.ON_OFF_frames, [ON OFF] frame number per stimulus
% ops.smf, spatial smooth sigma
% ops.ON_idx, index of frames treated as stim ON
% ops.OFF_idx, index of frames treated as stim OFF or baseline

% written by Ruix.Li in Jun, 2021

% using cartesian coordinates, where 0 deg is north and increase ccw
% the grating deg in psychopy increases cw, so add minus sign
ops.t0 = tic;
dataDir = ops.full_path;
saveDir = ops.save_dir;
dirList = wrapTo360(-ops.dir_list);
trialList = wrapTo360(-ops.trial_list);
trialNum = ops.trial_num;
onOffFrames = ops.ON_OFF_frames;
trialLength = sum(onOffFrames);
expectT = trialNum * length(dirList) * trialLength;

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
        data = load(dataDir);
        stack = data.stack;
    case '.bin'
        stack = load_bin_data(dataDir);
    case '.raw'
        stack = load_raw_data(dataDir);
    case '.tif'
        stack = load_tif_data(dataDir);
    case '.tiff'
        stack = load_tif_data(dataDir);
    case '.nd2'
        stack = load_nd2_data(dataDir);
end

%% normalize and average stack over trial
stack = single(stack);
% stack = stack - min(stack(:)) + 1;
[d1,d2,Loops] = size(stack);
fprintf(' %d frames loaded\n',Loops)

assert(Loops == ops.Loops, 'data length not match, check parameters are correct');
nDir = length(dirList);

% calculate some statistic maps
disp('computing pixel statistics ...')
fovMap = median(stack,3);
maxMap = max_val_map(stack);
skewMap = skew_val_map(stack);
sdMap = sd_val_map(stack) + 1;
zMap = z_val_map(stack);
ccMap = double(corr_val_map(stack,8,1)); % per-normalize

% trail average data
disp('trial averaging data ...')
[avg,dirType] = sort_and_avg_rand_trial(stack,trialLength,trialList);
dirList = double(sort(dirList,'ascend'));
if ~isequal(dirType, dirList')
    warning('angles of directions may be wrong');
    warning('orientation selectivity may be wrong')
end

avgcat = reshape(avg,d1,d2,[]);
nFr = size(avgcat,3);

%% compute response in each direction
enhancer1 = imgaussfilt(sdMap,ops.smf);
enhancer2 = imgaussfilt(ccMap,ops.smf);

disp('compute directional responses for all pixels ...')

dirMap = zeros(d1,d2,nDir);
dirMapE = zeros(d1,d2,nDir);
dirMapE2 = zeros(d1,d2,nDir);

for iD = 1:length(dirList)
    dirA = avg(:,:,:,iD); % averaged amplitude of direction iD
    
    dirResp = calculate_response(dirA, ops);
    dirResp = imgaussfilt(dirResp,ops.smf);
    
    dirMap(:,:,iD) = dirResp;
    dirMapE(:,:,iD) = dirResp.*enhancer1;
    dirMapE2(:,:,iD) = dirResp.*enhancer2;
end

%% compute response in each orientation
% assume dirList has paired directions of the same orientation
% orientational response is the MIN response within the 2 directions
disp('compute orientational responses for all pixels ...')

nOri = nDir/2;

oriMap = zeros(d1,d2,nOri);
oriMapE1 = zeros(d1,d2,nOri);
oriMapE2 = zeros(d1,d2,nOri);

for iD = 1:nOri
    oriMap(:,:,iD) = min(cat(3,dirMap(:,:,iD),dirMap(:,:,iD+nOri)),[],3,'omitnan');
    oriMapE1(:,:,iD) = oriMap(:,:,iD).*enhancer1;
    oriMapE2(:,:,iD) = oriMap(:,:,iD).*enhancer2;
end

%% estimate direction and orientation tunning
disp('estimate directional tunning for all pixels ...')

oriList = dir_to_ori(dirList);
oriList = oriList(1:nOri);

dirList = dirList(:);
oriList = 2*oriList(:);

% change the degree to polar vectors
dirPolarVector = [cos(deg2rad(dirList)) sin(deg2rad(dirList))];
oriPolarVector = [cos(deg2rad(oriList)) sin(deg2rad(oriList))];

dirMat = reshape(dirMap,d1*d2,nDir);
dirRspVector = dirMat*dirPolarVector;

oriMat = reshape(oriMap,d1*d2,nOri);
oriRspVector = oriMat*oriPolarVector;

% polar vector angle
dirPrefAngle = wrapTo360(rad2deg(atan2(dirRspVector(:,2),dirRspVector(:,1))));
oriPrefAngle = wrapTo360(rad2deg(atan2(oriRspVector(:,2),oriRspVector(:,1))))/2;

% polar vector magnitude
dirVectorNorm = sqrt(dirRspVector(:,1).^2 + dirRspVector(:,2).^2);
oriVectorNorm = sqrt(oriRspVector(:,1).^2 + oriRspVector(:,2).^2);

% polar vector tune = sum(z) / sum(|z|)
dirTune = dirVectorNorm./sum(dirMat,2);
oriTune = oriVectorNorm./sum(oriMat,2);

%% HSV color code maps
dirAngleHue = hue_KO2HSV(dirPrefAngle/360);
oriAngleHue = hue_KO2HSV(oriPrefAngle/180);

if ops.rescale_tune
    dirTuneSat = rescale_map(dirTune,[0.4,1]);
    oriTuneSat = rescale_map(oriTune,[0.4,1]);
else
    dirTuneSat = dirTune;
    oriTuneSat = oriTune;
end

% HSV1: Hue = angle, Sat = dirPolarTune, Val = meanValue
meanVal = rescale_map(fovMap(:));
dirHSV1 = cat(2, dirAngleHue, dirTuneSat, meanVal);
oriHSV1 = cat(2, oriAngleHue, oriTuneSat, meanVal);
dirHSVMap1 = hsv2rgb(reshape(dirHSV1, d1, d2, 3));
oriHSVMap1 = hsv2rgb(reshape(oriHSV1, d1, d2, 3));

% HSV2: H = angle, S = dirPolarTune, V = maxValue
maxVal = rescale_map(maxMap(:));
dirHSV2 = cat(2, dirAngleHue, dirTuneSat, maxVal);
oriHSV2 = cat(2, oriAngleHue, oriTuneSat, maxVal);
dirHSVMap2 = hsv2rgb(reshape(dirHSV2, d1, d2, 3));
oriHSVMap2 = hsv2rgb(reshape(oriHSV2, d1, d2, 3));

% HSV3: H = angle, S = dirPolarTune, V = skewness
skewVal = rescale_map(skewMap(:));
dirHSV3 = cat(2, dirAngleHue, dirTuneSat, skewVal);
oriHSV3 = cat(2, oriAngleHue, oriTuneSat, skewVal);
dirHSVMap3 = hsv2rgb(reshape(dirHSV3, d1, d2, 3));
oriHSVMap3 = hsv2rgb(reshape(oriHSV3, d1, d2, 3));

% HSV4: H = angle, S = dirPolarTune, V = stdMap
sdVal = rescale_map(sdMap(:));
dirHSV4 = cat(2, dirAngleHue, dirTuneSat, sdVal);
oriHSV4 = cat(2, oriAngleHue, oriTuneSat, sdVal);
dirHSVMap4 = hsv2rgb(reshape(dirHSV4, d1, d2, 3));
oriHSVMap4 = hsv2rgb(reshape(oriHSV4, d1, d2, 3));

% HSV5: H = angle, S = dirPolarTune, V = ccMap
ccVal = rescale_map(ccMap(:));
dirHSV5 = cat(2, dirAngleHue, dirTuneSat, ccVal);
oriHSV5 = cat(2, oriAngleHue, oriTuneSat, ccVal);
dirHSVMap5 = hsv2rgb(reshape(dirHSV5, d1, d2, 3));
oriHSVMap5 = hsv2rgb(reshape(oriHSV5, d1, d2, 3));

%% HV color code, discard saturation values
s = ones(size(dirTune));

% HV1: H = angle, S = 1, V = meanValue
meanVal = rescale_map(fovMap(:));
dirHV1 = cat(2, dirAngleHue, s, meanVal);
oriHV1 = cat(2, oriAngleHue, s, meanVal);
dirHVMap1 = hsv2rgb(reshape(dirHV1, d1, d2, 3));
oriHVMap1 = hsv2rgb(reshape(oriHV1, d1, d2, 3));

% HV2: H = angle, S = 1, V = maxValue
maxVal = rescale_map(maxMap(:));
dirHV2 = cat(2, dirAngleHue, s, maxVal);
oriHV2 = cat(2, oriAngleHue, s, maxVal);
dirHVMap2 = hsv2rgb(reshape(dirHV2, d1, d2, 3));
oriHVMap2 = hsv2rgb(reshape(oriHV2, d1, d2, 3));

% HV3: H = angle, S = 1, V = skewness
skewVal = rescale_map(skewMap(:));
dirHV3 = cat(2, dirAngleHue, s, skewVal);
oriHV3 = cat(2, oriAngleHue, s, skewVal);
dirHVMap3 = hsv2rgb(reshape(dirHV3, d1, d2, 3));
oriHVMap3 = hsv2rgb(reshape(oriHV3, d1, d2, 3));

% HV4: H = angle, S = 1, V = stdMap
sdVal = rescale_map(sdMap(:));
dirHV4 = cat(2, dirAngleHue, s, sdVal);
oriHV4 = cat(2, oriAngleHue, s, sdVal);
dirHVMap4 = hsv2rgb(reshape(dirHV4, d1, d2, 3));
oriHVMap4 = hsv2rgb(reshape(oriHV4, d1, d2, 3));

% HV5: H = angle, S = 1, V = ccMap
ccVal = ccMap(:);
dirHV5 = cat(2, dirAngleHue, s, ccVal);
oriHV5 = cat(2, oriAngleHue, s, ccVal);
dirHVMap5 = hsv2rgb(reshape(dirHV5, d1, d2, 3));
oriHVMap5 = hsv2rgb(reshape(oriHV5, d1, d2, 3));

%% Angle Maps
dirAngleZHue = hue_KO2HSV(dirPrefAngle/360);
dirAngleZMap = cat(2, dirAngleZHue, s, s);
dirAngleZMap = hsv2rgb(reshape(dirAngleZMap,d1,d2,3));

oriAngleZHue = hue_KO2HSV(oriPrefAngle/180);
oriAngleZMap = cat(2, oriAngleZHue, s, s);
oriAngleZMap = hsv2rgb(reshape(oriAngleZMap,d1,d2,3));

%% Tune Maps
dirZTuneMap = reshape(dirTune,d1,d2);
oriZTuneMap = reshape(oriTune,d1,d2);
oriList = oriList/2;

%% plot figures
%mkdir(saveDir);
disp('saving files ...')
imwrite(rescale_map(fovMap),fullfile(saveDir,'01_FOV.tif'));

fH = figure('Position',[300,100,600,600]); set(fH, 'color', 'w', 'visible', 'off');
im_view(maxMap)
export_fig(fH, fullfile(fullfile(saveDir,'02_maxValueMap.tif')))

im_view(skewMap)
export_fig(fH, fullfile(fullfile(saveDir,'03_skewnessMap.tif')))

im_view(sdMap)
export_fig(fH, fullfile(fullfile(saveDir,'04_stdMap.tif')))

im_view(zMap)
export_fig(fH, fullfile(fullfile(saveDir,'05_zscoreMap.tif')))

im_view(ccMap)
export_fig(fH, fullfile(fullfile(saveDir,'06_ccMap.tif')))

im_view(dirZTuneMap)
export_fig(fH, fullfile(fullfile(saveDir,'37_dirTune_gDSI.tif')))

im_view(oriZTuneMap)
export_fig(fH, fullfile(fullfile(saveDir,'38_oriTune_gOSI.tif')))

set(gcf,'Position',[300,100,900,300]);
Ma = im_T(avgcat);
plot(mean(Ma,2),'LineWidth',2,'Color','k')
plot_stim(ops.ON_idx(1):trialLength:nFr,ops.ON_idx(end)-ops.ON_idx(1),[0.95,0.71,0.71])
plot_stim(ops.OFF_idx(1):trialLength:nFr,ops.OFF_idx(end)-ops.OFF_idx(1),[0.8,0.8,0.8])
xlim([0,size(Ma,1)])
ylabel('all pixel timecourse')
set(gca,'XTick',0:trialLength:nFr)
export_fig(fH, fullfile(fullfile(saveDir,'07_all_pixel_timecourse.tif')))

plot(mean(Ma(:,ccMap>0.2),2),'LineWidth',2,'Color','k')
plot_stim(ops.ON_idx(1):trialLength:nFr,ops.ON_idx(end)-ops.ON_idx(1),[0.95,0.71,0.71])
plot_stim(ops.OFF_idx(1):trialLength:nFr,ops.OFF_idx(end)-ops.OFF_idx(1),[0.8,0.8,0.8])
xlim([0,size(Ma,1)])
ylabel('pixel high Corr')
set(gca,'XTick',0:trialLength:nFr)
export_fig(fH, fullfile(fullfile(saveDir,'07_highCorr_pixel_timecourse.tif')))

plot(mean(Ma(:,sdMap>median(sdMap(:))),2),'LineWidth',2,'Color','k')
plot_stim(ops.ON_idx(1):trialLength:nFr,ops.ON_idx(end)-ops.ON_idx(1),[0.95,0.71,0.71])
plot_stim(ops.OFF_idx(1):trialLength:nFr,ops.OFF_idx(end)-ops.OFF_idx(1),[0.8,0.8,0.8])
xlim([0,size(Ma,1)])
ylabel('pixel high SD')
set(gca,'XTick',0:trialLength:nFr)
export_fig(fH, fullfile(fullfile(saveDir,'07_highSD_pixel_timecourse.tif')))

plot(mean(Ma(:,skewMap>1),2),'LineWidth',2,'Color','k')
plot_stim(ops.ON_idx(1):trialLength:nFr,ops.ON_idx(end)-ops.ON_idx(1),[0.95,0.71,0.71])
plot_stim(ops.OFF_idx(1):trialLength:nFr,ops.OFF_idx(end)-ops.OFF_idx(1),[0.8,0.8,0.8])
xlim([0,size(Ma,1)])
ylabel('pixel high skewness')
set(gca,'XTick',0:trialLength:nFr)
export_fig(fH, fullfile(fullfile(saveDir,'07_highSkewness_pixel_timecourse.tif')))

plot(mean(Ma(:,zMap>1),2),'LineWidth',2,'Color','k')
plot_stim(ops.ON_idx(1):trialLength:nFr,ops.ON_idx(end)-ops.ON_idx(1),[0.95,0.71,0.71])
plot_stim(ops.OFF_idx(1):trialLength:nFr,ops.OFF_idx(end)-ops.OFF_idx(1),[0.8,0.8,0.8])
xlim([0,size(Ma,1)])
ylabel('pixel high zscore')
set(gca,'XTick',0:trialLength:nFr)
export_fig(fH, fullfile(fullfile(saveDir,'07_highZdrift_pixel_timecourse.tif')))

set(gcf,'Position',[300,100,600,400]);
histogram(sdMap(:),'EdgeAlpha',0,'FaceColor','k')
title('pixel std')
export_fig(fH, fullfile(fullfile(saveDir,'11_hist_std.tif')))

histogram(skewMap(:),'EdgeAlpha',0,'FaceColor','k')
title('pixel skewness')
export_fig(fH, fullfile(fullfile(saveDir,'12_hist_skewness.tif')))

histogram(zMap(:),'EdgeAlpha',0,'FaceColor','k')
title('pixel Zscore')
export_fig(fH, fullfile(fullfile(saveDir,'14_hist_Zscore.tif')))

polarhistogram(deg2rad(dirPrefAngle(:)),8*nDir,'EdgeAlpha',0,'FaceColor','k')
title('pixel polar dir')
export_fig(fH, fullfile(fullfile(saveDir,'15_hist_prefDir.tif')))

histogram(dirTune(:),'EdgeAlpha',0,'FaceColor','k')
title('pixel tune dir')
export_fig(fH, fullfile(fullfile(saveDir,'16_hist_gDSI.tif')))

polarhistogram(deg2rad(oriPrefAngle(:)),4*nDir,'EdgeAlpha',0,'FaceColor','k')
title('pixel polar ori')
export_fig(fH, fullfile(fullfile(saveDir,'17_hist_prefOri.tif')))

histogram(oriTune(:),'EdgeAlpha',0,'FaceColor','k')
title('pixel tune ori')
export_fig(fH, fullfile(fullfile(saveDir,'18_hist_gOSI.tif')))

nrow = nDir/4;
set(fH, 'Position', [100,100,800,200*nrow+40])

r1 = quantile(dirMap(:),0.98);
r2 = quantile(dirMapE(:),0.98);
r3 = quantile(dirMapE2(:),0.98);

for iD = 1:nDir
    sub_tight_plot(nrow,4,iD)
    imagesc(dirMap(:,:,iD),[0,r1])
    axis tight; axis equal; axis off;
    colormap gray
    title(sprintf('%d deg',dirList(iD)))
end
export_fig(fH, fullfile(fullfile(saveDir,'21_resp_by_dir.tif')))

for iD = 1:nDir
    sub_tight_plot(nrow,4,iD)
    imagesc(dirMapE(:,:,iD),[0,r2])
    axis tight; axis equal; axis off;
    colormap gray
    title(sprintf('%d deg',dirList(iD)))
end
export_fig(fH, fullfile(fullfile(saveDir,'22_resp_by_dir_e1.tif')))

for iD = 1:nDir
    sub_tight_plot(nrow,4,iD)
    imagesc(dirMapE2(:,:,iD),[0,r3])
    axis tight; axis equal; axis off;
    colormap gray
    title(sprintf('%d deg',dirList(iD)))
end
export_fig(fH, fullfile(fullfile(saveDir,'22_resp_by_dir_e2.tif')))

clf(fH)

for iD = 1:nOri
    sub_tight_plot(nrow,4,iD)
    imagesc(oriMap(:,:,iD),[0,r1])
    axis tight; axis equal; axis off;
    colormap gray
    title(sprintf('%d deg',oriList(iD)))
end
export_fig(fH, fullfile(fullfile(saveDir,'23_resp_by_ori.tif')))

for iD = 1:nOri
    sub_tight_plot(nrow,4,iD)
    imagesc(oriMapE1(:,:,iD),[0,r2])
    axis tight; axis equal; axis off;
    colormap gray
    title(sprintf('%d deg',oriList(iD)))
end

export_fig(fH, fullfile(fullfile(saveDir,'24_resp_by_ori_e1.tif')))

for iD = 1:nOri
    sub_tight_plot(nrow,4,iD)
    imagesc(oriMapE2(:,:,iD),[0,r3])
    axis tight; axis equal; axis off;
    colormap gray
    title(sprintf('%d deg',oriList(iD)))
end

export_fig(fH, fullfile(fullfile(saveDir,'24_resp_by_ori_e2.tif')))

close(fH)

% write pesudocolor maps

imwrite(dirHSVMap1,fullfile(saveDir,'31_dirHSV1_meanValue.tif'));
imwrite(dirHSVMap2,fullfile(saveDir,'31_dirHSV2_maxValue.tif'));
imwrite(dirHSVMap3,fullfile(saveDir,'31_dirHSV3_skewValue.tif'));
imwrite(dirHSVMap4,fullfile(saveDir,'31_dirHSV4_stdValue.tif'));
imwrite(dirHSVMap5,fullfile(saveDir,'31_dirHSV5_ccValue.tif'));

imwrite(oriHSVMap1,fullfile(saveDir,'32_oriHSV1_meanValue.tif'));
imwrite(oriHSVMap2,fullfile(saveDir,'32_oriHSV2_maxValue.tif'));
imwrite(oriHSVMap3,fullfile(saveDir,'32_oriHSV3_skewValue.tif'));
imwrite(oriHSVMap4,fullfile(saveDir,'32_oriHSV4_stdValue.tif'));
imwrite(oriHSVMap5,fullfile(saveDir,'32_oriHSV5_ccValue.tif'));

imwrite(dirHVMap1,fullfile(saveDir,'33_dirHV1_meanValue.tif'));
imwrite(dirHVMap2,fullfile(saveDir,'33_dirHV2_maxValue.tif'));
imwrite(dirHVMap3,fullfile(saveDir,'33_dirHV3_skewValue.tif'));
imwrite(dirHVMap4,fullfile(saveDir,'33_dirHV4_stdValue.tif'));
imwrite(dirHVMap5,fullfile(saveDir,'33_dirHV5_ccValue.tif'));

imwrite(oriHVMap1,fullfile(saveDir,'34_oriHV1_meanValue.tif'));
imwrite(oriHVMap2,fullfile(saveDir,'34_oriHV2_maxValue.tif'));
imwrite(oriHVMap3,fullfile(saveDir,'34_oriHV3_skewValue.tif'));
imwrite(oriHVMap4,fullfile(saveDir,'34_oriHV4_stdValue.tif'));
imwrite(oriHVMap5,fullfile(saveDir,'34_oriHV5_ccValue.tif'));

imwrite(dirAngleZMap,fullfile(saveDir,'35_prefDirMap.tif'));
imwrite(oriAngleZMap,fullfile(saveDir,'36_prefOriMap.tif'));

generate_color_wheel(nDir,saveDir)

%% save data
DSresult.ops = ops;

DSresult.trialAvgData = avg;
DSresult.FOV = fovMap;
DSresult.maxValueMap = maxMap;
DSresult.stdMap = sdMap;
DSresult.skewMap = skewMap;
DSresult.ccMap = ccMap;

DSresult.dirRespMap = dirMap;
DSresult.oriRespMap = oriMap;

DSresult.dirRspVector = dirRspVector;
DSresult.oriRspVector = oriRspVector;

DSresult.dirPrefAngle = dirPrefAngle;
DSresult.dirTune = dirTune;

DSresult.oriPrefAngle = oriPrefAngle;
DSresult.oriTune = oriTune;

fprintf('finished in %.2f s\n',toc(ops.t0))
end
