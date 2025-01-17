ops.full_path = 'N:\imaging\Hieda\241127\241018AF03-FOV2\DFT_006.nd2'; % full path to the data
[dataDir, fileName, ext] = fileparts(ops.full_path);
if isempty(dataDir);dataDir = pwd; end
ops.data_dir = dataDir; % full path to the folder
ops.save_dir = fullfile(dataDir,[fileName '-DSstat']); % the dir of save folder
ops.file_name = fileName;
ops.file_type = ext;
ops.d1 = 512; % can be empty if not known
ops.d2 = 512;
ops.dir_list = [0 45 90 135 180 225 270 315]; % direction list
ops.trial_list = []; % in case if the trials are in random order
ops.trial_num = 10; % number of trials
ops.ON_OFF_frames = [75,75];
ops.ON_idx = 76:150; % index of ON frames as the response
ops.OFF_idx = 50:75; % index of OFF frames as the baseline activity
ops.Loops = sum(ops.ON_OFF_frames)*ops.trial_num*numel(ops.dir_list); % number of frames in the data


if isempty(ops.trial_list)
    ops.trial_list = repmat(ops.dir_list,1,ops.trial_num);
end

if exist(ops.save_dir, 'dir') == 0
    mkdir(ops.save_dir);
else
    warning('The folder already exists, the data will be overwritten');
end

DSstat = dsmap.get_stat(FP,ops);

save(fullfile(ops.save_dir,[ops.file_name '-DSstat.mat']),'DSstat');