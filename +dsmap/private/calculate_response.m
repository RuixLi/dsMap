function R = calculate_response(A, ops)

if ~isfield(ops, 'movmean_win')
    ops.movmean_win = 0;
end

if ~isfield(ops, 'norm_method')
    ops.norm_method = 'dF/Fb';
end

if ~isfield(ops, 'resp_type')
    ops.resp_type = 'AUC';
end

% data denoising with moving average filter
if ops.movmean_win > 0
    A = movmean(A, ops.movmean_win, 3);
end

% data normalization

% method3: z-score normalization
switch ops.norm_method
    case 'raw' % method1: raw data
        R = A;
        
    case 'dF/Fa' % method2: dF/F normalization, F0 is the mean of the whole movie
        R = (A - mean(A, 3)) ./ mean(A, 3);
        
    case 'dF/Fb' % method3: dF/F0 normalization, F0 is the mean of the OFF frames
        R = (A - mean(A(:, :, ops.OFF_idx), 3)) ./ mean(A(:, :, ops.OFF_idx), 3);
        
    case 'z-score'
        R = (A - mean(A, 3)) ./ std(A, 0, 3);
end

% response calculation

switch ops.resp_type
    
    case 'average'
        % method1: averaged response during ON frames - averaged response during OFF frames
        R = mean(R(:, :, ops.ON_idx), 3) - mean(R(:, :, ops.OFF_idx), 3);
        
    case 'max'
        % method2: max response during ON frames - averaged response during OFF frames
        R = max(R(:, :, ops.ON_idx), [], 3) - mean(R(:, :, ops.OFF_idx), 3);
        
    case 'AUC'
        % method3: AUC during ON frames - AUC during OFF frames
        R = sum(R(:, :, ops.ON_idx), 3) - sum(R(:, :, ops.OFF_idx), 3);
end

% thresholding
R(R<=0) = 1e-6;

end






