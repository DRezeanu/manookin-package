function sta_tmp = analyzeSpatialNoise(response, frame_times, noiseClass, parameters, seed, varargin)


validNoiseTypes = {'jittered','fast','spatial','spatial_legacy','pink','binary','gaussian','uniform','ternary'};

ip = inputParser();
addRequired(ip, 'response', @(x) isnumeric(x) && isscalar(x));
addRequired(ip, 'frame_times', @(x) isnumeric(x));
addRequired(ip, 'parameters',@(x) isa(x, 'containers.Map'));
addRequired(ip, 'seed', @(x) isfloat(x) && isscalar(x));
addParameter(ip,'noiseClass','jittered',...
     @(x) any(validatestring(x,validNoiseTypes)));
addParameter(ip, 'method', 'fft', @(x)isstring(x));
addParameter(ip, 'timeBins', 31, @(x) isfloat(x) && isscalar(x));
addParameter(ip, 'binsPerFrame', 1, @(x) isfloat(x) && isscalar(x));
addParameter(ip, 'responseType', 'spikeTimes', @(x) isstring(x))
parse(ip, response, frame_times, noiseClass, parameters, seed, varargin{:});
% 
ops = ip.Results;

% Regenerate the frame sequence.
frameValues = manookinlab.util.regenerate_spatial_noise(ops.noiseClass, ops.parameters, ops.seed);

% Upsample frames, if needed.
if ops.binsPerFrame > 1
    frameValues = up_frames(frameValues, ops.binsPerFrame);
    f_times = [];
    for jj = 1:length(frame_times)-1
        f_tmp = linspace(frame_times(jj),frame_times(jj+1),binsPerFrame+1);
        f_times = [f_times, f_tmp(1:end-1)]; %#ok<AGROW>
    end
    f_times = [f_times, frame_times(end)]; 
    frame_times = f_times;
end
frameValues = frameValues(:,:,1:length(frame_times)-1,:);

switch ops.method
    case 'fft'
        sta_tmp = compute_sta_fft(ops.response, frameValues, ops.timeBins);
end
end

function upFrames = up_frames(frames, multiple)
% 
% frames 
%   4-D: [x,y,t,color]
% multiple: integer multiple (>1)

if multiple <= 1
    upFrames = frames;
    return;
end

% Determine the dimensions.
n = ndims(frames);

if n == 2
    upFrames = zeros(size(frames,1), size(frames,2)*multiple);
    
    for frameIndex = 1 : size(frames,2)*multiple
        fIndex = ceil(frameIndex / multiple);
        upFrames(:,frameIndex) = frames(:,fIndex);
    end
    
elseif n == 3
    upFrames = zeros(size(frames,1), size(frames,2), size(frames,3)*multiple);
    
    for frameIndex = 1 : size(frames,3)*multiple
        fIndex = ceil(frameIndex / multiple);
        upFrames(:,:,frameIndex) = frames(:,:,fIndex);
    end
elseif n == 4
    upFrames = zeros(size(frames,1), size(frames,2), size(frames,3)*multiple, size(frames,4));
    
    for frameIndex = 1 : size(frames,3)*multiple
        fIndex = ceil(frameIndex / multiple);
        upFrames(:,:,frameIndex,:) = frames(:,:,fIndex,:);
    end
else
    error('Number of dimensions must be either 2, 3 or 4!');
end
end


function sta_tmp = compute_sta_fft(response, frameValues, nkt)
sta_tmp = zeros(size(frameValues,1),size(frameValues,2),nkt,3);
for jj = 1 : size(frameValues,1)
    for kk = 1 : size(frameValues,2)
        for mm = 1:3
            stim_tmp = squeeze(frameValues(jj,kk,:,mm))';
            foo = ifft(fft(response) .* conj(fft(stim_tmp)));
            sta_tmp(jj,kk,:,mm) = foo(1:nkt);
        end
    end
end
end


