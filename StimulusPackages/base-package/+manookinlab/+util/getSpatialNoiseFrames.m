function stimulus = getSpatialNoiseFrames(numXStixels, numYStixels, numXChecks, numYChecks, chromaticClass, numFrames, stepsPerStixel, seed, frameDwell, gaussianFilter, filterSdStixels)
% GETSPATIALNOISEFRAMES Regenerate the frame sequence for the SpatialNoise
% (FastNoise) stimulus. MATLAB port of the Python getFastNoiseFrames.
%
%   stimulus = getSpatialNoiseFrames(numXStixels, numYStixels, numXChecks,
%       numYChecks, chromaticClass, numFrames, stepsPerStixel, seed,
%       [frameDwell], [gaussianFilter], [filterSdStixels])
%
% Parameters:
%   numXStixels    : number of stixels in the x direction.
%   numYStixels    : number of stixels in the y direction.
%   numXChecks     : number of checks in the x direction.
%   numYChecks     : number of checks in the y direction.
%   chromaticClass : chromatic class of the stimulus ('achromatic','BY','RGB',...)
%   numFrames      : number of frames in the stimulus.
%   stepsPerStixel : number of jitter steps per stixel.
%   seed           : seed for the random number generator.
%   frameDwell     : number of monitor frames to dwell on each frame (default 1).
%   gaussianFilter : whether to blur the stixels (default false).
%   filterSdStixels: Gaussian filter SD in stixels (default 1.0).
%
% Returns:
%   stimulus : 4-D single array of contrasts in [-1,1] with size
%              (numYChecks x numXChecks x nFrames x 3), i.e. (y, x, t, color).
%              Note the Python original returns (t, y, x, color); the pixel
%              values are identical, only the dimension order differs.
%
% The random sequence matches numpy.random.seed(seed)/rand draw-for-draw:
% both MATLAB's mt19937ar and numpy's legacy RandomState initialize with
% MT19937 init_genrand(seed) and generate 53-bit doubles (genrand_res53).
% CAVEAT: the single exception is seed == 0, which MATLAB silently remaps
% to the MT19937 default seed (5489) while numpy uses 0 literally, so a
% seed of exactly 0 will NOT reproduce the Python sequence.

if nargin < 9 || isempty(frameDwell)
    frameDwell = 1;
end
if nargin < 10 || isempty(gaussianFilter)
    gaussianFilter = false;
end
if nargin < 11 || isempty(filterSdStixels)
    filterSdStixels = 1.0;
end

frameDwell = double(frameDwell);
stepsPerStixel = double(stepsPerStixel);

% Time expansion factor for the chromatic classes.
if strcmpi(chromaticClass, 'BY')
    tfactor = 2;
elseif strcmpi(chromaticClass, 'RGB')
    tfactor = 3;
else % Black/white checks
    tfactor = 1;
end

% Get the size of the time dimension; expands for BY/RGB.
tsize = ceil(numFrames * tfactor / frameDwell);
if tfactor == 2 && mod(tsize, 2) ~= 0
    tsize = tsize + 1;
end

% Seed the random number generator.
noiseStream = RandStream('mt19937ar', 'Seed', seed);

% Generate the random grid of stixels.
% numpy draws rand(tsize, nStixels) in C (row-major) order: all stixels of
% frame 0, then frame 1, ... MATLAB fills column-major, so draw an
% (nStixels x tsize) matrix: column t holds the draws for frame t. Within
% a frame, numpy's reshape+transpose maps draw j to x = floor(j/numYStixels),
% y = mod(j, numYStixels), which is exactly MATLAB's column-major reshape
% to (numYStixels x numXStixels).
nStixels = numXStixels * numYStixels;
gridValues = reshape(rand(noiseStream, nStixels, tsize), ...
    [numYStixels, numXStixels, tsize]);
gridValues = single(2 * round(gridValues) - 1); % Convert to contrast

% Filter the stixels if indicated (scipy.ndimage.gaussian_filter,
% mode='wrap', truncate=4.0 -> radius = floor(4*sd + 0.5)).
if gaussianFilter
    r = floor(4 * filterSdStixels + 0.5);
    kx = (-r : r) / filterSdStixels;
    kernel = exp(-0.5 * kx.^2);
    kernel = kernel / sum(kernel);
    ny = size(gridValues, 1);
    nx = size(gridValues, 2);
    rowIdx = mod((1 - r : ny + r) - 1, ny) + 1; % circular padding
    colIdx = mod((1 - r : nx + r) - 1, nx) + 1;
    for t = 1 : tsize
        fp = gridValues(rowIdx, colIdx, t);
        f2 = conv2(kernel(:), kernel(:)', fp, 'valid');
        % numpy std is the population std (normalized by N).
        gridValues(:,:,t) = 0.5 * f2 / std(f2(:), 1);
    end
    gridValues(gridValues > 1.0) = 1.0;
    gridValues(gridValues < -1.0) = -1.0;
end

% Translate to the full grid of checks.
fullGrid = repelem(gridValues, stepsPerStixel, stepsPerStixel, 1);
fullY = size(fullGrid, 1);
fullX = size(fullGrid, 2);

% Generate the motion trajectory of the larger stixels. numpy re-seeds the
% same generator; draws are row-major pairs (x, y) per frame.
positionStream = RandStream('mt19937ar', 'Seed', seed);
posDraws = rand(positionStream, 2, tsize);
xSteps = round((stepsPerStixel - 1) * posDraws(1, :));
ySteps = round((stepsPerStixel - 1) * posDraws(2, :));

% Compute crop amounts so fullGrid is centered on the frame canvas.
% numpy rounds half to even (banker's rounding).
cropY = round_half_even((fullY - numYChecks) / 2);
cropX = round_half_even((fullX - numXChecks) / 2);

% Get the frame values for the finer grid. Indices below are 0-based to
% mirror the Python source; +1 applied at indexing time.
frameValues = zeros(numYChecks, numXChecks, tsize, 'single');
for k = 0 : tsize - 1
    xOff = xSteps(floor(k / tfactor) + 1);
    yOff = ySteps(floor(k / tfactor) + 1);
    % Y jitter moves up and adds to crop from the top.
    yOff = yOff + cropY;
    % X jitter moves right, so subtracts from crop from the left.
    xOff = cropX - xOff;

    % X start and ends
    gx0 = xOff;         gx1 = xOff + numXChecks;
    fx0 = 0;            fx1 = numXChecks;
    % Y start and ends
    gy0 = yOff;         gy1 = yOff + numYChecks;
    fy0 = 0;            fy1 = numYChecks;

    % If the grid starts left of the canvas, leave the left edge gray.
    if gx0 < 0
        fx0 = -gx0;
        gx0 = 0;
        gx1 = numXChecks - fx0;
    end
    % If the grid ends beyond the canvas, leave the bottom/right edge gray.
    if gy1 > fullY
        nBottomGray = gy1 - fullY;
        fy1 = numYChecks - nBottomGray;
        gy1 = fullY;
    end
    if gx1 > fullX
        nRightGray = gx1 - fullX;
        fx1 = numXChecks - nRightGray;
        gx1 = fullX;
    end

    frameValues(fy0 + 1 : fy1, fx0 + 1 : fx1, k + 1) = ...
        fullGrid(gy0 + 1 : gy1, gx0 + 1 : gx1, k + 1);
end

% Sort the pixel values into the proper color channels. (y, x, t, color)
if strcmpi(chromaticClass, 'BY')
    rg = frameValues(:, :, 1 : 2 : end);
    b  = frameValues(:, :, 2 : 2 : end);
    nF = min(size(rg, 3), size(b, 3));
    stimulus = zeros(numYChecks, numXChecks, nF, 3, 'single');
    stimulus(:, :, :, 1) = rg(:, :, 1 : nF);
    stimulus(:, :, :, 2) = rg(:, :, 1 : nF);
    stimulus(:, :, :, 3) = b(:, :, 1 : nF);
elseif strcmpi(chromaticClass, 'RGB')
    rr = frameValues(:, :, 1 : 3 : end);
    gg = frameValues(:, :, 2 : 3 : end);
    bb = frameValues(:, :, 3 : 3 : end);
    nF = min([size(rr, 3), size(gg, 3), size(bb, 3)]);
    stimulus = zeros(numYChecks, numXChecks, nF, 3, 'single');
    stimulus(:, :, :, 1) = rr(:, :, 1 : nF);
    stimulus(:, :, :, 2) = gg(:, :, 1 : nF);
    stimulus(:, :, :, 3) = bb(:, :, 1 : nF);
else % Black/white checks
    stimulus = repmat(reshape(frameValues, ...
        [numYChecks, numXChecks, tsize, 1]), [1, 1, 1, 3]);
end

% Deal with the frame dwell.
if frameDwell > 1
    numStimFrames = size(stimulus, 3);
    nOut = numStimFrames * frameDwell;
    idx = min(floor((0 : nOut - 1) / frameDwell), numStimFrames - 1) + 1;
    stimulus = stimulus(:, :, idx, :);
end

end

function r = round_half_even(x)
% Round half to even (numpy/banker's rounding) for scalar x.
f = floor(x);
if (x - f) == 0.5
    if mod(f, 2) == 0
        r = f;
    else
        r = f + 1;
    end
else
    r = round(x);
end
end
