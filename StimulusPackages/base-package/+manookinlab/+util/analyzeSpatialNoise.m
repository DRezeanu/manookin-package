function [strf, spaceFilter] = analyzeSpatialNoise(response, varargin)
% ANALYZESPATIALNOISE Compute the spatiotemporal receptive field for the
% SpatialNoise (FastNoise) stimulus by reverse correlation.
%
%   [strf, spaceFilter] = analyzeSpatialNoise(response, ...)
%
% Inputs:
%   response : response binned at the frame rate (one bin per stimulus
%              frame, e.g. spike counts from BinSpikeRate).
%
% Named parameters:
%   numXStixels, numYStixels : stixel grid dimensions (epoch parameters).
%   numXChecks, numYChecks   : check dimensions (epoch parameters).
%   chromaticClass  : 'achromatic', 'BY', 'RGB', ... (default 'achromatic')
%   numFrames       : total number of stimulus frames (default numel(response)).
%   stepsPerStixel  : jitter steps per stixel (default 1).
%   seed            : random seed for the unique sequence.
%   frameDwell      : monitor frames per stimulus frame (default 1).
%   uniqueFrames    : number of frames of unique (seeded) noise. Frames
%                     beyond this come from the repeating sequence.
%                     (default: all frames unique)
%   repeatFrames    : number of repeating-sequence frames (default 0).
%   repeatingSeed   : seed of the repeating sequence (default 1).
%   gaussianFilter  : whether the stixels were blurred (default false).
%   filterSdStixels : Gaussian filter SD in stixels (default 1.0).
%   frameRate       : monitor frame rate in Hz (default 60).
%   filterFrames    : number of time lags in the filter
%                     (default floor(frameRate/2)).
%
% Outputs:
%   strf        : (numYChecks x numXChecks x filterFrames x nChannels)
%                 spatiotemporal receptive field. Channels are: achromatic
%                 and single-color classes -> 1; 'BY' -> [yellow, blue];
%                 'RGB' -> [red, green, blue].
%   spaceFilter : (numYChecks x numXChecks x nChannels) spatial RF averaged
%                 over the dominant time lobe.

ip = inputParser();
ip.addParameter('numXStixels', [], @(x)isfloat(x));
ip.addParameter('numYStixels', [], @(x)isfloat(x));
ip.addParameter('numXChecks', [], @(x)isfloat(x));
ip.addParameter('numYChecks', [], @(x)isfloat(x));
ip.addParameter('chromaticClass', 'achromatic', @(x)ischar(x));
ip.addParameter('numFrames', [], @(x)isfloat(x));
ip.addParameter('stepsPerStixel', 1, @(x)isfloat(x));
ip.addParameter('seed', 1, @(x)isfloat(x));
ip.addParameter('frameDwell', 1, @(x)isfloat(x));
ip.addParameter('uniqueFrames', [], @(x)isfloat(x));
ip.addParameter('repeatFrames', 0, @(x)isfloat(x));
ip.addParameter('repeatingSeed', 1, @(x)isfloat(x));
ip.addParameter('gaussianFilter', false, @(x)islogical(x) || isfloat(x));
ip.addParameter('filterSdStixels', 1.0, @(x)isfloat(x));
ip.addParameter('frameRate', 60.0, @(x)isfloat(x));
ip.addParameter('filterFrames', [], @(x)isfloat(x));
ip.parse(varargin{:});

numXStixels = ip.Results.numXStixels;
numYStixels = ip.Results.numYStixels;
numXChecks = ip.Results.numXChecks;
numYChecks = ip.Results.numYChecks;
chromaticClass = ip.Results.chromaticClass;
numFrames = ip.Results.numFrames;
stepsPerStixel = double(ip.Results.stepsPerStixel);
seed = ip.Results.seed;
frameDwell = double(ip.Results.frameDwell);
uniqueFrames = ip.Results.uniqueFrames;
repeatFrames = ip.Results.repeatFrames;
repeatingSeed = ip.Results.repeatingSeed;
gaussianFilter = logical(ip.Results.gaussianFilter);
filterSdStixels = ip.Results.filterSdStixels;
frameRate = ip.Results.frameRate;
filterFrames = ip.Results.filterFrames;

y = response(:);

if isempty(numFrames)
    numFrames = numel(y);
end
if isempty(filterFrames)
    filterFrames = floor(frameRate * 0.5);
end

% Regenerate the stimulus: (numYChecks x numXChecks x t x 3) contrasts.
if ~isempty(uniqueFrames) && repeatFrames > 0 && uniqueFrames < numFrames
    % Unique (seeded) segment followed by the repeating segment.
    stim = manookinlab.util.getSpatialNoiseFrames(numXStixels, numYStixels, ...
        numXChecks, numYChecks, chromaticClass, uniqueFrames, ...
        stepsPerStixel, seed, frameDwell, gaussianFilter, filterSdStixels);
    stim = stim(:, :, 1 : min(end, uniqueFrames), :);
    stimRep = manookinlab.util.getSpatialNoiseFrames(numXStixels, numYStixels, ...
        numXChecks, numYChecks, chromaticClass, numFrames - uniqueFrames, ...
        stepsPerStixel, repeatingSeed, frameDwell, gaussianFilter, filterSdStixels);
    stimRep = stimRep(:, :, 1 : min(end, numFrames - uniqueFrames), :);
    stim = cat(3, stim, stimRep);
else
    stim = manookinlab.util.getSpatialNoiseFrames(numXStixels, numYStixels, ...
        numXChecks, numYChecks, chromaticClass, numFrames, ...
        stepsPerStixel, seed, frameDwell, gaussianFilter, filterSdStixels);
end

% Match the response and stimulus lengths.
nT = min(numel(y), size(stim, 3));
y = y(1 : nT);
stim = stim(:, :, 1 : nT, :);

% Zero out the first second while the cell is adapting to the stimulus,
% and the last 15 (padding) frames.
adaptFrames = min(floor(frameRate), nT);
y(1 : adaptFrames) = 0;
stim(:, :, 1 : adaptFrames, :) = 0;
if nT > 15
    y(end - 14 : end) = 0;
    stim(:, :, end - 14 : end, :) = 0;
end

% Pick the informative color channels.
if strcmpi(chromaticClass, 'BY')
    channels = [1, 3]; % yellow (R=G), blue
elseif strcmpi(chromaticClass, 'RGB')
    channels = [1, 2, 3];
else
    channels = 1; % all three channels are identical
end
nChannels = numel(channels);

% Reverse correlation via FFT for each check and channel.
strf = zeros(numYChecks, numXChecks, filterFrames, nChannels);
fftY = fft([y; zeros(60, 1)]);
for c = 1 : nChannels
    for m = 1 : numYChecks
        for n = 1 : numXChecks
            s = double(squeeze(stim(m, n, :, channels(c))));
            tmp = ifft(fftY .* conj(fft([s; zeros(60, 1)])));
            strf(m, n, :, c) = tmp(1 : filterFrames);
        end
    end
end

% Spatial RF averaged over the dominant time lobe.
lobePts = 2 : min(4, filterFrames);
spaceFilter = reshape(mean(strf(:, :, lobePts, :), 3), ...
    [numYChecks, numXChecks, nChannels]);

end
