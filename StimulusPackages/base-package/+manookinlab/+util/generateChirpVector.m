function stim = generateChirpVector(sampleRate, varargin)
%GENERATECHIRPVECTOR  Build the intensity time-course of the Chirp
%   stimulus (Baden et al. 2016 style: pre-time + light step + dark
%   step + frequency sweep + contrast sweep + tail-time) at a given
%   sample rate.
%
%   stim = generateChirpVector(sampleRate, 'Name', Value, ...)
%
%   The returned vector is in absolute intensity units in [0, 1]
%   centred on ``backgroundIntensity``. Its length is exactly
%   ``round((preTime + tailTime + interTime*3 + stepTime*2 +
%   frequencyTime + contrastTime) * 1e-3 * sampleRate)`` samples;
%   time inputs are in MILLISECONDS.
%
%   Segment order, in wall-clock samples:
%
%     preTime          -- backgroundIntensity only
%     stepTime         -- + stepContrast (light step)
%     interTime        -- backgroundIntensity
%     stepTime         -- - stepContrast (dark step)
%     interTime        -- backgroundIntensity
%     frequencyTime    -- sinusoid, freq swept linearly from
%                         frequencyMin to frequencyMax at
%                         frequencyContrast
%     interTime        -- backgroundIntensity
%     contrastTime     -- sinusoid at contrastFrequency, contrast
%                         swept linearly from contrastMin to
%                         contrastMax
%     tailTime         -- backgroundIntensity only
%
%   Required input
%   --------------
%     sampleRate   Scalar acquisition rate in Hz. All time-based
%                  options below are then converted to sample counts
%                  via round(time_ms * 1e-3 * sampleRate).
%
%   Name-value options (defaults in parens; every option is a scalar)
%   -----------------------------------------------------------------
%     'preTime'             (500 ms)    Pre-stimulus grey window.
%     'tailTime'            (500 ms)    Post-stimulus grey window.
%     'stepTime'            (500 ms)    Duration of each step
%                                        (light and dark).
%     'interTime'           (500 ms)    Inter-segment grey gap.
%     'frequencyTime'       (15000 ms)  Duration of the frequency
%                                        sweep segment.
%     'contrastTime'        (8000 ms)   Duration of the contrast
%                                        sweep segment.
%     'stepContrast'        (1.0)       Weber contrast of the two
%                                        steps. Multiplied by
%                                        backgroundIntensity and
%                                        added / subtracted from it.
%     'frequencyContrast'   (1.0)       Contrast (peak amplitude) of
%                                        the frequency-sweep sinusoid.
%     'frequencyMin'        (0.0 Hz)    Starting frequency of the
%                                        linear sweep.
%     'frequencyMax'        (10.0 Hz)   Ending frequency of the
%                                        linear sweep.
%     'contrastMin'         (0.02)      Starting contrast of the
%                                        contrast-sweep sinusoid.
%     'contrastMax'         (1.0)       Ending contrast.
%     'contrastFrequency'   (2.0 Hz)    Fixed carrier frequency of
%                                        the contrast-sweep sinusoid.
%     'backgroundIntensity' (0.5)       Background luminance level
%                                        (0-1). All segments are
%                                        expressed relative to this.
%
%   Output
%   ------
%     stim         (1 x nSamples) intensity vector in [0, 1], with
%                  values centred on backgroundIntensity. nSamples =
%                  round(totalTimeMs * 1e-3 * sampleRate).
%
%   Example
%   -------
%     % Default 25-second chirp at 10 kHz -> 250000-sample vector.
%     stim = generateChirpVector(10000);
%
%     % Custom: shorter sweep, higher contrast, 20 kHz sampling.
%     stim = generateChirpVector(20000, ...
%                                'frequencyTime', 10000, ...
%                                'contrastMax',   0.8);
%
%   See also: manookinlab.util.getJitteredNoiseFrames,
%             manookinlab.util.regenerate_spatial_noise

ip = inputParser();
addRequired(ip, 'sampleRate', @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'preTime', 500, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'tailTime', 500, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'stepTime', 500, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'interTime', 500, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'frequencyTime', 15000, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'contrastTime', 8000, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'stepContrast', 1.0, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'frequencyContrast', 1.0, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'frequencyMin', 0.0, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'frequencyMax', 10.0, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'contrastMin', 0.02, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'contrastMax', 1.0, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'contrastFrequency', 2.0, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'backgroundIntensity', 0.5, @(x) isnumeric(x) && isscalar(x));
parse(ip, sampleRate, varargin{:});
% 
ops = ip.Results;

timeToPts = @(t)(round(t * 1e-3 * ops.sampleRate));
ptsToTime = @(p)(p / ops.sampleRate); % in sec

stimTime = ops.interTime*3 + ops.stepTime*2 + ops.frequencyTime + ops.contrastTime;
           
totTime = ops.preTime + stimTime + ops.tailTime;
totPts = timeToPts(totTime);
stim = ones(1, totPts);
stim(1:totPts) = ops.backgroundIntensity;

prePts = timeToPts(ops.preTime);
interPts = timeToPts(ops.interTime);
stepPts = timeToPts(ops.stepTime);
freqPts = timeToPts(ops.frequencyTime);
contrastPts = timeToPts(ops.contrastTime);

frequencyDelta = (ops.frequencyMax - ops.frequencyMin)/freqPts/2; % not sure why factor of 2 needed but gets frequencies right
contrastDelta = (ops.contrastMax - ops.contrastMin)/contrastPts;

% increment and decrement steps
stim(prePts+(1:stepPts)) = stim(prePts+(1:stepPts)) + ops.stepContrast * ops.backgroundIntensity;
stim(prePts+interPts+stepPts+(1:stepPts)) = stim(prePts+interPts+stepPts+(1:stepPts)) - ops.stepContrast * ops.backgroundIntensity;

% frequency sweep
for t = 1:freqPts
    stim(t + prePts+interPts*2+stepPts*2) = ops.frequencyContrast*ops.backgroundIntensity*sin(2*pi*ptsToTime(t)*(ops.frequencyMin+frequencyDelta*t)) + ops.backgroundIntensity;
end

% contrast sweep
for t = 1:contrastPts
    stim(t + prePts+interPts*3+stepPts*2+freqPts) = (ops.contrastMin+t*contrastDelta)*ops.backgroundIntensity*sin(2*pi*ptsToTime(t)*ops.contrastFrequency) + ops.backgroundIntensity;
end     