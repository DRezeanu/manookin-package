function frameValues = regenerate_spatial_noise(noiseClass, parameters, seed)
%REGENERATE_SPATIAL_NOISE  Reconstruct one epoch's frame sequence from
%   the saved parameters + seed, dispatching to the right lab helper
%   based on the noise-class label.
%
%   frameValues = regenerate_spatial_noise(noiseClass, parameters, seed)
%
%   Validates that every parameter the target helper needs is present
%   before dispatching, and raises a descriptive error listing what's
%   missing and what was available if not. Handles contrast scaling and
%   chromatic-class expansion (RGB / BY / yellow / blue) uniformly for
%   helpers that return achromatic frames, so downstream code always
%   sees frames in the same units regardless of noise class.
%
%   Inputs
%   ------
%     noiseClass  Character vector or string, case-insensitive. One of:
%
%         'jittered'  -> manookinlab.util.getJitteredNoiseFrames
%                        (fine grid + per-frame motion jitter)
%         'fast'      -> inline binary rand (FastNoise protocol; no
%                        standalone helper exists in the util package)
%         'spatial'   -> manookinlab.util.getSpatialNoiseFrames_legacy
%                        (binary / gaussian / ternary spatial noise)
%         'pink'      -> manookinlab.util.getPinkNoiseFrames
%                        (spatio-temporal 1/f noise)
%         'binary'    -> manookinlab.util.getBinaryNoiseFrames
%                        (temporal-only binary noise)
%         'gaussian'  -> manookinlab.util.getGaussianNoiseFrames
%                        (temporal-only Gaussian noise)
%         'uniform'   -> manookinlab.util.getUniformNoiseFrames
%                        (temporal-only uniform noise)
%         'ternary'   -> manookinlab.util.getTernaryNoiseFrames
%
%     parameters  containers.Map or struct carrying the epoch's saved
%                 parameters. See PARAM_KEYS below for the keys each
%                 noiseClass requires. Both types are accepted so the
%                 same function works from a Symphony session (where
%                 parameters come from blockParams.map) and from
%                 Python (where matlab.engine marshals dicts into
%                 structs).
%
%     seed        Integer RNG seed used to generate this epoch's
%                 frames. Come from epoch.parameters('seed') in
%                 Symphony, or from the DataJoint Epoch row in
%                 Python.
%
%   Output
%   ------
%     frameValues (numYChecks x numXChecks x numFrames [x 3])
%                 Frame stack. Achromatic classes return a 3-D array;
%                 chromatic classes return a 4-D array with a trailing
%                 color dimension. Contrast scaling has already been
%                 applied (multiplied by parameters.contrast when
%                 present) so callers see final stimulus values.
%
%   PARAM_KEYS -- required parameters per noiseClass
%   -----------------------------------------------
%     jittered : numXStixels, numYStixels, numXChecks, numYChecks,
%                numFrames, stepsPerStixel, frameDwell
%                (optional: contrast, chromaticClass, colorWeights)
%     fast     : numXStixels, numYStixels, numFrames, frameDwell
%                (optional: contrast, chromaticClass, colorWeights)
%     spatial_legacy  : numXChecks, numYChecks, stimTime, noiseClass,
%                chromaticClass  (helper handles chromatic itself)
%     pink     : numXChecks, numYChecks, numFrames, noiseContrast,
%                spatialAmplitude, temporalAmplitude, chromaticClass
%     binary   : numFrames, frameDwell, contrast
%     gaussian : numFrames, frameDwell, stdev
%     uniform  : numFrames, frameDwell, contrast
%     ternary  : (varargin-based; parameters passed through)
%
%   Example -- from a Python matlab.engine caller
%   ---------------------------------------------
%       # In Python:
%       row = df.iloc[0]                              # DJ epoch row
%       params = row['epoch_parameters']              # dict
%       seed   = int(params['seed'])
%       frames = eng.feval('manookinlab.util.regenerate_spatial_noise',
%                          'jittered', params, seed, nargout=1)
%
%   See also: manookinlab.util.getJitteredNoiseFrames,
%             manookinlab.util.getSpatialNoiseFrames_legacy,
%             manookinlab.util.getPinkNoiseFrames,
%             manookinlab.util.getBinaryNoiseFrames,
%             manookinlab.util.getGaussianNoiseFrames,
%             manookinlab.util.getUniformNoiseFrames,
%             manookinlab.util.getTernaryNoiseFrames

    % Normalise the parameters input to a plain struct so we can use
    % dot access uniformly regardless of whether the caller handed us
    % a containers.Map, a struct, or (from Python) a struct-like dict.
    p = to_struct(parameters);

    % Fill in the two defaults that every Symphony protocol also
    % defaults to when the epoch didn't record them.
    p = default_if_missing(p, 'frameDwell',     1);
    p = default_if_missing(p, 'contrast',       1.0);
    p = default_if_missing(p, 'chromaticClass', 'achromatic');
    p = default_if_missing(p, 'frameRate',     60.0);

    % Dispatch to the right underlying helper. Each branch (a) lists
    % the parameters it needs, (b) checks them up-front so we get a
    % clear error naming any missing key rather than a cryptic
    % downstream failure, and (c) calls the helper with positional
    % args in the exact order its signature expects.
    key = lower(char(noiseClass));
    switch key
        case 'jittered'
            required = {'numXStixels','numYStixels','numXChecks', ...
                        'numYChecks','numFrames','stepsPerStixel'};
            check_params(p, required, key);
            frameValues = manookinlab.util.getJitteredNoiseFrames( ...
                p.numXStixels, p.numYStixels, ...
                p.numXChecks,  p.numYChecks, ...
                p.numFrames,   p.stepsPerStixel, ...
                seed,          p.frameDwell);

        case 'fast'
            % FastNoise generates frames inline in the protocol code;
            % no standalone helper exists in util. Reproduce here.
            required = {'numXStixels','numYStixels','numXChecks', ...
                        'numYChecks','numFrames','stepsPerStixel'};
            check_params(p, required, key);
            frameValues = generate_fast_noise(p, seed);

        case 'spatial_legacy'
            % Note: the inner 'noiseClass' parameter (binary/gaussian/
            % ternary) is a per-epoch property distinct from the outer
            % noiseClass switch we're inside of. Symphony's SpatialNoise
            % protocol stores it under 'noiseClass'.
            required = {'numXChecks','numYChecks', 'stimTime', ...
                        'noiseClass','chromaticClass','frameRate'};
            check_params(p, required, key);
            p.numFrames = floor(p.stimTime*1e-3 * p.frameRate / p.frameDwell);
            frameValues = manookinlab.util.getSpatialNoiseFrames_legacy( ...
                p.numXChecks, p.numYChecks, p.numFrames, ...
                p.noiseClass, p.chromaticClass, seed);
            % Move time dimension to third dimension.
            frameValues = permute(frameValues,[2,3,1]);
            frameValues = repmat(frameValues,[1,1,1,3]);

        case 'pink'
            required = {'numXChecks','numYChecks','numFrames', ...
                        'noiseContrast','spatialAmplitude', ...
                        'temporalAmplitude','chromaticClass'};
            check_params(p, required, key);
            frameValues = manookinlab.util.getPinkNoiseFrames( ...
                p.numXChecks, p.numYChecks, p.numFrames, ...
                p.noiseContrast, p.spatialAmplitude, ...
                p.temporalAmplitude, p.chromaticClass, seed);

        case 'binary'
            required = {'numFrames'};
            check_params(p, required, key);
            frameValues = manookinlab.util.getBinaryNoiseFrames( ...
                p.numFrames, p.frameDwell, p.contrast, seed);

        case 'gaussian'
            required = {'numFrames','stdev'};
            check_params(p, required, key);
            frameValues = manookinlab.util.getGaussianNoiseFrames( ...
                p.numFrames, p.frameDwell, p.stdev, seed);

        case 'uniform'
            required = {'numFrames'};
            check_params(p, required, key);
            frameValues = manookinlab.util.getUniformNoiseFrames( ...
                p.numFrames, p.frameDwell, p.contrast, seed);

        case 'ternary'
            % Ternary's varargin signature accepts a struct/kv list
            % directly; forward the whole parameters struct.
            frameValues = manookinlab.util.getTernaryNoiseFrames(p);

        otherwise
            error('regenerate_spatial_noise:badNoiseClass', ...
                ['Unknown noiseClass ''%s''. Supported: ' ...
                 'jittered, fast, spatial, pink, binary, ' ...
                 'gaussian, uniform, ternary.'], noiseClass);
    end

    % Contrast scaling for helpers that return unit-magnitude output.
    % Spatial / pink / binary / gaussian / uniform / ternary already
    % apply their own scaling internally, so they're skipped here.
    if any(strcmp(key, {'jittered','fast'}))
        frameValues = p.contrast * frameValues;
    end

    % Chromatic-class expansion for helpers that return achromatic
    % frames only. Spatial / pink handle color themselves; ternary
    % has its own convention; binary / gaussian / uniform are
    % temporal-only. Only jittered / fast need expansion.
    if any(strcmp(key, {'jittered','fast'}))
        frameValues = apply_chromatic(frameValues, p, key, seed);
    end
end


% =====================================================================
% Local helpers
% =====================================================================

function p = to_struct(parameters)
%TO_STRUCT  Normalise containers.Map / struct / MATLAB-Engine dict to
%   a plain struct with dot access. Keys are sanitised through
%   makeValidName so keys with unusual characters don't break struct
%   field access.
    if isa(parameters, 'containers.Map')
        ks = parameters.keys;
        p = struct();
        for i = 1:length(ks)
            p.(matlab.lang.makeValidName(ks{i})) = parameters(ks{i});
        end
    elseif isstruct(parameters)
        p = parameters;
    else
        error('regenerate_spatial_noise:badParameters', ...
            'parameters must be a containers.Map or struct (got %s).', ...
            class(parameters));
    end
end


function p = default_if_missing(p, name, value)
%DEFAULT_IF_MISSING  Return p with p.(name) = value only if the field
%   is absent or empty. Leaves existing values untouched.
    if ~isfield(p, name) || isempty(p.(name))
        p.(name) = value;
    end
end


function check_params(p, required, noiseClass)
%CHECK_PARAMS  Verify each required key is present and non-empty in p.
%   Raises a single error listing every missing key plus the fields
%   that WERE available, so the caller can spot naming mismatches
%   (e.g. camelCase vs snake_case) immediately.
    missing = {};
    for i = 1:length(required)
        k = required{i};
        if ~isfield(p, k) || isempty(p.(k))
            missing{end+1} = k; %#ok<AGROW>
        end
    end
    if ~isempty(missing)
        error('regenerate_spatial_noise:missingParam', ...
            ['Missing parameter(s) for noiseClass=%s: %s.\n' ...
             'Available: %s.'], ...
            noiseClass, strjoin(missing, ', '), ...
            strjoin(fieldnames(p), ', '));
    end
end


function frameValues = generate_fast_noise(p, seed)
%GENERATE_FAST_NOISE  Reproduce the inline binary-rand frame generation
%   from the FastNoise protocol (no standalone helper exists).
    noiseStream = RandStream('mt19937ar', 'Seed', seed);
    nEffective  = ceil(p.numFrames / p.frameDwell);
    frameValues = 2 * double( ...
        noiseStream.rand(p.numYStixels, p.numXStixels, nEffective) > 0.5 ...
    ) - 1;
    % Repeat each frame frameDwell times back up to numFrames.
    if p.frameDwell > 1
        frameValues = repelem(frameValues, 1, 1, p.frameDwell);
        frameValues = frameValues(:, :, 1:p.numFrames);
    end
end


function frameValues = apply_chromatic(frameValues, p, noiseClass, seed)
%APPLY_CHROMATIC  Expand achromatic frames along a colour dimension
%   according to p.chromaticClass. Achromatic passthrough leaves the
%   3-D stack alone.
    if strcmpi(p.chromaticClass, 'achromatic')
        return
    end

    switch upper(p.chromaticClass)
        case 'BY'
            % Yellow (R,G) uses seed, Blue uses seed+1. Regenerate the
            % blue channel with an independent noise stream so the two
            % channels are decorrelated but reproducible.
            blue = regenerate_helper_for_blue(p, noiseClass, seed + 1);
            out  = zeros(size(frameValues, 1), size(frameValues, 2), ...
                         size(frameValues, 3), 3);
            out(:, :, :, 1) = frameValues;   % R = yellow
            out(:, :, :, 2) = frameValues;   % G = yellow
            out(:, :, :, 3) = blue;
            frameValues = out;

        otherwise
            % RGB / yellow / blue: expand achromatic frames along a
            % colour axis, scale by colorWeights if present, then
            % suppress the appropriate channel(s) for yellow / blue.
            expanded = repmat(frameValues, [1, 1, 1, 3]);
            if isfield(p, 'colorWeights') && ~isempty(p.colorWeights)
                cw = p.colorWeights;
                for k = 1:3
                    expanded(:, :, :, k) = cw(k) * expanded(:, :, :, k);
                end
            end
            switch lower(p.chromaticClass)
                case 'yellow'
                    expanded(:, :, :, 3) = -1;
                case 'blue'
                    expanded(:, :, :, 1:2) = -1;
            end
            frameValues = expanded;
    end
end


function frames = regenerate_helper_for_blue(p, noiseClass, seed)
%REGENERATE_HELPER_FOR_BLUE  Re-invoke the same underlying frame helper
%   used by the top-level dispatch, but with seed+1, to produce the
%   blue-channel stack for the 'BY' chromatic class.
    switch noiseClass
        case 'jittered'
            frames = manookinlab.util.getJitteredNoiseFrames( ...
                p.numXStixels, p.numYStixels, ...
                p.numXChecks,  p.numYChecks, ...
                p.numFrames,   p.stepsPerStixel, ...
                seed,          p.frameDwell);
        case 'fast'
            frames = generate_fast_noise(p, seed);
        otherwise
            error('regenerate_spatial_noise:BYNotImplemented', ...
                'BY chromaticClass regeneration is only defined for jittered / fast noise.');
    end
    frames = p.contrast * frames;
end
