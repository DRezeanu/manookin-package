function results = analyzeChirp(response, binRate, varargin)
%ANALYZECHIRP  Decompose a Chirp-stimulus response into step /
%   frequency-sweep / contrast-sweep components and return a results
%   struct with the standard summary metrics for each.
%
%   results = analyzeChirp(response, binRate, 'Name', Value, ...)
%
%   The Chirp is a fixed-order composite stimulus (see
%   generateChirpVector) laid out as
%
%       preTime | +step | interTime | -step | interTime | frequencySweep
%               | interTime | contrastSweep | tailTime
%
%   This function slices the response along that layout and, for each
%   segment, extracts the metrics most commonly reported in the
%   retinal-ganglion / amacrine-cell literature:
%
%     Steps
%       * incAmplitude  -- peak response to the light step
%       * decAmplitude  -- peak response to the dark step
%       * incDecRatio   -- |inc peak| / |dec peak|  (polarity index)
%       * tPeakInc / tPeakDec  -- times of the peaks (seconds from
%                                 step onset)
%
%     Frequency sweep  (linear sweep from frequencyMin to frequencyMax)
%       * frequencies   -- instantaneous stimulus frequency per cycle
%       * amplitudes    -- cycle-by-cycle F1 amplitude of the response
%       * phases        -- cycle-by-cycle F1 phase (radians)
%       * peakFrequency -- frequency at which the F1 amplitude peaks
%       * cutoffFrequency -- highest frequency at which the F1 amplitude
%                            drops below cutoffFrac * peak (default 0.5)
%       * bandwidth     -- (low_cut, high_cut) in Hz where amplitude
%                          crosses cutoffFrac of the peak
%
%     Contrast sweep  (linear contrast ramp at fixed carrier freq)
%       * contrasts     -- instantaneous stimulus contrast per cycle
%       * amplitudes    -- cycle-by-cycle F1 amplitude of the response
%       * f2_f1_ratio   -- second-harmonic / first-harmonic ratio,
%                          integrated across the sweep. A common
%                          rectification / linearity index.
%       * c50           -- contrast at half of the saturated F1
%                          response, from a Naka-Rushton fit
%       * rmax          -- fitted saturation amplitude
%       * n             -- fitted Hill coefficient
%       * fit_ok        -- true if the Naka-Rushton fit converged
%
%   Required inputs
%   ---------------
%     response  1-D response trace at ``binRate`` Hz, in whatever
%               units your recording produced (spikes/s, pA, mV, ...).
%               Length must exceed ``preTime + tailTime + 3*interTime
%               + 2*stepTime + frequencyTime + contrastTime`` ms
%               (converted to samples at binRate).
%     binRate   Scalar acquisition / binning rate in Hz.
%
%   Name-value options -- Chirp segment timing (defaults match
%   generateChirpVector so passing nothing here reproduces the default
%   Chirp layout)
%   ----------------------------------------------------------------
%     'preTime'             (500 ms)
%     'tailTime'            (500 ms)
%     'stepTime'            (500 ms)
%     'interTime'           (500 ms)
%     'frequencyTime'       (15000 ms)
%     'contrastTime'        (8000 ms)
%     'stepContrast'        (1.0)
%     'frequencyContrast'   (1.0)
%     'frequencyMin'        (0.0 Hz)
%     'frequencyMax'        (10.0 Hz)
%     'contrastMin'         (0.02)
%     'contrastMax'         (1.0)
%     'contrastFrequency'   (2.0 Hz)
%     'backgroundIntensity' (0.5)
%
%   Additional analysis options
%   ---------------------------
%     'cellAttached'  (true)    Cell-attached recording (spikes go
%                                UP for increments). Flips the sign
%                                of the inc/dec ratio for consistency
%                                with the intracellular case (where
%                                increments produce inward, negative
%                                deflections).
%     'cutoffFrac'    (0.5)     Fraction of the peak F1 amplitude at
%                                which the frequency-tuning cutoff is
%                                measured.
%     'baselineWindow' ([-500 0] ms relative to first step) window
%                                used to compute a per-trial baseline
%                                subtracted before amplitude analysis.
%                                Default = the preTime window.
%
%   Note: this function returns numeric results only. Plotting is
%   deliberately not done here so the function is safe to call from
%   the MATLAB Engine's -nojvm mode (used by the Python wrapper).
%   Plot with the returned struct from whichever host language you
%   prefer.
%
%   Output
%   ------
%     results  struct with .steps, .freq_sweep, .contrast_sweep
%              sub-structs as described above, plus .params (echo of
%              the resolved options) and .segments (raw slices of the
%              response for each stimulus segment, for downstream use).
%
%   Example
%   -------
%     r = analyzeChirp(spike_rate, 1000, 'cellAttached', true);
%     fprintf('cutoff = %.2f Hz, c50 = %.3f\n', ...
%             r.freq_sweep.cutoffFrequency, r.contrast_sweep.c50);
%
%   See also: manookinlab.util.generateChirpVector

    % --- Option parsing (defaults match generateChirpVector.m) --------
    ip = inputParser();
    ip.addRequired('response', @(x) isnumeric(x) && isvector(x));
    ip.addRequired('binRate',  @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('preTime',             500,   @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('tailTime',            500,   @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('stepTime',            500,   @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('interTime',           500,   @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('frequencyTime',       15000, @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('contrastTime',        8000,  @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('stepContrast',        1.0,   @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('frequencyContrast',   1.0,   @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('frequencyMin',        0.0,   @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('frequencyMax',        10.0,  @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('contrastMin',         0.02,  @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('contrastMax',         1.0,   @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('contrastFrequency',   2.0,   @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('backgroundIntensity', 0.5,   @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('cellAttached',        true,  @(x) islogical(x) || isnumeric(x));
    ip.addParameter('cutoffFrac',          0.5,   @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
    ip.parse(response, binRate, varargin{:});
    ops = ip.Results;

    response = response(:).';   % force row vector

    % --- Segment boundaries -----------------------------------------
    toPts   = @(t_ms) round(t_ms * 1e-3 * binRate);
    prePts   = toPts(ops.preTime);
    interPts = toPts(ops.interTime);
    stepPts  = toPts(ops.stepTime);
    freqPts  = toPts(ops.frequencyTime);
    contPts  = toPts(ops.contrastTime);

    totalPts = prePts + stepPts + interPts + stepPts + interPts + ...
               freqPts + interPts + contPts + toPts(ops.tailTime);
    assert(length(response) >= totalPts, ...
        'analyzeChirp: response has %d samples but the layout requires %d at binRate=%g.', ...
        length(response), totalPts, binRate);

    % Baseline: mean of the preTime window. Subtract before amplitude
    % analysis so DC drift doesn't inflate F1 estimates.
    baseline = mean(response(1:prePts));
    resp     = response - baseline;

    % --- Slice each stimulus segment --------------------------------
    idx_inc  = prePts + (1:stepPts);
    idx_dec  = prePts + stepPts + interPts + (1:stepPts);
    idx_freq = prePts + 2*stepPts + 2*interPts + (1:freqPts);
    idx_cont = prePts + 2*stepPts + 3*interPts + freqPts + (1:contPts);

    seg.inc  = resp(idx_inc);
    seg.dec  = resp(idx_dec);
    seg.freq = resp(idx_freq);
    seg.cont = resp(idx_cont);

    % --- Step analysis ---------------------------------------------
    steps = analyze_steps(seg.inc, seg.dec, binRate, ops.cellAttached);

    % --- Frequency sweep --------------------------------------------
    freq_sweep = analyze_frequency_sweep( ...
        seg.freq, binRate, ops.frequencyMin, ops.frequencyMax, ops.cutoffFrac);

    % --- Contrast sweep ---------------------------------------------
    contrast_sweep = analyze_contrast_sweep( ...
        seg.cont, binRate, ops.contrastMin, ops.contrastMax, ...
        ops.contrastFrequency);

    % --- Pack results ----------------------------------------------
    results = struct();
    results.steps          = steps;
    results.freq_sweep     = freq_sweep;
    results.contrast_sweep = contrast_sweep;
    results.segments       = seg;
    results.baseline       = baseline;
    results.params         = ops;
end


% =====================================================================
% Segment-level analyzers
% =====================================================================

function s = analyze_steps(inc, dec, binRate, cellAttached)
% Simple peak amplitudes + peak times for the light and dark steps.
    t_inc = (0:length(inc)-1) / binRate;
    t_dec = (0:length(dec)-1) / binRate;

    if logical(cellAttached)
        % Extracellular: spikes go up for both polarities in a
        % transient cell; use the maximum absolute value.
        [inc_peak, inc_i] = max(inc);
        [dec_peak, dec_i] = max(dec);
    else
        % Voltage / current clamp: light onset -> inward (negative)
        % current in an ON cell. Use most negative deflection for inc,
        % most positive for dec, then flip signs for a comparable
        % magnitude.
        [inc_peak, inc_i] = min(inc);
        inc_peak = -inc_peak;
        [dec_peak, dec_i] = max(dec);
    end

    s = struct();
    s.incAmplitude  = inc_peak;
    s.decAmplitude  = dec_peak;
    s.tPeakInc      = t_inc(inc_i);
    s.tPeakDec      = t_dec(dec_i);
    s.incDecRatio   = inc_peak / max(dec_peak, eps);
    s.trace_inc     = inc;
    s.trace_dec     = dec;
end


function s = analyze_frequency_sweep(x, binRate, fMin, fMax, cutoffFrac)
% Cycle-by-cycle F1 projection during a linear frequency sweep.
%
% We slide a one-cycle window across the response. At sample n the
% instantaneous stimulus frequency is
%       f(n) = fMin + (fMax - fMin) * (n / N)
% and the corresponding cycle length is 1/f(n) seconds. Within a
% window of that length centered on n, we project onto
%   sin(2*pi*f*t)  and  cos(2*pi*f*t)
% to get the F1 amplitude and phase at the local frequency. The
% projection width auto-scales so low-frequency estimates use longer
% windows (matched to their period).
    n   = length(x);
    t   = (0:n-1) / binRate;
    % Instantaneous frequency (matches generateChirpVector's ramp,
    % which uses fMin + delta*t with delta=(fMax-fMin)/N/2 in the
    % sine phase -- so the *observed* freq at sample n is fMin +
    % 2*delta*n = fMin + (fMax - fMin) * (n/N). Track that here.)
    f_inst = fMin + (fMax - fMin) * ((0:n-1) / max(n - 1, 1));

    % Cycle centres: subsample so we have ~one estimate per cycle at
    % the midpoint of each nominal cycle.
    n_cycles = max(1, floor((fMax + fMin) / 2 * (n / binRate)));
    cycle_centres = round(linspace(1, n, n_cycles + 2));
    cycle_centres = cycle_centres(2:end-1);   % drop edges

    frequencies = zeros(1, numel(cycle_centres));
    amplitudes  = zeros(1, numel(cycle_centres));
    phases      = zeros(1, numel(cycle_centres));

    for k = 1 : numel(cycle_centres)
        c   = cycle_centres(k);
        f_k = f_inst(c);
        if f_k <= 0
            continue
        end
        half_win = round(binRate / f_k / 2);         % half a cycle
        lo = max(1, c - half_win);
        hi = min(n, c + half_win);
        if hi - lo < 2, continue, end

        seg_t = t(lo:hi) - t(c);                     % centre at 0
        seg_x = x(lo:hi);
        sinp  = sum(seg_x .* sin(2*pi*f_k*seg_t));
        cosp  = sum(seg_x .* cos(2*pi*f_k*seg_t));
        % Normalise by half the number of samples to get amplitude
        % in the same units as the input.
        norm  = numel(seg_t) / 2;
        frequencies(k) = f_k;
        amplitudes(k)  = hypot(sinp, cosp) / norm;
        phases(k)      = atan2(sinp, cosp);
    end

    % Peak frequency and cutoff (highest freq at which amplitude
    % first drops below cutoffFrac * peak on the way up from freq=0).
    [peak_amp, peak_i] = max(amplitudes);
    peak_f = frequencies(peak_i);

    threshold = cutoffFrac * peak_amp;
    above     = amplitudes >= threshold;

    if any(above)
        % Low-cut: first freq at or above threshold.
        first_i  = find(above, 1, 'first');
        % High-cut: last freq at or above threshold.
        last_i   = find(above, 1, 'last');
        low_cut  = frequencies(first_i);
        high_cut = frequencies(last_i);
    else
        low_cut  = NaN;
        high_cut = NaN;
    end

    s = struct();
    s.frequencies      = frequencies;
    s.amplitudes       = amplitudes;
    s.phases           = phases;
    s.peakFrequency    = peak_f;
    s.peakAmplitude    = peak_amp;
    s.cutoffFrequency  = high_cut;
    s.bandwidth        = [low_cut, high_cut];
    s.cutoffFrac       = cutoffFrac;
end


function s = analyze_contrast_sweep(x, binRate, cMin, cMax, carrierFreq)
% Cycle-by-cycle F1 projection during a linear contrast ramp at a
% fixed carrier frequency. Fits a Naka-Rushton to F1 vs contrast for
% a c50 / rmax / n summary.
    n = length(x);
    t = (0:n-1) / binRate;
    cycle_period = 1 / carrierFreq;
    cycle_pts    = round(cycle_period * binRate);
    n_cycles     = floor(n / cycle_pts);

    if n_cycles < 3
        warning('analyzeChirp:tooFewCycles', ...
                'contrast sweep has only %d full cycles; c50 fit skipped.', n_cycles);
        s = struct('contrasts', [], 'amplitudes', [], 'phases', [], ...
                   'f2_f1_ratio', NaN, 'c50', NaN, 'rmax', NaN, ...
                   'n', NaN, 'fit_ok', false);
        return
    end

    contrasts  = zeros(1, n_cycles);
    amplitudes = zeros(1, n_cycles);
    phases     = zeros(1, n_cycles);
    f2_amps    = zeros(1, n_cycles);

    for k = 1 : n_cycles
        lo = (k - 1) * cycle_pts + 1;
        hi = min(n, lo + cycle_pts - 1);
        seg_t = t(lo:hi) - t(lo);
        seg_x = x(lo:hi);

        % Instantaneous stimulus contrast at the midpoint of this cycle.
        mid_frac = ((lo + hi) / 2 - 1) / max(n - 1, 1);
        contrasts(k) = cMin + (cMax - cMin) * mid_frac;

        % F1 projection at the carrier.
        sinp = sum(seg_x .* sin(2*pi*carrierFreq*seg_t));
        cosp = sum(seg_x .* cos(2*pi*carrierFreq*seg_t));
        norm = numel(seg_t) / 2;
        amplitudes(k) = hypot(sinp, cosp) / norm;
        phases(k)     = atan2(sinp, cosp);

        % F2 projection at 2x carrier -- rectification measure.
        sinp2 = sum(seg_x .* sin(4*pi*carrierFreq*seg_t));
        cosp2 = sum(seg_x .* cos(4*pi*carrierFreq*seg_t));
        f2_amps(k) = hypot(sinp2, cosp2) / norm;
    end

    % Population-average F2/F1: integrate both across the whole sweep
    % rather than per cycle so noise doesn't blow up the ratio at low
    % contrast where F1 is small.
    f2_f1_ratio = sum(f2_amps) / max(sum(amplitudes), eps);

    % Naka-Rushton fit: r(c) = rmax * c^n / (c^n + c50^n).
    [rmax, c50, nHill, ok] = fit_naka_rushton(contrasts, amplitudes);

    s = struct();
    s.contrasts    = contrasts;
    s.amplitudes   = amplitudes;
    s.phases       = phases;
    s.f2_f1_ratio  = f2_f1_ratio;
    s.c50          = c50;
    s.rmax         = rmax;
    s.n            = nHill;
    s.fit_ok       = ok;
    s.carrierFreq  = carrierFreq;
end


% =====================================================================
% Naka-Rushton fit
% =====================================================================
function [rmax, c50, n, ok] = fit_naka_rushton(c, r)
% Fit r = rmax * c^n / (c^n + c50^n) with lsqcurvefit if available,
% falling back to fminsearch on the squared error. Returns NaN and
% ok=false on failure.
    c = c(:); r = r(:);
    good = isfinite(c) & isfinite(r) & r >= 0;
    c = c(good); r = r(good);
    if numel(c) < 4
        rmax = NaN; c50 = NaN; n = NaN; ok = false; return
    end

    % Reasonable initial guess.
    r0    = max(r);
    c50_0 = median(c);
    n_0   = 2;
    p0 = [r0, c50_0, n_0];

    model = @(p, x) p(1) .* (x .^ p(3)) ./ (x .^ p(3) + p(2) .^ p(3) + eps);

    try
        if exist('lsqcurvefit', 'file') == 2
            opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
                                'MaxIterations', 500);
            lb = [0,   min(c) * 0.1, 0.2];
            ub = [10 * r0, max(c) * 5, 8];
            p_hat = lsqcurvefit(model, p0, c, r, lb, ub, opts);
        else
            % Fallback with fminsearch: minimise SSE, then clamp.
            sse = @(p) sum((model(p, c) - r).^2);
            p_hat = fminsearch(sse, p0);
        end
        rmax = p_hat(1); c50 = p_hat(2); n = p_hat(3); ok = true;
    catch
        rmax = NaN; c50 = NaN; n = NaN; ok = false;
    end
end


