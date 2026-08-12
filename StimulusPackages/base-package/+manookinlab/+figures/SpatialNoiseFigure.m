classdef SpatialNoiseFigure < symphonyui.core.FigureHandler
    properties (SetAccess = private)
        device
        recordingType
        numXChecks
        numYChecks
        noiseClass
        chromaticClass
        preTime
        stimTime
        frameRate
        numFrames
        stixelSize
        stepsPerStixel
        frameDwell
        gaussianFilter
        filterSdStixels
    end

    properties (Access = private)
        axesHandle
        imgHandle
        strf
        spaceFilter
        nChannels
        xaxis
        yaxis
    end

    methods

        function obj = SpatialNoiseFigure(device, varargin)
            ip = inputParser();
            ip.addParameter('recordingType', 'extracellular', @(x)ischar(x));
            ip.addParameter('stixelSize', [], @(x)isfloat(x));
            ip.addParameter('numXChecks', [], @(x)isfloat(x));
            ip.addParameter('numYChecks', [], @(x)isfloat(x));
            ip.addParameter('noiseClass', 'binary', @(x)ischar(x));
            ip.addParameter('chromaticClass', 'achromatic', @(x)ischar(x));
            ip.addParameter('preTime', 0.0, @(x)isfloat(x));
            ip.addParameter('stimTime', 0.0, @(x)isfloat(x));
            ip.addParameter('frameRate', 60.0, @(x)isfloat(x));
            ip.addParameter('numFrames', [], @(x)isfloat(x));
            ip.addParameter('stepsPerStixel', 1, @(x)isfloat(x));
            ip.addParameter('frameDwell', 1, @(x)isfloat(x));
            ip.addParameter('gaussianFilter', false, @(x)islogical(x) || isfloat(x));
            ip.addParameter('filterSdStixels', 1.0, @(x)isfloat(x));

            ip.parse(varargin{:});

            obj.device = device;
            obj.recordingType = ip.Results.recordingType;
            obj.stixelSize = ip.Results.stixelSize;
            obj.numXChecks = ip.Results.numXChecks;
            obj.numYChecks = ip.Results.numYChecks;
            obj.noiseClass = ip.Results.noiseClass;
            obj.chromaticClass = ip.Results.chromaticClass;
            obj.preTime = ip.Results.preTime;
            obj.stimTime = ip.Results.stimTime;
            obj.frameRate = ip.Results.frameRate;
            obj.numFrames = ip.Results.numFrames;
            obj.stepsPerStixel = ip.Results.stepsPerStixel;
            obj.frameDwell = ip.Results.frameDwell;
            obj.gaussianFilter = logical(ip.Results.gaussianFilter);
            obj.filterSdStixels = ip.Results.filterSdStixels;

            % Number of analyzed color channels.
            if strcmpi(obj.chromaticClass, 'BY')
                obj.nChannels = 2; % yellow, blue
            elseif strcmpi(obj.chromaticClass, 'RGB')
                obj.nChannels = 3;
            else
                obj.nChannels = 1;
            end

            % Set the x/y axes
            obj.xaxis = linspace(-obj.numXChecks/2,obj.numXChecks/2,obj.numXChecks)*obj.stixelSize;
            obj.yaxis = linspace(-obj.numYChecks/2,obj.numYChecks/2,obj.numYChecks)*obj.stixelSize;

            obj.createUi();
        end

        function createUi(obj)
            import appbox.*;

            for k = 1 : 4
            obj.axesHandle(k) = subplot(2, 2, k, ...
                'Parent', obj.figureHandle, ...
                'FontUnits', get(obj.figureHandle, 'DefaultUicontrolFontUnits'), ...
                'FontName', get(obj.figureHandle, 'DefaultUicontrolFontName'), ...
                'FontSize', get(obj.figureHandle, 'DefaultUicontrolFontSize'), ...
                'XTickMode', 'auto');
            end

            obj.strf = zeros(obj.numYChecks, obj.numXChecks, floor(obj.frameRate*0.5), obj.nChannels);
            obj.spaceFilter = [];

            obj.setTitle([obj.device.name ' receptive field']);
        end

        function setTitle(obj, t)
            set(obj.figureHandle, 'Name', t);
        end

        function clear(obj)
            cla(obj.axesHandle);
            obj.strf = zeros(obj.numYChecks, obj.numXChecks, floor(obj.frameRate*0.5), obj.nChannels);
            obj.spaceFilter = [];
            % Set the x/y axes
            obj.xaxis = linspace(-obj.numXChecks/2,obj.numXChecks/2,obj.numXChecks)*obj.stixelSize;
            obj.yaxis = linspace(-obj.numYChecks/2,obj.numYChecks/2,obj.numYChecks)*obj.stixelSize;
        end

        function handleEpoch(obj, epoch)
            if ~epoch.hasResponse(obj.device)
                error(['Epoch does not contain a response for ' obj.device.name]);
            end

            response = epoch.getResponse(obj.device);
            [quantities, ~] = response.getData();
            sampleRate = response.sampleRate.quantityInBaseUnits;
            prePts = round((obj.preTime*1e-3 - 1/60)*sampleRate);

            if numel(quantities) > 0
                % Parse the response by type.
                y = manookinlab.util.responseByType(quantities, obj.recordingType, obj.preTime, sampleRate);

                if strcmp(obj.recordingType,'extracellular') || strcmp(obj.recordingType, 'spikes_CClamp')
                    y = BinSpikeRate(y(prePts+1:end), obj.frameRate, sampleRate);
                else
                    % Bandpass filter to get rid of drift.
                    y = bandPassFilter(y, 0.2, 500, 1/sampleRate);
                    if prePts > 0
                        y = y - median(y(1:prePts));
                    else
                        y = y - median(y);
                    end
                    y = binData(y(prePts+1:end), obj.frameRate, sampleRate);
                end

                % Pull the epoch parameters.
                p = epoch.parameters;
                seed = p('seed');
                numXStixels = p('numXStixels');
                numYStixels = p('numYStixels');
                nFrames = getParameter(p, 'numFrames', obj.numFrames);
                sps = getParameter(p, 'stepsPerStixel', obj.stepsPerStixel);
                dwell = getParameter(p, 'frameDwell', obj.frameDwell);
                uniqueFrames = getParameter(p, 'unique_frames', []);
                repeatFrames = getParameter(p, 'repeat_frames', 0);
                repeatingSeed = getParameter(p, 'repeating_seed', 1);

                % Make the response the same size as the stim frames.
                y = y(1 : min(numel(y), nFrames));

                % Columate.
                y = y(:);

                % Regenerate the stimulus and compute the STRF by reverse
                % correlation. Channels: achromatic -> 1; BY -> [yellow,
                % blue]; RGB -> [red, green, blue].
                strfTmp = manookinlab.util.analyzeSpatialNoise(y, ...
                    'numXStixels', numXStixels, 'numYStixels', numYStixels, ...
                    'numXChecks', obj.numXChecks, 'numYChecks', obj.numYChecks, ...
                    'chromaticClass', obj.chromaticClass, ...
                    'numFrames', nFrames, 'stepsPerStixel', double(sps), ...
                    'seed', seed, 'frameDwell', double(dwell), ...
                    'uniqueFrames', uniqueFrames, 'repeatFrames', repeatFrames, ...
                    'repeatingSeed', repeatingSeed, ...
                    'gaussianFilter', obj.gaussianFilter, ...
                    'filterSdStixels', obj.filterSdStixels, ...
                    'frameRate', obj.frameRate, ...
                    'filterFrames', floor(obj.frameRate*0.5));

                obj.strf = obj.strf + strfTmp;
                lobePts = 2 : 4;
                obj.spaceFilter = reshape(mean(obj.strf(:,:,lobePts,:),3), ...
                    [obj.numYChecks, obj.numXChecks, obj.nChannels]);

                % Display the spatial RF at four time lags.
                if obj.nChannels == 1
                    for k = 1 : 4
                        imagesc('XData',obj.xaxis,'YData',obj.yaxis,...
                            'CData', obj.strf(:,:,k+2), 'Parent', obj.axesHandle(k));
                        axis(obj.axesHandle(k),'image');
                        colormap(obj.axesHandle(k), 'gray');
                        title(obj.axesHandle(k),['t: -',num2str(round((k+1)/obj.frameRate*1000)),' ms']);
                    end
                else
                    % Combined color maps, jointly normalized across
                    % channels so relative amplitudes are preserved.
                    maxAbs = max(abs(obj.strf(:)));
                    if maxAbs == 0
                        maxAbs = 1;
                    end
                    for k = 1 : 4
                        img = zeros(obj.numYChecks, obj.numXChecks, 3);
                        if obj.nChannels == 2 % BY: yellow -> R+G, blue -> B
                            img(:,:,1) = 0.5 + 0.5*obj.strf(:,:,k+2,1)/maxAbs;
                            img(:,:,2) = img(:,:,1);
                            img(:,:,3) = 0.5 + 0.5*obj.strf(:,:,k+2,2)/maxAbs;
                        else % RGB
                            for c = 1 : 3
                                img(:,:,c) = 0.5 + 0.5*obj.strf(:,:,k+2,c)/maxAbs;
                            end
                        end
                        img = min(max(img, 0), 1);
                        imagesc('XData',obj.xaxis,'YData',obj.yaxis,...
                            'CData', img, 'Parent', obj.axesHandle(k));
                        axis(obj.axesHandle(k),'image');
                        title(obj.axesHandle(k),['t: -',num2str(round((k+1)/obj.frameRate*1000)),' ms']);
                    end
                end
            end
        end

    end
end

function v = getParameter(p, key, defaultValue)
% Pull an epoch parameter with a fallback default.
if isKey(p, key)
    v = p(key);
else
    v = defaultValue;
end
end
