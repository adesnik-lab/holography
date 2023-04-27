classdef bascam < handle

    properties
        vid
        src
        camera_model = 'Y800_2048x1536';
        driver = 'winvideo';
        is_running = 0;
        preview_bit_depth = 'auto';
    end

    properties (Constant)
        bit_depth = 8;
        castFun = @uint8;
        castAs = 'uint8';
        camMax = 2^8-1;
    end

    properties (Dependent)
        exposure
        resolution
    end

    methods

        function obj = bascam(model, driver)
            if nargin == 2
                obj.driver = driver;
                obj.camera_model = model;
            elseif nargin == 1
                obj.camera_model = model;
            end
            fprintf('Initializing Basler... ')
            obj.vid = videoinput(obj.driver, 1, obj.camera_model);
            obj.vid.ReturnedColorspace = 'grayscale';
            obj.vid.FramesPerTrigger = 1;
            obj.vid.TriggerRepeat = Inf;
            triggerconfig(obj.vid, 'manual');
            obj.src = getselectedsource(obj.vid);
            obj.src.Exposure = -3;
            obj.src.Brightness = 33;
            fprintf('done.\r')
        end

        function start(obj)
            if ~obj.is_running
                start(obj.vid);
                obj.is_running = 1;
            end
        end

        function stop(obj)
            if obj.is_running
                stop(obj.vid);
                obj.is_running = 0;
            end
        end

        function preview(obj)
            f = figure('Name','Basler Preview', 'NumberTitle','off');
            f.Position = [817 712 1000 800];
            movegui(f, 'center')
            while isgraphics(f)
                frame = obj.grab(1);
                imagesc(frame)
                axis image
                if isa(obj.preview_bit_depth, 'double')
                    caxis([0 obj.preview_bit_depth])
                end
                colorbar
                drawnow
            end
        end

        function frames = grab(obj, nframes)
            if obj.vid.FramesPerTrigger ~= nframes
                obj.stop()
                obj.vid.FramesPerTrigger = nframes;
                obj.start()
            end

            trigger(obj.vid);
            x = getdata(obj.vid);

            if nframes == 1
                frames = x;
            else
                frames = squeeze(x(:,:,1,:));
            end
        end

        function set.exposure(obj, exposure)
            obj.stop()
            obj.src.exposure = exposure;
            obj.start()
        end

        function exposure = get.exposure(obj)
            exposure = obj.src.Exposure;
        end

        function res = get.resolution(obj)
            res = obj.vid.VideoResolution;
        end

        function reset(obj)
            fprintf('Restarting camera... ')
            obj.stop()
            pause(1)
            obj.start()
            fprintf('OK.\r')
        end

        function delete(obj)
            delete(obj.vid)
        end

    end

end