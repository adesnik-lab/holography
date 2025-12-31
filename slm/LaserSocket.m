classdef LaserSocket < handle
    properties
        port
        sock
        srvsock
    end

    methods
        function obj = LaserSocket(port)
            fprintf('Waiting for msocket communication From DAQ... ')
            obj.port = port;
            obj.srvsock = mslisten(port);
            obj.sock = msaccept(obj.srvsock, 15);
            msclose(obj.srvsock);
            mssend(obj.sock, 'A');
            
            % complete handshake
            invar = [];
            while ~strcmp(invar, 'B')
                invar = msrecv(obj.sock, .5);
            end
            fprintf('done.\r')
        end

        function set_power(obj, power)
            % -- array is [power, DE, nTargets]
            % note: we typically do not adjust the power by changing the 
            % number of targets (at least in calibration), we just 
            % increase the total request power in the main calibration code
            % (eg. multi_target_power = pwr * ntargs). though technically, 
            % the code in alignCodeDAQ.m could compute it as 
            % (power*nTargets)/DE, but we aren't using that here, so the
            % last value is set to 1, DE is also set to 1 because we want
            % all holograms to be the same power (again, at in 
            mssend(obj.sock, [power/1000 1 1]);

            % complete handshake
            invar=[];
            while ~strcmp(invar, 'gotit')
                invar = msrecv(obj.sock,0.01);
            end
        end

        function flush(obj)
            invar = 'flush';
            while ~isempty(invar)
                invar = msrecv(obj.sock,0.01);
            end
        end
    end
end

