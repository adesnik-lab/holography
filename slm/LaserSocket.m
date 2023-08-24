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
            mssend(obj.sock, [power/1000 1 1]);

            % complete handshake
            invar=[];
            while ~strcmp(invar, 'gotit')
                invar = msrecv(obj.sock,0.01);
            end
        end
    end
end

