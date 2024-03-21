classdef SIsocket < handle
    properties
        port
        sock
        srvsock
    end

    methods
        function obj = SIsocket(port)
            fprintf('Waiting for msocket communication with ScanImage... ')
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

        function send(obj, msg)
            mssend(obj.sock, msg);
        end
        % does this need flush functionality?
    end

end