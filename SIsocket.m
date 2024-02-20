classdef SIsocket < handle
    properties
        port
        sock
        srvsock
    end

    methods
        function obj = SIsocket(port)
            fprintf('Waiting for msocket communication From ScanImage... ')
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

        function set_optotune(obj, arr)
            mssend(obj.sock, arr);

            % complete handshake
            invar=[];
            while ~strcmp(invar, 'gotit')
                invar = msrecv(obj.sock,0.01);
            end
        end

        function check_eval(obj)
            fprintf('Checking socket connection with ScanImage...\r')
            mssend(obj.sock,'1+2');
            
            invar=[];
            while ~strcmp(num2str(invar),'3')
                invar = msrecv(obj.sock,0.01);
            end
            fprintf('ScanImage computer linked and is in eval mode for hole burn.\r')
        end

        function cmd(obj, str)
            mssend(obj.socket, str)
        end

        function flush(obj)
            invar = 'flush';
            while ~isempty(invar)
                invar = msrecv(obj.sock,0.01);
            end
        end
    end
end

