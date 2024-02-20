classdef RESTio < handle
    properties
        server
        port
        ops
        debug = false
    end

    methods
        function obj = RESTio(server, port)
            obj.server = server;
            obj.port = port;
            obj.ops = weboptions('MediaType', 'application/json');
        end

        function send(obj, data, target, sender)
                recv = webwrite(obj.get_url('msg', target), obj.encode(data, sender), obj.ops);
            if obj.debug
                disp(recv)
            end
        end

        function recv = scan(obj, target)
            try
                recv = webread(obj.get_url('msg', target), obj.ops); % cause error 404 when nothing
            catch ME
                if ~strcmp(ME.identifier, 'MATLAB:webservices:HTTP404StatusCodeError')
                    warning(ME.message)
                end

                recv = [];

                return
            end

            if strcmp(recv.message_status, 'read')
                recv = [];
            end
        end

        function [out, recv] = read(obj, target, timeout)
            tic
            recv = [];
            while isempty(recv) && toc < timeout
                recv = scan(obj, target);
            end

            out = obj.decode(recv);

            if obj.debug
                disp(recv)
            end
        end

        function config(obj, data, target)
            webwrite(obj.get_url('config', target), obj.encode(data), obj.ops);
        end

        function flush(obj, target)
            try
                % delete everything on the server?? (bad?)
                webread(obj.get_url('msg', target), weboptions('RequestMethod', 'delete'));
            catch ME
                warning(ME.message)
            end
        end

        function out = encode(obj, data, sender)
            out = struct();
            out.sender = sender;
            out.message = mps.json.encode(data, 'Format', 'large');
            out = mps.json.encode(out, 'Format', 'small'); % explicitly short
        end

        function out = decode(obj, data)
            if isempty(data)
                out = [];
                return
            end
            out = mps.json.decode(data.message);
        end

        function reset(obj, target)
            try
                webread(obj.get_url('db', target), weboptions('RequestMethod', 'delete'));
            end
        end

        function url = get_url(obj, path, target)
            if nargin < 3 || isempty(target)
                target = [];
            end
            url = sprintf('%s/%s/%s', obj.server, path, target) ;
        end

    end
end