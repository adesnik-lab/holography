function flushSocket(socket)
% Flush the specified socket of all vars.
%
% flushSocket(socket)
%   socket (var): ref to an int32 address of a socket
%   
%   usage-> flushSocket(masterSocket)

invar='flush';
while ~isempty(invar)
    invar = msrecv(socket,0.01);
end