function flushMSocket(masterSocket)

invar='flush';
while ~isempty(invar)
    invar = msrecv(masterSocket,0.01);
end