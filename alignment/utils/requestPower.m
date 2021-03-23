function requestPower(pwr,socket)
% Send the DAQ a new laser power so it shoots or kills a
% holo. Waits for confirmation from the master/DAQ computer.
%
% requestPower(pwr)
%   pwr (float): laser power in mW to request
flushSocket(socket)
mssend(socket,[pwr/1000 1 1]);

% Wait for confirmation
invar=[];
c=1;
while ~strcmp(invar,'gotit') && c<1e6
    invar = msrecv(socket,0.01);
    c=c+1;
end