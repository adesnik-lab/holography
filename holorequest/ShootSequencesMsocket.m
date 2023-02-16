function order = ShootSequencesMsocket(Setup,sequences,masterSocket);
%updated 1/19/21 to inlcude output;

flushMSocket(masterSocket);

sendVar = 'C';
mssend(masterSocket,sendVar);


order = [];
disp('waiting for socket to send sequence number')
while isempty(order)
    order = msrecv(masterSocket,.5);
end
disp(['received sequence of length ' num2str(length(order))]);


if any(order>size(sequences,3))
    disp('ERROR: Sequence error. blanking SLM...')
    blank = zeros(size(sequences,1),size(sequences,2));
    outcome = Function_Feed_SLM(Setup.SLM, blank);
    return
end

T=zeros([1 10E5]);
T2=zeros([1 10E5]);

O = zeros([1 10E5]);
% t=tic;

timeout = false;
counter = 1;

useSmallOrder =1;

if useSmallOrder
    [itemsUsed, ~, smallOrder] = unique(order);
    smallSeq = sequences(:,:,itemsUsed);
    
    order = smallOrder;
    sequences = smallSeq;
end



saveDetails =1;
if saveDetails
    
    SLM = Setup.SLM;
    while ~timeout && counter<=length(order)
        %disp(['now queuing hologram ' num2str(order(counter))])
        t=tic;
        %     outcome = Function_Feed_SLM(Setup.SLM, sequences(:,:,order(counter)));
%         calllib('Blink_C_wrapper', 'Write_image', 1, sequences(:,:,order(counter)), 1920*1152, SLM.wait_For_Trigger, SLM.external_Pulse, SLM.timeout_ms);
        if SLM.is_onek
            calllib('Blink_C_wrapper', 'Write_image', 1, sequences(:,:,order(counter)), SLM.Nx*SLM.Ny, SLM.wait_For_Trigger,0, 1, 0, SLM.timeout_ms);
        else
            calllib('Blink_C_wrapper', 'Write_image', 1, sequences(:,:,order(counter)), 1920*1152, SLM.wait_For_Trigger, SLM.external_Pulse, SLM.timeout_ms);
        end
        T(counter)=toc(t);
        
        t = tic;
        outcome = calllib('Blink_C_wrapper', 'ImageWriteComplete', 1, SLM.timeout_ms);
        
        T2(counter)=toc(t);
        O(counter) = outcome;
        if outcome == -1
            timeout = true;
        end
        counter = counter+1;
    end
    
else
    while ~timeout && counter<=length(order)
        %disp(['now queuing hologram ' num2str(order(counter))])
        outcome = Function_Feed_SLM(Setup.SLM, sequences(:,:,order(counter)));
        if outcome == -1
            timeout = true;
        end
        counter = counter+1;
    end
end

if ~timeout
    disp('completed sequence to the end')
else
    disp(['timeout while waiting to display hologram order ' num2str(counter-1)]);
end

%t;