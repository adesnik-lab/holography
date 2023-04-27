% function matchSIfov()

figure('Name','ScanImage FOV Match', 'NumberTitle','off');
clf

s = scatter(SImatchXY(1,:),SImatchXY(2,:),[],SImatchProb,'filled');
p = plot(SIxboundary,SIyboundary);
r = rectangle('position',...
    [SImatchRangeXforFine(1) 
     SImatchRangeYforFine(1) 
     SImatchRangeXforFine(2)-SImatchRangeXforFine(1) 
     SImatchRangeYforFine(2)-SImatchRangeYforFine(1)]);
r.EdgeColor='g';

