function newfig(name)

h = findobj(0,'Name',name);
if isempty(h)
    h = figure('Name', name, 'NumberTitle','off');
end
figure(h);