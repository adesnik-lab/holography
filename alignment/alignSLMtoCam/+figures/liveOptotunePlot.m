function liveOptotunePlot(dataUZ)

newfig('Optotune Live Images')
clf
for i=1:size(dataUZ, 3)
    subplot(6,6,i); 
    imagesc(dataUZ(:,:,i));
    colorbar;
end
pause(0.2)