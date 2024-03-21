function checkHoloPower(frame_crop)

newfig('Hologram Test')
clf
imagesc(frame_crop)
axis square
colorbar
% clim([0 255])
xL=xlim;
yL=ylim;
mx = max(frame_crop,[],"all");
str = ['Peak Intensity: ' num2str(mx) ' A.U.'];
text(0.03*xL(2),0.03*yL(2),str,'HorizontalAlignment','left','VerticalAlignment','top', 'Color','w')