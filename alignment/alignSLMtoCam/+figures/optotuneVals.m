function optotuneVals(fits)

camXYZ = fits.camXYZ;
camPower = fits.camPower;

newfig('Optotune Fluorescence')
f = gcf();
f.Units = 'normalized';
f.Position = [0 0.4 0.4 0.4];

scatter3(camXYZ(1,:), camXYZ(2,:), camXYZ(3,:), [], camPower, 'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Measured Fluorescence intensity by space')
c = colorbar;
c.Label.String = 'Fluorescent Intensity';
axis square