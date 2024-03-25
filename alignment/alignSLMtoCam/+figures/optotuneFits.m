function optotuneFits(fits, optotunePlanes)

camXYZ = fits.camXYZ;
camToOpto = fits.camToOpto;

newfig('Optotune Fits')
clf
f = gcf();
f.Units = 'normalized';
f.Position = [0 0.4 0.4 0.4];

subplot(2,2,1)
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],optotunePlanes,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Measured Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

subplot(2,2,2)
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],polyvaln(camToOpto,camXYZ(1:3,:)'),'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Estimated Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

subplot(2,2,3)
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],polyvaln(camToOpto,camXYZ(1:3,:)')-optZ,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Error (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

subplot(2,2,4)
c = sqrt((polyvaln(camToOpto,camXYZ(1:3,:)')-optotunePlanes).^2);
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],c,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Error RMS (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

