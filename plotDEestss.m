[X,Y] = meshgrid(0.20:0.01:0.8, 0.20:0.01:0.8);
test_pts = [X(:) Y(:)];
test_pts(:,3)=0;
de_list = DEfromSLMCoords(test_pts);
figure(1);
clf
scatter(test_pts(:,1), test_pts(:,2), [], de_list)
colorbar
