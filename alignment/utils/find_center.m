function [ x,y ] = find_center( h )
h = imgaussfilt(h, 2);
[~,y] = max(max(h));
[~,x] = max(h(:,y));

%[LX,LY] = size(h);
%UX = 1:LX;
%UY = 1:LY;
%[XX,YY] = ndgrid(UX,UY);
%mass = sum(sum(double(h)));
%x = sum(sum(XX.*double(h)))/mass;
%y = sum(sum(YY.*double(h)))/mass;
end