function frame = mean_img(im, bgd)

frame = mean(im, 3);

if nargin == 2
    frame = max(frame-bgd, 0);
end
