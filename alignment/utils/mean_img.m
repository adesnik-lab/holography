function frame = mean_img(im, bgd)

frame = mean(im, 3);

if nargin == 2
    frame = max(data-bgd, 0);
end
