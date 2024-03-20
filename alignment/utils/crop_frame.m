function frame = crop_frame(im, x, y, amt)

xr = x-amt:x+amt;
yr = y-amt:y+amt;
frame = im(xr, yr);