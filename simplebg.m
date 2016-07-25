function [I2_bgsubtract] = simplebg(Lcytofin,Lnuc,I2)

allmask = Lcytofin | Lnuc;
%figure; imshow(allmask)
imgNoCells = I2;
imgNoCells(allmask) = 0;
imgNoCells = imfill(imgNoCells);

%imshow(imgNoCells)
%imshow(imgNoCells,[])

I2_bgsubtract = imsubtract(I2,imgNoCells);

%imshow(I2_bgsubtract,[])
end