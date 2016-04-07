% AW
function showMaskWithNumber(mask)

bwmask = mask > 0;
stats = regionprops(bwmask,mask,'MeanIntensity','Centroid');
xy = stats2xy(stats);
m_int = [stats.MeanIntensity];

imshow(bwmask); hold on;
for ii=1:length(m_int)
    text(xy(ii,1),xy(ii,2),int2str(m_int(ii)),'Color','m');
end