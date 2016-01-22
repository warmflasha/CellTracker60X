function  MaxProj_alldataseparateAN(direc,direc2,tg,chan)
% direc = directory with raw images

% tg = number fo time points to be used from the dataset
% chan = corresponds to ff.w(chan), e.g. ff.w(1) = nuclear channel
% need to run this function separately for each channel
% direc2 = where to save the output projections
%

ff=readAndorDirectory(direc);
if isempty(tg)
tg = size(ff.t,2);
end
Npos = length(ff.p);

for j= 1:(Npos)
for xx = 1:(tg)
max_img =andorMaxIntensity(ff,j,xx,ff.w(chan));

   if j < 10
    imwrite(max_img,[direc2 'Maxprojection_f000' num2str(j) '_t000' num2str(1) '_w000' num2str(ff.w(chan)) '.tif'],'writemode','append','Compression','none');%
   end
   if j >= 10
    imwrite(max_img,[direc2 'Maxprojection_f00' num2str(j) '_t000' num2str(1) '_w000' num2str(ff.w(chan)) '.tif'],'writemode','append','Compression','none');%
   end
end 
end
end
