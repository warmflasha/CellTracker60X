
%function [maskz] = 3Dlblmask(CC,nucleilist)

%%

[PILsn,PILsSourcen, masterCCn, stats, nucleilist, zrange,CC] = traceobjectsz(smasks, userParam.matchdistance, zrange, userParam.zmatch);
%%

size(CC,2); % number of zplanes
size(nucleilist,2); % oved how many plane the niclei are spread
size(nucleilist,1); % howmany objects were found in plane 1
badind = cellfun(@isempty,CC);
CC(badind) = [];

 %%
   % from AW
   % get the masks based on nuclei list
  for ii = 1:length(CC)
      newplane = zeros(1024);
      for jj = 1:length(CC{ii}.PixelIdxList)
            newplane(CC{ii}.PixelIdxList{jj}) = jj;
      end
      trymask(:,:,ii) = newplane;
  end
  %  plot these
for k =1:4
figure(1),subplot(2,2,k); showMaskWithNumber(trymask(:,:,k));
end
  newmask = getGoodMask(trymask,nucleilist);
for k =1:4
figure(2),subplot(2,2,k); showMaskWithNumber(newmask(:,:,k));
end
  
 %%
 % get the mean intensitied from regionprops
 % newmask = the 3dlabel mask
 stats3d = regionprops(newmask,inuc(:,:,1:4),'MeanIntensity','Centroid');
 xyz_fin = cat(1,stats3d.Centroid);
% stats = regionprops(Lnuc,'Centroid','PixelIdxList');   % from 2D code
% xy = [stats.Centroid];
% xx = xy(1:2:end);
% yy=xy(2:2:end);

%vImg = mkVoronoiImageFromPts([xx' yy'],[1024 1024]);

%[V,C] = voronoin(X)

[V,C] = voronoin(xyz_fin);



%cc = bwconncomp(LcytoIl+Lnuc & ~vImg);    % from 2D code
 
 
 
 
 
 %