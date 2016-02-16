
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
maskz = zeros(1024,1024,size(zrange,2));
Inew = zeros(1024,1024);

for k=1:size(nucleilist,2)
    bw = regionprops(CC{k},'PixelIdxList')
    for xx = 1:size(bw,1)
        
       % if isfinite(nucleilist(xx,k))
            ind = find(~isnan(nucleilist(:,k)));
            Inew(bw(xx).PixelIdxList) =nucleilist(ind(xx),k); % label of the object as it is in the nuclei list, same label for this object in different planes
            
        %end
    end
    maskz(:,:,k) = Inew;
    Inew = zeros(1024,1024);
end
%%

for l=1:size(maskz,3)
    figure, imshow(maskz(:,:,l),[])
end

        