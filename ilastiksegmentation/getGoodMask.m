% AW, get the labeled mask in 3D
function newmask = getGoodMask(maskz,nucleilist)

ncell = size(nucleilist,1);
nimage = size(nucleilist,2);
newmask = zeros(size(maskz));

for ii = 1:nimage
    newplane = zeros(size(newmask,1),size(newmask,2));
    for jj = 1:ncell
        if ~isnan(nucleilist(jj,ii))
            newplane(maskz(:,:,ii)==nucleilist(jj,ii)) = jj;
        end
    end
    newmask(:,:,ii) = newplane;
end
