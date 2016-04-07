function [img_nuc_reader] = getrawimgfiles(imagedir,nzslices, pos,timegroup,chan)
%[pnuc, inuc] = readmaskfiles1(maskno, segfiledir, rawfiledir, dirinfo, dirinfo1, nzslices, imageno);
% if all the z slices are separate files


% reading raw data
ff=readAndorDirectory(imagedir);
nz = nzslices;
filename = cell(1,nz);

imgs = cell(1,nz);
img_nuc_reader = cell(1,nz);
%dirinfo1(start1).name
for j=1:nz
    filename{j} = getAndorFileName(ff,pos,ff.t(timegroup),ff.z(j),ff.w(chan));
end
for m = 1:nz
img_nuc_reader{m} = bfGetReader(filename{m});
end

% plane1 = img_nuc_reader.getIndex(0,0, k - 1) + 1;
% nuc_img = bfGetPlane(img_nuc_reader,plane1);%reader.getIndex(z, c, t)

% for m = 1:nz
%     imgs{m} = bfopen(filename{m});  %
%  
% end


end
