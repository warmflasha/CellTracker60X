function [pnuc, inuc] = readmaskfilesAN(ilastikNucAll,imagedir, pos,tpt, timegroup,chan)%(maskno, segfile, rawfile, dirinfo, dirinfo1,  nzslices, pos)
%[pnuc, inuc] = readmaskfiles1(maskno, segfiledir, rawfiledir, dirinfo, dirinfo1, nzslices, imageno);
% if all the z slices are separate files

% reading masks

 nzslices = size(ilastikNucAll,2);
for m=1:nzslices
ilastikfile=ilastikNucAll(m).name;
io = h5read(ilastikfile,'/exported_data');
io = io(2,:,:,:);
io1 = squeeze(io);
io1_t = io1(:,:,tpt);
pnuc(:,:,m) = io1_t'; % flipped here
end

% reading raw data
ff=readAndorDirectory(imagedir);
nz = nzslices;
filename = cell(1,nz);

imgs = cell(1,nz);
%dirinfo1(start1).name
for j=1:nz
    filename{j} = getAndorFileName(ff,pos,ff.t(timegroup),ff.z(j),ff.w(chan));
    
end


for m = 1:nz
    imgs{m} = bfopen(filename{m});  %
    img_now = imgs{m}{1}{tpt,1};    % open each timegroup and each z slice
    inuc(:,:,m) = img_now;
    
end


end

