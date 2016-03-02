function [pnuc] = readmaskfiles1(ilastikNucAll,tpt)
%[pnuc, inuc] = readmaskfiles1(maskno, segfiledir, rawfiledir, dirinfo, dirinfo1, nzslices, imageno);
% if all the z slices are separate files

% reading masks

 nzslices = size(ilastikNucAll,2);
for m=1:nzslices
ilastikfile=ilastikNucAll{m};
io = h5read(ilastikfile,'/exported_data');
io = io(2,:,:,:);
io1 = squeeze(io);
io1_t = io1(:,:,tpt);
pnuc(:,:,m) = io1_t'; % flipped here
end

end

