function [pnuc] = readmaskfiles1(ilastikNucAll,tpt,lblN)
%[pnuc, inuc] = readmaskfiles1(maskno, segfiledir, rawfiledir, dirinfo, dirinfo1, nzslices, imageno);
% if all the z slices are separate files
% if don't have time point will have less-dim .h5 file

% reading masks

 %nzslices = size(ilastikNucAll,2);
 
for m=1:size(ilastikNucAll,2)
    
ilastikfile=ilastikNucAll{m};
io = h5read(ilastikfile,'/exported_data');
io = io(lblN,:,:,:);% make this into a variable ( which ilastik label to use as signal and which to use as background
io1 = squeeze(io);
io1_t = io1(:,:,tpt);
pnuc(:,:,m) = io1_t'; % flipped here
end

end

