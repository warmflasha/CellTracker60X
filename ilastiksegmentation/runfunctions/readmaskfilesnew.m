function [pnuc, inuc] = readmaskfilesnew(segfiledir, rawfiledir, imnum, pos,  nzslices, nuc_ch)


% reading raw data and masks

for j = 1:nzslices
    file2read = sprintf('fish%01d_f%04d_z%04d_w%04d.tif', imnum, pos, j-1, nuc_ch);
    filename = strcat(rawfiledir,'/', file2read);
    
    inuc(:,:,j) = imread(filename);
    
    mask2read = sprintf('fish%01d_f%04d_z%04d_w%04d_Probabilities.h5', imnum, pos, j-1, nuc_ch);
    maskname = strcat(segfiledir, '/', mask2read);
    
    mask = h5read(maskname, '/exported_data');
    
    nmask1 = squeeze(mask(1,:,:));
    nmask2 = squeeze(mask(2,:,:));
    
    pnuc(:,:,j) = nmask1;
    pnuc(:,:,j) = pnuc(:,:,j)';
end


end
