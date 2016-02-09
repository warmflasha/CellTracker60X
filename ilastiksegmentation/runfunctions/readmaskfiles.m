function [pnuc, inuc] = readmaskfiles(maskno, segfile, rawfile, dirinfo, dirinfo1,  nzslices, imageno)

    % reading masks
    m = 1;
    for i = maskno:maskno+nzslices-1
        filename = strcat(segfile,'/', dirinfo(i).name);
        file = h5read(filename, '/exported_data');
        pnuc(:,:,m) = squeeze(file(1,:,:));
        pnuc(:,:,m) = pnuc(:,:,m)';
        m = m+1;
    end
    
    
    % reading raw data
    m1 = 1;
    for j = imageno:imageno+nzslices-1
        filename = strcat(rawfile, '/', dirinfo1(j).name);
        inuc(:,:,m1) = imread(filename);
        m1 = m1+1;
    end
end