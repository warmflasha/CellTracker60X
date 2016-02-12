function [pnuc, inuc] = readmaskfiles1(maskno, segfile, rawfile, dirinfo, dirinfo1,  nzslices, pos)

    % reading masks
    
    m = 1;
    for i = maskno:maskno+nzslices-1
        filename = strcat(segfile,'/', dirinfo(i).name);
        file = h5read(filename, '/exported_data');
%         pnuc(:,:,m) = squeeze(file(1,:,:));
%         pnuc(:,:,m) = pnuc(:,:,m)';
%         m = m+1;
        
    pnuc(:,:,m) = file(2,:,:,pos);% for probabilities exported
    pnuc(:,:,m) = squeeze(pnuc(:,:,pos));
    pnuc(:,:,m) = pnuc(:,:,pos)';
    m = m+1;
    end
    
    
    % reading raw data
    ff=readAndorDirectory(direc);
timegroups = size(ff.t,2);
nz = size(ff.z,2);
filename = cell(1,nz);
imgs = cell(1,nz);
for xx = 1:size(tg,2)
for j=1:nz
    
    filename{j} = getAndorFileName(ff,pos,ff.t(tg(xx)),ff.z(j),ff.w(chan));
end

% filename{1} = plane z0000; now if open it with bfopen will uncover all
% the timepoints within first time group taken at z0000

for x = 1:nz
    imgs{x} = bfopen(filename{x});% open each timegroup and each z slice
    
end
    
    
%     m1 = 1;
%     for j = imageno:imageno+nzslices-1
%         filename = strcat(rawfile, '/', dirinfo1(j).name);
%         inuc(:,:,m1) = imread(filename);
%         m1 = m1+1;
%     end
end