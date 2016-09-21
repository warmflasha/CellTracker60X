function [img_nuc_reader] = getrawimgfilesselectZ(imagedir,selectZ, pos,timegroup,chan)
%[pnuc, inuc] = readmaskfiles1(maskno, segfiledir, rawfiledir, dirinfo, dirinfo1, nzslices, imageno);
% if all the z slices are separate files


% reading raw data
ff=readAndorDirectoryANmod(imagedir);
nz = selectZ;
filename = cell(1,1);

imgs = cell(1,1);
img_nuc_reader = cell(1,size(selectZ,2));
%dirinfo1(start1).name
if (isempty(timegroup) == 1) && (isempty(ff.z) == 1)
    for j=1
        filename{j} = getAndorFileName(ff,pos,[],[],ff.w(chan));%%%
    end
else if isempty(timegroup) == 1
        for j=1:size(selectZ,2)
            filename{j} = getAndorFileName(ff,pos,[],ff.z(selectZ(j)),ff.w(chan));%%%
        end
    else if (isempty(ff.z) == 1)
            for j=1
                filename{j} = getAndorFileName(ff,pos,ff.t(timegroup),[],ff.w(chan));%%%
            end
            
        else
            for j=1:size(selectZ,2)
                filename{j} = getAndorFileName(ff,pos,ff.t(timegroup),ff.z(selectZ(j)),ff.w(chan));%%%
            end
            
            
        end
    end
end

    for m = 1:size(selectZ,2)
        img_nuc_reader{m} = bfGetReader(filename{m});
    end

    % plane1 = img_nuc_reader.getIndex(0,0, k - 1) + 1;
    % nuc_img = bfGetPlane(img_nuc_reader,plane1);%reader.getIndex(z, c, t)
    
    % for m = 1:nz
    %     imgs{m} = bfopen(filename{m});  %
    %
    % end
    
    
end
