function ff = readFISHdir(direc, diffsamples)

%%

%direc: nuc chnnel image directory
%diffsamples: number of different samples/conditions. Each separate
%chip/dish refers to a separate sample.

%%

for i = 1:diffsamples
    ff.samples(i) = i;
    
    
    prefix = strcat('fish', int2str(i));
    
    clear filesnew;
    
    filesnew = dir([direc filesep strcat(prefix, '_*')]);
    nImagesnew = length(filesnew);
    
    lastImagenew = filesnew(nImagesnew).name;
    
    f_index = regexp(lastImagenew, '_f');
    f_indexn = str2num(lastImagenew(f_index+2:f_index+5));
    
    ff.positions(i) = f_indexn+1;
    
    for j = 1:f_indexn+1
        posname  = sprintf('_f%04d', j-1);
        posnew = dir([direc filesep prefix strcat(posname, '*')]);
        nimages = length(posnew);
        pos_lastimage = posnew(nimages).name;
        
        z_index = regexp(pos_lastimage, '_z');
        z_indexn = str2num(pos_lastimage(z_index+2: z_index+5));
        
        z_indices(j) = z_indexn+1;
        
        
    end
    
    ff.zslices{i} = z_indices;
    
end


end