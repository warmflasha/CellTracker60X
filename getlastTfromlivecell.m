function getlastTfromlivecell(direc2,imagedir1,pl,positions,N)
% get the images for the montage from the last time points of the live cell
% imaging before fixing

%imagedir1 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/nuc_raw');%rawimages_nuc%
timegroup = [];%1 []
ff = readAndorDirectory(imagedir1);
chan = [1];
chanal = 1;
%pl =2;
% positions vector needs to start from zero
for kk=0:(size(positions,2)-1)
    pos = kk;
    [imgsnuc_reader]   =  getrawimgfiles(imagedir1,pl, (pos),timegroup,chanal(1));
    nT = imgsnuc_reader{1}.getSizeT;
    
    for jj =(nT-N)    %loop over time points within a given time group
        
        for m = 1:pl %
            planenuc = imgsnuc_reader{m}.getIndex(0,0, jj - 1) + 1;
            inuc(:,:,m) = bfGetPlane(imgsnuc_reader{m},planenuc);
            
            %icyto(:,:,m) = bfGetPlane(imgscyto_reader{m},planecyto);
        end
        % make the max projection
        for ii=1:pl
            img_now = inuc(:,:,ii);
            if ii==1
                max_img=img_now;
            else
                max_img=max(img_now,max_img);
            end
        end
    end
    if pos <10
        imwrite(max_img,[direc2 'lastTmaxprojection_f000' num2str(pos) '_w000' num2str(ff.w(chan)) '.tif'],'writemode','append','Compression','none');%
    end
    if pos>=10
        imwrite(max_img,[direc2 'lastTmaxprojection_f00' num2str(pos) '_w000' num2str(ff.w(chan)) '.tif'],'writemode','append','Compression','none');%
        
    end
end

end

 
 
 
 