function pmasks = primaryfilterAN(pnuc, probthresh,areafilter)


    pmasks = false(size(pnuc));
%     logim2 = zeros(size(pnuc));
%     s = logfilter;
%     h  = fspecial('log',s);
   se = strel('disk',2);
    
    for z= 1:size(pnuc,3)
        im = pnuc(:,:,z);
        im2 = imfill(im > probthresh,'holes');% for probabilities exported
        im3 = bwareafilt(im2,[areafilter areafilter*10]);
        im4 = imopen(im3,se);
        pmasks(:,:,z) = im4; 
       
    end
    
   
%     t = graythresh(logim2);
%     
    
%     for z = 1:size(pnuc,3)
%         
%         logimsl = adapthisteq(logim2(:,:,z));
%         im1 = logimsl>bthreshfilter*t;
%         im2 = imfill(im1, 'holes');
%         im2b = bwareafilt(im2,[areafilter areafilter*1000]); % an filer out the junk that is smaller than average nucleu area
%         im3 = imerode(im2b, se);
%        % im2b = bwareafilt(im2,[areafilter areafilter*1000]); % an filer out the junk that is smaller than average nucleu area
%         tmp = bwareaopen(im3,areafilter);
%         
%         pmasks(:,:,z) = tmp;
%     end
    

end
