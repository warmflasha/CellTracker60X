function [colonies] = runMicroColonyGrouping(dir,fr_stim,pos,isshift)

% HERE the peaks from all different time groups are merged into one big
% outfile, then this peaks go throus tracker and then through colony
% analysis, output peaks have 9 columns
% fr_stim =  frame at which stimulation was added(bmp4, etc.) or frame at
% which the shift happened (still need to incorporate the shift correction)
% dir = directory with matfiles, resulting from segementation of ilastik
% masks
% pos = same, the frame number that is being procesed (starts from 0)
% 
% TO DO:    need to add the condition that if dt = 1, then no need to
% concatenate anything and the peaks is already containing all frames, so
% need to go directly to running the tracker of them (end of this function)

[nums, files]=folderFilesFromKeyword(dir,['Outfile_000' num2str(pos) '_t']);%'Outfile_0000_t' 

peaksall = cell(1,length(nums));
imgfiles_all = cell(1,length(nums));
imgfilescyto_all= cell(1,length(nums));
for j=1:length(nums)

matfile = files(j).name;
load(matfile,'peaks','peaks','dims','imgfiles','imgfilescyto');

peaksall{j} = peaks;
imgfiles_all{j} = imgfiles;
imgfilescyto_all{j} = imgfilescyto;

end
% TO DO 2: cat in less stupid way
for xx=1%:size(peaksall,2)
peaks = cat(2,peaksall{xx},peaksall{xx+1},peaksall{xx+2},peaksall{xx+3});%,
imgfiles = cat(2,imgfiles_all{xx},imgfiles_all{xx+1},imgfiles_all{xx+2},imgfiles_all{xx+3});
imgfilescyto = cat(2,imgfilescyto_all{xx},imgfilescyto_all{xx+1},imgfilescyto_all{xx+2},imgfilescyto_all{xx+3});
end
% introducing the shift, if it occured during addition of bmp4
if isshift == 1
for k=fr_stim+1:length(peaks)             % to ensure that the tracker gets input all nonempty peaks (otherwise if any of the peaks is empty the tracker does not run
    if ~isempty(peaks{k})
        n1 = uncompressBinaryImg(imgfiles(fr_stim-1).compressNucMask); % this chunck is to calculate the actual shift vector, and need not be performed in the loop
        n2 = uncompressBinaryImg(imgfiles(fr_stim+1).compressNucMask);
        bw1 = bwconncomp(n1);
        bw2 = bwconncomp(n2);
        stats1 = regionprops(bw1,'Centroid');
        stats2 = regionprops(bw2,'Centroid');
        xy1 = [stats1(1).Centroid];
        xy2 = [stats2(1).Centroid];
        diffx = sqrt(power(xy1(1)-xy2(1),2));
        diffy = sqrt(power(xy1(2)-xy2(2),2));
        shift = [diffx,diffy];
        
        %            %here need to introduce the shift
        peaks{k}(:,1) = round(peaks{k}(:,1)+shift(1)); % zero the point where there was shift
        peaks{k}(:,2) = round(peaks{k}(:,2)+shift(2));
        
    end
    peaks{fr_stim} = [];
end
end

peaks = peaks(~cellfun(@isempty, peaks)); % to ensure that the tracker gets input all nonempty peaks (otherwise if any of the peaks is empty the tracker does not run

save(['Outfile_000' num2str(pos) '_tps'],'peaks','dims','imgfiles','imgfilescyto');% saves the new matfile, containing peaks data for all the timepoints
 matfile = ['Outfile_000' num2str(pos) '_tps'];                      % TO DO: save all the data in this file, not only peaks, cat the imgfiles and imgfilescyto
load(matfile);
 runTrackerEDS(matfile,'newTrackParam');
 load(matfile);

[colonies,peaks]=peaksToMicroColoniesANadjusted(peaks);% for each time frame % here the colonies is a cell array : each cell is a colony object

 save(['Outfile_000' num2str(pos) '_tps'],'peaks','dims','imgfiles','imgfilescyto','colonies','cells');
% save('Outfile_0000_tps','peaks','dims','imgfiles','imgfilescyto','colonies','cells');



end
