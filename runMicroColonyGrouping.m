function [colonies] = runMicroColonyGrouping(dir,fr_stim,pos)
% fr_stim frame at which stimulation was added(bmp4, etc.) or frame at
% which the shift happened
% dir = directory with matfiles, resulting from segementation of ilastik
% masks
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
% TO DO 1: here if there was a shift need to shift the coordinates of nuclei before
% the tracker is run on peaks; also introduce as a parameter, where the bmp
% was added, so that to calculate stats for bf and after
% peaks = cat(1,peaks(j))

end
% TO DO 2: cat in less stupid way
for xx=1%:size(peaksall,2)
peaks = cat(2,peaksall{xx},peaksall{xx+1},peaksall{xx+2},peaksall{xx+3});
imgfiles = cat(2,imgfiles_all{xx},imgfiles_all{xx+1},imgfiles_all{xx+2},imgfiles_all{xx+3});
imgfilescyto = cat(2,imgfilescyto_all{xx},imgfilescyto_all{xx+1},imgfilescyto_all{xx+2},imgfilescyto_all{xx+3});
end

for k=1:length(peaks)             % to ensure that the tracker gets input all nonempty peaks (otherwise if any of the peaks is empty the tracker does not run
        if isempty(peaks{k})
          peaks{k} = zeros(3,7);
          peaks{fr_stim} = zeros(3,7); % zero the point where there was shift
    end
end
save(['Outfile_000' num2str(pos) '_tps'],'peaks','dims','imgfiles','imgfilescyto');% saves the new matfile, containing peaks data for all the timepoints
 matfile = ['Outfile_000' num2str(pos) '_tps'];                      % TO DO: save all the data in this file, not only peaks, cat the imgfiles and imgfilescyto
load(matfile);
 runTrackerEDS(matfile,'newTrackParam');
 load(matfile);

[colonies,peaks]=peaksToMicroColoniesANadjusted(peaks);% for each time frame % here the colonies is a cell array : each cell is a colony object

 save(['Outfile_000' num2str(pos) '_tps'],'peaks','dims','imgfiles','imgfilescyto','colonies','cells');
% save('Outfile_0000_tps','peaks','dims','imgfiles','imgfilescyto','colonies','cells');



end
