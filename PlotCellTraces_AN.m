
function [datcell] = PlotCellTraces_AN(dir,col,N,fr_stim)

%plot cell traces for specific colonies using the output from the runTracker
%EDS and runcolonygrouping
% in that output last column of peaks = has the trajectory number that the
% cell belongs to
% for this peaks should have 9 columns: 9th has the size of the colony that
% the cell belongs to
% 

[nums, files]=folderFilesFromKeyword(dir,'Outfile_');


for j=1:length(nums)

matfile = files(j).name;

load(matfile,'peaks','cells'); % AN

for k=1:length(peaks)
if ~isempty(peaks{k})
trN(k) = max(peaks{k}(:,8));
end
end
trN = max(trN);
trN =(1:trN);

 %trN = max(peaks{1}(:,4)); % number of trajectories within the frame

alldata=zeros(length(peaks),length(trN));%length(trN)
 
for xx = 1:length(trN)
    
    for k=1:length(peaks)
        if ~isempty(peaks{k})
            
            nc = size(peaks{k},1);% number of cells found within each frame k
            for i=1:nc
                
                if peaks{k}(i,9) == N && peaks{k}(i,8) == trN(xx);%
                    alldata(k,xx) = peaks{k}(i,col(1))./peaks{k}(i,col(2));
                    
                end
                
            end
            
        end
    end
end


datcell{j} = alldata;
end

bkgsign = zeros(length(datcell),3);
for j=1:length(datcell)
    for k=1:size(datcell{j},2)
        if ~isempty(datcell{j}(:,k))
            plot(datcell{j}(:,k),'-*');
            hold on
            bkgsign(j,1) = mean(datcell{j}(1:fr_stim,k));% mean signaling before stimulation
            bkgsign(j,2) = mean(datcell{j}(fr_stim:end,k));% mean signaling after stimulation
            bkgsign(j,3) = abs(datcell{j}(fr_stim-1,k)-max(datcell{j}(fr_stim:end,k)));% value of the jump upon stimulation
            %bkgsign(j,3) = abs(bkgsign(j,1)-bkgsign(j,2));
        end
    end
    ylim([0 2.4])
    title(['CellTraces for colonies of size ' num2str(N) ]);
    
    
end

%save(['CellTraces_' num2str(N) ],'datcell');


end






