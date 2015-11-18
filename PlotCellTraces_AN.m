
function [datcell] = PlotCellTraces_AN(dir, col ,N )

%plot cell traces for specific colonies using the output from the runTracker
%EDS
% in that output last column of peaks = has the trajectory number that the
% cell belongs to
% in the colonies cell array( saved in the otfile) the last column is the
% number of colony within that frma that the cells belong to

% 1. read in the files in the directory, find the ones that have the
% 'outfile' srting in them
% 2. setup a loop to precess each file (it's peaks colonies and cells
% 3. group the output based on cell belonging to the same trajectory and
% the a colony of size N ( add the number of cells within the colony to the
% cells structure.
% plot the output for different colony sizes

[nums, files]=folderFilesFromKeyword(dir,'Outfile_');


for j=1:length(nums)

matfile = files(j).name;

load(matfile,'peaks','colonies','cells'); % AN

% check the number of trajectories withn the frame and add the loop over
% same trajectories
% 


 trN = max(peaks{1}(:,4)); % number of trajectories within the frame
 trN =(1:trN);


alldata=zeros(length(peaks),size(peaks{1},1));%length(trN)
 
for xx = 1:length(trN)
    
    %q=1;
    for k=1:length(colonies)
        
        if ~isempty(colonies{k})
            nc = size(colonies{k},2);% number of colonies found within each frame k
            for i=1:nc
                if N == 1
                    if colonies{k}(i).ncells == N && colonies{k}(i).data(:,8) == trN(xx);% 
                        alldata(k,xx) = colonies{k}(i).data(:,col(1))./colonies{k}(i).data(:,col(2));
                        
                    end
                end
                if N > 1
                    if colonies{k}(i).ncells == N && colonies{k}(i).data(i,8) == trN(xx);
                        alldata(k,xx) = colonies{k}(i).data(i,col(1))./colonies{k}(i).data(i,col(2));
                        
                    end
                end
                
            end
            %q=q+nlines(k);
        end
    end
    
end

datcell{j} = alldata;
end

%save(['CellTraces_' num2str(N) ],'datcell');


end






