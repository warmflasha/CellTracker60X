function [data_perframe] = GetDynamicColonyStats(matfile,fr_stim,delta_t)

load(matfile,'colonies','ncells','peaks');

p = fr_stim*delta_t/60;
colors = colorcube(50);

colgr = size(colonies,2);% how many colonies were found 
ntimes = length(peaks); 

for ii = 1:colgr;
    %if colSZ == ncells{ii}(1); % how many cells were there in the first frame of the i-th colony
    
    Ntr = size(colonies(ii).cells,2); % number of trajectories
    %onframesall = zeros(length(peaks),size(colonies(ii).cells,2));
    data_perframe = zeros(length(peaks),size(colonies(ii).cells,2));
    
    
    for j = 1:Ntr
        one = (colonies(ii).cells(j).fluorData(:,2)./colonies(ii).cells(j).fluorData(:,3));
        tmp = (colonies(ii).cells(j).onframes)';
        data_perframe(tmp(1):tmp(end),j) = one;
        
        for k = 1:ntimes
            data_perframe(k,j) = data_perframe(k,j);
            
            
        end
        
        
    end
    
    data_perframe{ii} = data_perframe;
    %end
end

end