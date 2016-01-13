function [datafin] = GetDynamicColonyStats(matfile,fr_stim,delta_t,flag,colSZ,resptime)
% resptime = how many frames after stimulation want to look at the response
load(matfile,'colonies','ncells','peaks');

p = fr_stim*delta_t/60;
colors = colorcube(50);
ntimes = length(peaks);

colgr = size(colonies,2);% how many colonies were found
datafin = cell(colgr,1); % preallocate , more than necessary col 1 - means before, col2  mean after

for ii = 1:colgr;
    if colSZ == ncells{i}(1); % how many cells were there in the frame before stimulation of the i-th colony,
    
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
    
    datafin{ii}(j,1) = mean(nonzeros(data_perframe(1:fr_stim,j)));% mean before stimulation
    datafin{ii}(j,2) = mean(nonzeros(data_perframe(fr_stim:fr_stim+resptime,j)));% mean after stimulation, 20 frames past stimulation only
    %end
    
    if flag == 1
        p2 = ((resptime)*delta_t)/60;
        figure(10),subplot(1,3,1),plot(nonzeros(datafin{ii}(:,1)),'*','color',colors(ii,:),'markersize',15);
        hold on
        ylim([0 2.5]);
        %xlim([0 lenth(datcell)]);% number traces
        legend('before');
        ylabel('nuc/cyto ratio, mean over time');
        xlabel('Positions');
        title(['microCol of size ' num2str(colSZ) ],'fontsize',15);
        figure(10),subplot(1,3,2),plot(nonzeros(datafin{ii}(:,2)),'*','color',colors(ii,:),'markersize',15);
        hold on
        legend([num2str(p2) 'hours after stimulation']);
        ylim([0 2.5])
        % xlim([0 length(datcell)]);
        ylabel('nuc/cyto ratio, mean over time');
        xlabel('Positions');
        title(['microCol of size ' num2str(colSZ) ],'fontsize',15);
        % ncells for each colony , per frame
            figure(10),subplot(1,3,3),plot(ncells{ii},'*','color',colors(ii,:),'markersize',5);
            hold on
            ylim([0 5])
            ylabel('Number of cells in the colony','fontsize',9);
        
    end
    end
end
end