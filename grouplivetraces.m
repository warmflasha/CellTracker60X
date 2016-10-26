%% split traces based on the last time point signaling level for the same colony size
% Get means over time but separately for different colony sizes
% new figure for each new colony size ( data drawn from multiple .mat
% files)

fr_stim = 22;        %22(jan8data)  %38    %16 (feb16 and     july7 data)   july 26 data (fr_stim = 12)
delta_t = 17;%15       % 12 min                 15min           17min        (17min)
p = fr_stim*delta_t/60;
trajmin = 30;%30
tpt = 99;%81 99 83
strdir = '*_out.mat';%_outFebsetBGan
ff = dir(strdir);%'*_60X_testparam_allT.mat'ff = dir('*60Xjan8_R*.mat');
C = {'b','r'};
signalbin = [];
signthresh = 1;
clear traces
clear traces_one
clear traces_two
clear traces_three
clear a
totalcol = zeros(6,1);
q = 1;
nc = 1;                       % look at nc-cell colonies
N = 20;                       % how many last time point to average in order to sort into beans
for k=1:length(ff)
    outfile = ff(k).name; %nms{k};
    load(outfile,'colonies');
    if ~exist('colonies','var');
        disp('does not contain colonies structure')
    end
    numcol = size(colonies,2); % how many colonies were grouped within the frame
    traces = cell(1,numcol);
    
    for j = 1:numcol
        if size(colonies(j).ncells_actual,1)>fr_stim  %&& (any(colonies(j).onframes) == tpt)             % new segmentation
            colSZ =colonies(j).ncells_actual(fr_stim) ;
            % new segmentation
            %colSZ =colonies(j).numOfCells(timecolSZ-1) ;                                                      % for old mat files analysis
            if colSZ == nc
                jj =1;
                
                traces{j} = colonies(j).NucSmadRatio;                                                        % colonies(j).NucSmadRatio(:)
                % traces{j} = colonies(j).NucSmadRatioOld;                                                     % for old mat files analysis
                sz = size(traces{j},2);
                for h = 1:size(traces{j},2)
                    [r,~] = find(isfinite(traces{j}(:,h)));                  %
                    dat = zeros(tpt,1);
                    dat(r,1) = traces{j}(r,h);
                    s = mean(nonzeros(dat(end-N:end)));                       % take the mean of the last chunk of N time points
                    if length(nonzeros(dat))>trajmin && (s < signthresh);      % mean(nonzeros(dat(end-15:end))))  FILTER OUT SHORT TRAJECTORIES  and select the ones that end low in signaling
                        jj =1;
                        disp(['filter trajectories below' num2str(trajmin)]);
                        disp(['use' num2str(length(nonzeros(dat)))]);
                        figure(jj), plot(dat,'-*','color',C{jj});hold on         % here plot the traces that met the condition
                        signalbin{jj}(:,q+sz-1) = dat;                           % here store the traces which meat condition
                        % disp(q+sz-1)
                        xx = size(dat,1)-1;
                        yy = dat(end-1);%traces{k}(end,h);
                        text(xx,yy,[num2str(s)],'color','m','fontsize',11);%
                        figure(jj) ,hold on
                        ylim([0 2.5]);
                        xlim([0 (tpt+10)]);
                        ylabel('mean Nuc/Cyto smad4');
                        xlabel('frames');
                    end
                     if length(nonzeros(dat))>trajmin && (s >= signthresh);      % FILTER OUT SHORT TRAJECTORIES  and select the ones that end low in signaling
                        jj =2;
                        disp(['filter trajectories below' num2str(trajmin)]);
                        disp(['use' num2str(length(nonzeros(dat)))]);
                        figure(jj), plot(dat,'-*','color',C{jj});hold on         % here plot the traces that met the condition
                        signalbin{jj}(:,q+sz-1) = dat;                          % here store the traces which meat condition
                        % disp(q+sz-1)
                        xx = size(dat,1)-1;
                        yy = dat(end-1);%traces{k}(end,h);
                        text(xx,yy,[num2str(s)],'color','m','fontsize',11);%
                       figure(jj) ,hold on
                       ylim([0 2.5]);
                       xlim([0 (tpt+10)]);
                       ylabel('mean Nuc/Cyto smad4  ');
                       xlabel('frames');
                    end
                    
                end
                q = q+sz;
                
            end
        end
    end                       % new segmentation
end
%end
% figure, plot(1:size(totalcol,1),totalcol,'r-*','markersize',18,'linewidth',3);
% xlabel('cells per colony','fontsize',20);
% ylabel('totla colonies','fontsize',20);
% title('colony size distribution','fontsize',20)
%% average those trajectoris

vect = (1:tpt)';
binmean = zeros(tpt,2);
err =zeros(tpt,2); 
for j =1:2                  % remove Nans
for k=1:size(binmean,1)
    for jj=1:size(signalbin{j},2)
   if (isfinite(signalbin{j}(k,jj))==0) || signalbin{j}(k,jj)< 0.5 || signalbin{j}(k,jj)>1.85 % to remove signaling values that come from long-traced junk
       signalbin{j}(k,jj) = 0;
   end
    end
end
end
% average over cells
for j =1:2% loop over the bins;         
for k=1:size(binmean,1)
    binmean(k,j) = mean(nonzeros(signalbin{j}(k,:)));   % mean over nonzero values of signaling at each time point
    err(k,j) = std(nonzeros(signalbin{j}(k,:)));
end
end
%% plot means for each bin in signaling
vect = (1:tpt)';
colormap = colorcube;
C = {'c','r'};
label = {'signal below','signal above'};
b = [signthresh];

for j = 1:2
figure(11), errorbar(binmean(:,j),err(:,j),'-.','color',C{j},'linewidth',1.5); hold on%colormap(j+5,:)
ylim([0.3 1.8]);
xlim([0 105])
ylabel('mean Nuc/Cyto smad4  ');
xlabel('frames');
text(vect(end)+0.2*j,binmean(end,j)+0.1*j,[ label(j) num2str(b(1))],'color',C{j},'fontsize',20);%['mean ColCdx2 ' num2str(colonies2(j).cells(h).fluorData(1,end))]
title('One-cell colonies','fontsize',20);
if nc == 2
title('Two-cell colonies','fontsize',20);
end

figure(12), plot(vect,binmean(:,j),'-.','color',C{j},'linewidth',2); hold on%colormap(j+5,:)
ylim([0.3 1.8]);
xlim([0 115])
ylabel('mean Nuc/Cyto smad4  ');
xlabel('frames');
text(vect(end)+0.2*j,binmean(end,j)+0.1*j,[ label(j) num2str(b(1))],'color',C{j},'fontsize',20);%['mean ColCdx2 ' num2str(colonies2(j).cells(h).fluorData(1,end))]
title('One-cell colonies','fontsize',20);
if nc == 2
title('Two-cell colonies','fontsize',20);
end
end
