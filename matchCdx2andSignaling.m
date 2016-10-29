%% get the centroids of the live data
% load('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/registeredDAPI.mat','colfixall','xyall','fluordata'); %load the data with coordinates of the fixed colonies (centroid of the colony)
% load('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/alignWithDapi.mat');              % load acoords which matched live and fixed data after reimaging with dapi 
% load('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/tformwithDAPI.mat');
load('/Volumes/data2/Anastasiia/LiveCellImagingGFPs4RFPh2b/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/registeredDAPI.mat','colfixall','xyall','fluordata');
load('/Volumes/data2/Anastasiia/LiveCellImagingGFPs4RFPh2b/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/alignWithDapi.mat'); 
load('/Volumes/data2/Anastasiia/LiveCellImagingGFPs4RFPh2b/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/tformwithDAPI.mat');

dirlive = ('/Volumes/data2/Anastasiia/LiveCellImagingGFPs4RFPh2b/2016-10-19-LIVECELLanalysisLatest/new_outfiles_tiling1');

%dirlive = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-10-17-projections/new_outfiles_tiling1');%new_outfiles_tiling1anBG new_outfiles_tiling1
%t = imread('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/testwithDAPI.tif');
t = imread('/Volumes/data2/Anastasiia/LiveCellImagingGFPs4RFPh2b/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/testwithDAPI.tif');
figure(1), imshow(t(:,:,1),[0 5e3]); hold on;

mytform = fitgeotrans(movingPoints, fixedPoints, 'affine');
toTranslate = [-250,50];    
midpoint = floor(size(t)/2);
rotmat = mytform.T(1:2,1:2);

colormap = prism;
positions = (0:39);
strdir = '_out.mat';% _out tile1BGan 
last = 100;             % which time point in live dataset to check before matching

q = 1;

for k=1:size(positions,2)
    outfile = [ dirlive '/' num2str(k-1) strdir ];
    load(outfile,'colonies')
    if ~isempty(colonies)
        for j=1:size(colonies,2)
            posnow = [];
            for h=2:size(colonies(j).cells,2)                               % start from second, since the first one is empty , as returned by the new analysis
                if colonies(j).cells(h).onframes(end) == (last)
                                %figure(q); hold on;
                    posnow = [posnow; colonies(j).cells(h).position(end,:)];
                end
            end
            if ~isempty(posnow)
            datatogether(q).colony = colonies(j);
            posnow = bsxfun(@plus,posnow,acdapi(k).absinds([2 1]));
            posnow = bsxfun(@minus,posnow,mean(xyall));
            posnow = posnow*rotmat;
            posnow = bsxfun(@plus,posnow,midpoint([2 1])+toTranslate);
            
            colnow_center = mean(posnow,1);
            dists = bsxfun(@minus,colfixall,colnow_center);
            dists = sqrt(sum(dists.*dists,2));
            [~, ind]=min(dists);
            datatogether(q).fixedData = fluordata(ind,:);
            %disp(outfile);
            datatogether(q).outfiles = outfile; % save the outfile that was matched
            figure(1), hold on
            plot(colnow_center(1),colnow_center(2),'.','color',colormap(j+10,:,:),'markersize',10);%
            plot(colfixall(ind,1),colfixall(ind,2),'*','color',colormap(j+10,:,:),'markersize',10);%
            q= q + 1;
            end
        end
        
    end
    
end

%save('registeredDAPInewTraces','datatogether','-append');
%% plot traces with fixed data
close all
%load('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/registeredDAPInewTraces.mat');
%load('/Volumes/data2/Anastasiia/LiveCellImagingGFPs4RFPh2b/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/registeredDAPInewTraces.mat');
positions = (0:39);
clear tracesbycol
trajmin = 40; % 
fr_stim = 16; % needed to check the colony size at this point
traces = [];  % initialize
cdx2val = []; % initialize for Cdx2 values
cdx2todapi = [];
finsign=[];
tracesbybin = cell(1,2);   % two bins, same colony size
binSZ = [1];               
C = {'b','r'};
nc = 1;                                                              % colony size to be plotted
clear dat
clear jj
q = 1;
N = 15;
for k=1:size(datatogether,2)   %  loop over colonies
    
    if size(datatogether(k).colony.ncells_actual,1)>fr_stim;          % if the cell traced at least untill stimulation
        colSZ =datatogether(k).colony.ncells_actual(fr_stim) ;          % check colony size
        cdx2todapi(k) = datatogether(k).fixedData(:,3)/datatogether(k).fixedData(:,1); % get the value of cdx2dapi for this colony
        
        if colSZ == nc && (cdx2todapi(k)<= binSZ(1))
            jj =1;
            traces{k} = datatogether(k).colony.NucSmadRatio;              % get the nuclea2smad ratio for all traces in this colony
            traces{k}((traces{k} == 0)) = nan;                            % put nans instead of zeros, to avoid merging theminto consecutive time points and ahe plotting wrong
            sz = size(traces{k},2);                                       % how many traces are in the colony
            for h = 1:size(traces{k},2)                                  % loop over traces
                [r,~] = find(isfinite(traces{k}(:,h)));                  %
                dat = zeros(size(traces{k},1),1);
                dat(r,1) = traces{k}(r,h);                              % dat contains values of signling corresponding to their frame numbers, rest is zeros
                if length(nonzeros(dat))>trajmin                         % here filter out short trajectories
                    disp(['filter trajectories below' num2str(trajmin)]);
                    disp(['use' num2str(length(nonzeros(dat)))]);
                    figure(jj), plot(dat,'-*','color','b');hold on         % here plot the traces that met the condition
                    tracesbybin{jj}(:,q+sz-1) = dat;                      % here store the traces which meat condition
                    cdx2val{jj}(:,q+sz-1)= datatogether(k).fixedData(:,3)/datatogether(k).fixedData(:,1); % here store the cdx2todapi values which meat condition
                    finsign{jj}(:,q+sz-1) = mean(nonzeros(dat((end-N):end)));
                    % disp(q+sz-1)
                end
                
           end
            
            q = q+sz;
            xx = size(dat,1);
            yy = dat(end);%traces{k}(end,h);
            text(xx,yy,[num2str(cdx2todapi(k)) ', outfile' num2str(datatogether(k).outfiles(end-9:end-8))],'color',C{jj},'fontsize',11);%['mean ColCdx2 ' num2str(colonies2(j).cells(h).fluorData(1,end))]
            figure(jj) ,hold on
            ylim([0 2.5]);
            xlim([0 115])
            ylabel('mean Nuc/Cyto smad4  ');
            xlabel('frames');
             
        end
        
        if colSZ == nc && (cdx2todapi(k)> binSZ(1))                 % here do the same as above but for colonies of size nc and cdx2 values above the binsz(1)
          jj = 2;
            traces{k} = datatogether(k).colony.NucSmadRatio;              % get the nuclea2smad ratio for all traces in this colony
            traces{k}((traces{k} == 0)) = nan;                            % put nans instead of zeros, to avoid merging theminto consecutive time points and ahe plotting wrong
            sz = size(traces{k},2);                                       % how many traces are in the colony
            for h = 1:size(traces{k},2)                                  % loop over traces
                [r,~] = find(isfinite(traces{k}(:,h)));                  %
                dat = zeros(size(traces{k},1),1);
                dat(r,1) = traces{k}(r,h);                              % dat contains values of signling corresponding to their frame numbers, rest is zeros
                if length(nonzeros(dat))>trajmin                         % here filter out short trajectories
                    disp(['filter trajectories below' num2str(trajmin)]);
                    disp(['use' num2str(length(nonzeros(dat)))]);
                    figure(jj), plot(dat,'-*','color','r');hold on         % here plot the traces that met the condition
                    
                    tracesbybin{jj}(:,q+sz-1) = dat;                      % here store the traces which meat condition
                    cdx2val{jj}(:,q+sz-1)= datatogether(k).fixedData(:,3)/datatogether(k).fixedData(:,1); % here store the cdx2todapi values which meat condition
                    finsign{jj}(:,q+sz-1) = mean(nonzeros(dat((end-N):end)));
                    % disp(q+sz-1)
                end
                
            end
            
            q = q+sz;
            xx = size(dat,1);
            yy = dat(end);
            text(xx,yy,[num2str(cdx2todapi(k)) ', outfile' num2str(datatogether(k).outfiles(end-9:end-8))],'color',C{jj},'fontsize',11);%display the cdx2 value and the outfile that the trace came from
            figure(jj) ,hold on
            ylim([0 2.5]);
            xlim([0 115])
            ylabel('mean Nuc/Cyto smad4  ');
            xlabel('frames');
        end
        
        
    end
end
if nc == 2
    tracesbybin2 = tracesbybin;
   % save('registeredDAPInewTraces','tracesbybin2','-append');%% save the
    %data for the two-cell clonies

end
%save('registeredDAPInewTraces','tracesbybin','binSZ','cdx2val','-append');%% average the trajectories that end up with variaous values of Cdx2
%% get the mean trajectories corresponding to each CDX2 bin, if run only this section, need to set variable nc
%load('/Volumes/data2/Anastasiia/LiveCellImagingGFPs4RFPh2b/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/registeredDAPInewTraces.mat','tracesbybin','tracesbybin2','binSZ','datatogether');

%load('registeredDAPI.mat','tracesbybin','tracesbybin2','binSZ','datatogether');
%load('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/registeredDAPInewTraces.mat','tracesbybin','binSZ','datatogether','cdx2val');

% tracesbybin       cell array that contains traces for the one cell
                  % colonies but separated into cell arrays according to the final Cdx2 value
                  
% binSZ             defines the Cdx2 range
nc = 1;
if nc == 2
tracesbybin = tracesbybin2; % for the twocellcolonies
end
vect = (1:100)';
binmean = zeros(100,2);
err =zeros(100,2); 

for j =1:2                  % remove Nans
for k=1:size(binmean,1)
    for jj=1:size(tracesbybin{j},2)
   if (isfinite(tracesbybin{j}(k,jj))==0) || tracesbybin{j}(k,jj)>1.5 || tracesbybin{j}(k,jj)< 0.5  % to remove signaling values that come from long-traced junk
       tracesbybin{j}(k,jj) = 0;
   end
    end
end
end

 % average over cells
for j =1:2% loop over the bins;         
for k=1:size(binmean,1)
    
    binmean(k,j) = mean(nonzeros(tracesbybin{j}(k,:)));   % mean over nonzero values of signaling at each time point
    err(k,j) = std(nonzeros(tracesbybin{j}(k,:)));
end
end
if nc == 2
    binmean2 = binmean;
    err2 = err;
   % save('registeredDAPInewTraces','binmean2','err2','-append');
end
%save('registeredDAPInewTraces','binmean','err','-append');
%% plot averaged trajectories corresponding to each CDX2 bin, if run only this section, need to set variable nc
%load('/Volumes/data2/Anastasiia/LiveCellImagingGFPs4RFPh2b/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/registeredDAPInewTraces.mat','binmean','err','binmean2','err2');
%load('registeredDAPInewTraces.mat','tracesbybin','tracesbybin2','binmean','binmean2','err','err2','binSZ','datatogether');

nc = 2;
vect = (1:100)';
colormap = colorcube;
C = {'c','r'};
label = {'CDX2 below','CDX2 above'};
b = [1];
if nc == 2
    binmean = binmean2;
    err = err2;
end

for j = 1:2
figure(5), errorbar(binmean(:,j),err(:,j),'-.','color',C{j},'linewidth',1.5); hold on%colormap(j+5,:)
ylim([0.3 1.8]);
xlim([0 105])
ylabel('mean Nuc/Cyto smad4  ');
xlabel('frames');
text(vect(end)+0.2*j,binmean(end,j)+0.1*j,[ label(j) num2str(b(1))],'color',C{j},'fontsize',20);%['mean ColCdx2 ' num2str(colonies2(j).cells(h).fluorData(1,end))]
title('One-cell colonies','fontsize',20);
if nc == 2
title('Two-cell colonies','fontsize',20);
end

figure(6), plot(vect,binmean(:,j),'-.','color',C{j},'linewidth',2); hold on%colormap(j+5,:)
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