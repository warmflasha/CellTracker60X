%%
%ilastikfile = ('/Users/warmflashlab/Desktop/IlastikMasks_headless_PluriW0/NucMaskPluri_tg56.h5');
% ilastikfile = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/03-02-2016dataTrainingOutput/nuc_mask_{P}.h5');%('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/NucMsks3D/NucMasks3Djan8set1_z3.h5');
% ilastikfile2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/trainingset/cyto_mask_15_z1.h5');
% ilastikfile = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/trainingset/nuc_mask_15_z1.h5');

ilastikfile2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/July7TilingLCellilastik_CytoMasks/cytomask3DJuly7set1_z0.h5');%/cyto_mask_z1.h5'
ilastikfile = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/July7TilingLCellilastik_NucMasks/nucmask3DJuly7set1_z0.h5');%/nuc_mask_z1.h5'

nuc = h5read(ilastikfile,'/exported_data');
cyto = h5read(ilastikfile2,'/exported_data');

   k=95;% 41 15,14,16; ,42 % 43,44 , 52,53,54- no
    nuc1 = nuc(1,:,:,k);% for probabilities exported % for the nuc dataset (for the tiled expteimtne) the labels are reversed
    nuc1 = squeeze(nuc1);
    mask1 = nuc1;
    
    mask3 = imfill(mask1 > 0.9,'holes');
    mask3 = bwareafilt(mask3,[850 850*10]);
    %imshow(mask3);
    
     cyto = cyto(2,:,:,k);% for probabilities exported
     cyto = squeeze(cyto);
     mask2 = cyto;
     
    % look at probability masks output
    figure(1), subplot(1,2,1),imshow(mask1');
    %figure(1), subplot(1,3,2),imshow(mask2');
    figure(1), subplot(1,2,2),imshow(mask3');
   
    Lnuc = mask3';%im2bw(mask1,0.5);
    %Lcyto = im2bw(mask2,0.85);
    Lcyto = imfill(mask2>0.7,'holes');
   % Lcyto = bwareafilt(Lcyto1,[0 50000]);
    figure(2),subplot(1,2,1), imshow(Lnuc);
   figure(2),subplot(1,2,2), imshow(Lcyto'&~Lnuc);
    
   
   %%
   [MaskFin2] = Unmergetwonuclei(mask3);
   %[MaskFin3] = Unmergetwonuclei(MaskFin2);
   figure, imshow(MaskFin2);
   
%%
% obtain peaks for one frame 
% these functions are for the 3D segmentation  ( can have one plane only
% anyway)

% last data ( July 26, 2016, 40X)
ilastikdircyto = ('/Users/warmflashlab/Desktop/July26Tiling2_CytoMasks');
ilastikdirnuc =('/Users/warmflashlab/Desktop/July26Tiling2_NucMasks');
imagedir1 =('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-26-Tiling2_24hrtotalBMP410ngml/uColTiling2_20hrinBMP4_20160726_24717PM/nucrawdata_z2z3');
imagedir2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-26-Tiling2_24hrtotalBMP410ngml/uColTiling2_20hrinBMP4_20160726_24717PM/cytorawdata_z2z3');

timegroup = [];%1

chan = [1 2];
chanal = [1 2];
paramfile = 'setUserParamLiveImagingAN_40X';       %'setUserParamLiveImagingAN_40X   setUserParamLiveImagingAN';
paramfile3D = 'setUserParam3DsegmentationAN_40X';  %'setUserParam3DsegmentationAN    setUserParam3DsegmentationAN_40X'

outfile = 'testtest.mat';%3Dsegm_febdata

  positions = [1]; %july 26 tiling 2 position (40X)
  %positions = [0 2 3 4 5 8 11 15 16 17 18 20 21 22 25 26 29 30 32 33 10]; %  february 3 dataset
    %positions = [5 7 8 9 11 13 15 16 17 18 20 21 23 25 27 28 29 33 34 35 36 37 38 39 40];
 % July 7 dataset 40 poszitions therelastpart

pl =2;%2,3 5
strcyto = 'cytomask'; %CytoMasks3Djan8set  %cytomask3DFebset % nucmask3DJuly7setcyto_mask_
strnuc = 'nucmask';   %NucMasks3Djan8set   %nucmask3DFebset %nuc_mask_cytomask3DJuly7set
k = 1;
%for k = 1:size(positions,2)
    pos = positions(k);
    rundataset3D(ilastikdirnuc,ilastikdircyto,imagedir1,imagedir2,pos,paramfile,timegroup,outfile,paramfile3D,pl,strnuc,strcyto,chanal);
%end


%%
% run tracking and colony grouping on the live cell dataset
strdir = '39_40X_imprBGandSegm.mat';%[33 35 36 39 40] 
paramfiletrack = 'newTrackParamAN';
TrackGroupuCol(strdir,paramfiletrack);   %TrackGroupuCol60X
%%
outfile ='37_40X_imprBGandSegm.mat';%1_test_3DsegmJul26data40x
trajmin =1;
plotcelltraces(outfile,trajmin)
%%
n1 = 1;
n2 = 10;m = 1;matfile = '0_test2.mat';
coltoplot = 1;
MatchTrajectorywithColony(ilastikdirnuc,ilastikdircyto,imagedir1,imagedir2,positions,pl,strnuc,strcyto,m,n1,n2,matfile,coltoplot)

%%
% plot nuc and cyto masks, colorcoded for label matrix elements
N =20;
n =uncompressBinaryImg(imgfiles(N).compressNucMask);
nc = uncompressBinaryImg(imgfilescyto(N).compressNucMask);
%close all
bwn = bwlabel(n);
bwn2 = bwconncomp(bwn);
a2 = label2rgb(bwn);
bw2 = bwlabel(nc);
bw3 = bwconncomp(bw2);
a = label2rgb(bw2);
figure(28),subplot(1,2,1),imshow(a2);
hold on
subplot(1,2,2),imshow(a);


                 
%%
% Get means over time but separately for different colony sizes
% new figure for each new colony size ( data drawn from multiple .mat
% files)

fr_stim = 16 ;        %22(jan8data)  %38    %16 (feb16 and     july7 data)   july 26 data (fr_stim = 12)
delta_t = 15;%15       % 12 min                 15min           17min        (17min)
p = fr_stim*delta_t/60;
timecolSZ = fr_stim;%10
p2 = (timecolSZ)*delta_t/60;
cmap = colorcube;close all
cmap2 = hot;close all
sz = 81;%81 99 83 
 strdir = '*_60X_testparam_allT.mat';
ff = dir(strdir);%'*_60X_testparam_allT.mat'ff = dir('*60Xjan8_R*.mat');
clear traces 
clear traces_one
clear traces_two
clear traces_three
totalcol = zeros(6,1);
q = 1;
r = 1;
p = 1;

for k=1:length(ff)
    outfile = ff(k).name; %nms{k};
    
    % outfile = ('12_3D_20hr_test_xyz.mat');
    
    load(outfile,'colonies','peaks');
    if ~exist('colonies','var');
        continue
    end
    tps = length(peaks);
    numcol = size(colonies,2); % how many colonies were grouped within the frame
    traces = cell(1,numcol);
    
    for j = 1:numcol
        colSZ = colonies(j).numOfCells(timecolSZ-1); % colony size determined at the time of stimulation
        
        traces{j} = colonies(j).NucSmadRatio(:);
        traces{j}(traces{j}==0) = NaN;
        %(traces{j}(traces{j}(1:fr_stim,1)>1.2) = NaN;           % these cells already have very high ratio value in the first couple frames, most likely junk
       %traces{j}(traces{j}>1.4) = NaN;       % these outliers are dead cells ot junk
       % traces(cellfun(@isempty,traces)==1)=[];
       if colSZ>0 && colSZ<5  
           
           for h = 1:size(traces{j},2)
               if length(traces{j}(isnan(traces{j}(:,h))==0))>75       % FILTER OUT SHORT TRAJECTORIES
                 totalcol(colSZ) = totalcol(colSZ)+1;
                   figure(colSZ), plot(traces{j}(:,h),'*','color',cmap(k,:));hold on% cmap(j,:) 'r' traces
                   ylim([0 2]);
                   xlim([5 sz]);
                   ylabel('mean Nuc/Cyto smad4 ');
                   xlabel('time, hours');
                   title(['All microColonies of size ' num2str(colSZ) ]);
                   %text(60,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] );
                   
                   traces{j}(isnan(traces{j})==1) = 0;
                   d =  size(traces{j},2);
                   if colSZ == 1
                       traces_one{q}= traces{j}(:,h);
                   end
                   if colSZ == 2
                       traces_two{r}= traces{j}(:,h);
                   end
                   if colSZ >2% == 3 || colSZ == 4
                       traces_three{p}= traces{j}(:,h);
                   end
                   %             if colSZ == 4
                   %             traces_four{p}= traces{j};
                   %             end
               end
               q = q+1;
               p = p+1;
               r = r+1;
           end
       end
    end
end

figure, plot(1:size(totalcol,1),totalcol,'r-*','markersize',18,'linewidth',3);
xlabel('cells per colony','fontsize',20);
ylabel('totla colonies','fontsize',20);
title('colony size distribution','fontsize',20)

%%
%clear findatanew
%findatanew = cell(1,3);
traces_curr = traces_one;
colSZ =1;
traces_curr(cellfun(@isempty,traces_curr)==1)=[];

%findatanew = cell(1,3);
d = size(traces_curr,2);
clear replace
clear sm

sz = 81;
%sz = 81; %81                % see if any positions have traces end earlier than the latest time point in peaks (sz) 
for k=1:d
a =size(traces_curr{k},1);
if a < sz && a ~= sz;%a ~= sz
sm(k) = a;

 else
     sm(k) = sz;
 end
end
b = find(sm < sz);
if ~isempty(b)
    % need to do this only if above , some traces were ended earlier than sz
%     a = find(sm>0);
%     sm1 = nonzeros(sm);
    replace = cell(1,sz);
    for jj=1:size(b,2)%nonzeros(sm)
        replace{jj} = zeros(sz,size(traces_curr{b(jj)},2));
        replace{jj}(1:size(traces_curr{b(jj)},1),1) = traces_curr{b(jj)}(:,1);
        traces_curr{b(jj)} = replace{jj};
    end
end
clear traces_one_new
q = 1;
for k=1:d
s = size(traces_curr{k},2);
if s<=30
traces_one_new(:,q:q+s-1) = traces_curr{k}(:,:);%  smoothinditrace(traces_curr{k}(:,:),4) smoothtrace = smoothinditrace(vect,N);N = 4, smoothwindow
q = q+1;
end
end

fin_data = zeros(sz,2);

for j =1:size(traces_one_new,1)
    for k=1:size(traces_one_new,2)
        %if size(nonzeros(traces_one_new(:,k)),1)>60;% don't include short traces into the analysis
        fin_data(j,1) = mean(nonzeros(traces_one_new(j,:)));
        fin_data(j,2) = std(nonzeros(traces_one_new(j,:)));
        %end
    end
end
findatanew{colSZ} = fin_data;
%%

save('60X_long2.mat','findatanew');

figure(1), hold on
plot(vect1',fin_data(:,1),'-*r','linewidth',3)

%%
% smooth individual  trace
vect = traces_curr{8};
smoothtrace = smoothinditrace(vect,N);
figure, plot(vect); hold on; plot(smoothtrace,'*');


%%
load '60X_long2.mat'%meanTraj_SDrmjunktraces
fr_stim = 16 ;        %22(jan8data) %38 %16(feb16 and july7 data)    12july26 data
delta_t = 15;%15       % 12 min              15min           17min   17 min
p = fr_stim*delta_t/60;
timecolSZ = fr_stim;%10
p2 = (timecolSZ)*delta_t/60;
sz = 81;%81 

pfin = (sz)*delta_t/60;
cmap = colorcube;close all
cmap2 = hot;close all
 
 vect1 = (1:sz);
 startplot = 4;   %        START FROM TIME POITN5 SINCE THERE IS SOME JUNK IN THE FIRST TPT THAT MAKES THE SIGNAL JUMP UP
vect = (startplot:sz)*(delta_t)/60;
colN = {'b','g','r'};

for k=1:size(findatanew,2)
figure(10), errorbar(vect',findatanew{k}(startplot:end,1),findatanew{k}(startplot:end,2),'color',colN{k},'marker','*');hold on%findatanew{k}(:,2)
figure(10),title('All microColonies');
legend('1-cell colonies','2-cell colonies','3-cell colonies')
text(40,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] ); 
xlim([0.5 pfin]);
ylim([0.5 2]);
ylabel('mean Nuc/Cyto smad4 ');
xlabel('time, hours');
figure(5), plot(vect',findatanew{k}(startplot:end,1),'color',colN{k},'marker','*'); hold on
figure(6), plot(vect',findatanew{k}(startplot:end,2),'color',colN{k},'marker','*');  hold on%smoothVAR(findatanew{k}(startplot:end,2),2)    power(findatanew{k}(startplot:end,2),2)
end

figure(5),title('Mean Trajectories');
legend('1-cell colonies','2-cell colonies','3- and 4-cell colonies')
text(40,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] ); 
xlim([0.5 pfin]);
ylim([0.5 2]);
ylabel('mean Nuc/Cyto smad4 ');
xlabel('time, hours');
figure(6),title('Variance');
text(40,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] );
legend('1-cell colonies','2-cell colonies','3- and 4-cell colonies');
xlim([0.5 pfin]);
ylim([0 0.7]);
ylabel('Variance ');
xlabel('time, hours');
% figure(10+k),title('All microColonies');
% legend('1-cell colonies','2-cell colonies','3- and 4-cell colonies')
% text(40,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] ); 
% xlim([0.5 pfin]);
% ylim([0.5 1.8]);
% ylabel('mean Nuc/Cyto smad4 ');
% xlabel('time, hours');


  %%
  % relabel the x axis into the correct time units
%M = number of different colony sizes
M = 3;

for ii = 1:M
            h = figure(ii);
            mm = max(h.Children.XTick);
            for k=1:mm
            tnew(k) = k*delta_t/60;
            end
             h.Children.XTick = [(1:12:mm)];
             h.Children.XTickLabel = [(tnew(1:12:mm))];
end
%%
% NEW: get the before and after  plots
clear all;
fr_stim = 22;% 16 22 %38
delta_t = 12; % 12% in minutes
p = fr_stim*delta_t/60;
global userParam;
resptime =72;                            % in frames, count starts after (fr_stim+resptime) frames
range = [26 20];                         % if want to look at average signaling between specific frames
jumptime = 5;                            % number of frames to achieve response to ligand
p2 = (resptime+jumptime)*delta_t/60;
coloniestoanalyze = 3;
cmap = summer;
C = {'b','r','g','m','c'};
for jj = 1:5
q(jj) = 1;
w(jj) = 1;
s(jj) = 1; 
end
    
    ff = dir('*_test*.mat');       %jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs %Outfile  %dec31_set_Diff
    for k=1:size(ff,1)
        outfile = ff(k).name ;
        load(outfile,'colonies','peaks');
        tps = length(peaks);
        numcol = size(colonies,2); % how many colonies were grouped within the frame
        for j = 1:numcol
            if ~isempty(fr_stim)   % for the case when                  
                          
                colSZ = colonies(j).numOfCells(fr_stim); % how many cells within colony at time fr_stim
                              
                if colSZ >0 && (colSZ<4)
                    stats =  colonies(j).DynNucSmadRatio(tps,fr_stim,resptime,range,jumptime);%,resptime
                    bf = (stats(:,1));
                    aft =(stats(:,2));
                    window = (stats(:,3));
                    jump = (stats(:,4));
                    
                    %b = find(isnan(bf));
                    bf(isnan(bf))=0;
                    aft(isnan(aft)) = 0;
                                       
                    w = find(isnan(window));
                    window(w)=0;
                    
                    currlengthbf = size(nonzeros(bf),1);
                    currlengthaft = size(nonzeros(aft),1);
                    currlengthgap = size(nonzeros(window),1);
                    
                    if length(nonzeros(bf)) ~= currlengthbf
                        disp(['Here: file ' num2str(k) ' colony ' num2str(j)] );
                    end
                    
                    bf_fin{colSZ}((q(colSZ):(q(colSZ)+currlengthbf)-1),1) = nonzeros(bf);
                   % disp([num2str(q(colSZ)) '   ' num2str(currlengthbf) '   ' num2str(bf_fin{colSZ}')]);
                    aft_fin{colSZ}((s(colSZ):(s(colSZ)+currlengthaft)-1),1) = nonzeros(aft);
                    %disp([num2str(s) '   ' num2str(currlengthaft) '   ' num2str(aft_fin{colSZ}')]);
                    window_fin{colSZ }((w(colSZ):(w(colSZ)+currlengthgap)-1),1) = nonzeros(window);
                    %disp([num2str(w) '   ' num2str(currlengthaft) '   ' num2str(window_fin{colSZ }')]);
                                        
                    q(colSZ) = q(colSZ)+currlengthbf;
                    w(colSZ) = w(colSZ)+currlengthgap;
                    s(colSZ) = s(colSZ)+currlengthaft;
                    
                    
                    hold on,figure(4), plot(bf,aft,'*','color',C{colSZ},'markersize',15);% cmap(colSZ,:,:)
                    title(['mean Nuc/Cyto smad4 ' num2str(p2) 'hrs after bmp4']);
                    ylim([0 2]);
                    hold on,figure(5),subplot(1,2,1), plot(colSZ,bf,'*','color',C{colSZ},'markersize',15);
                    ylim([0 2]);
                    hold on,figure(5),subplot(1,2,2), plot(colSZ,aft,'*','color',C{colSZ},'markersize',15);%cmap(colSZ*16,:,:)
                    ylim([0 2]);
                    hold on,figure(10), plot(colSZ,window,'*','color',C{colSZ},'markersize',15);%bf,window
                    title(['mean Nuc/Cyto smad4 between ' num2str((range(1)*delta_t)/60) 'and' num2str((range(2)*delta_t)/60) 'hours']);
                    ylim([0 2]);
                    xlim([0 (coloniestoanalyze+1)]);
                    hold on,figure(11), plot(colSZ,jump,'*','color',C{colSZ},'markersize',15);%amplitude of the actual jump during the jumptime
                    title(['Amplitude of the jump, bf to ' num2str(((fr_stim+jumptime+range(2))*delta_t)/60) 'hours']);
                    ylim([0 2]);
                    xlim([0  (coloniestoanalyze+1)]);
                                       
                    %  end
                    save('meansdiff.mat','bf_fin','aft_fin','window_fin');
                end
            end
            
            if isempty(fr_stim)   % for the pluri case
                %if colSZ == colonies(j).numOfCells(10);%
                colSZ = colonies(j).numOfCells(10);        
                if colSZ > 0
                    stats =  colonies(j).DynNucSmadRatio(tps,fr_stim,resptime,range,jumptime);%,resptime
                    bf = (stats(:,1));
                    window = (stats(:,3));
                    
                    bf(isnan(bf))=0;
                    window(isnan(window))=0;
                    
                    currlengthbf = size(nonzeros(bf),1);
                    currlengthgap = size(nonzeros(window),1);
                    
                    bf_fin{colSZ}((q:(q+currlengthbf)-1),1) = nonzeros(bf);%sum(bf);
                  %  disp([num2str(q) '   ' num2str(currlengthbf) '   ' num2str(bf_fin{colSZ }')]);
                    window_fin{colSZ}((w:(w+currlengthgap)-1),1) = nonzeros(window);%sum(window);
                    
                    q = q+currlengthbf;
                    w = w+currlengthgap;
                    
                    hold on,figure(4), plot(colSZ,bf,'*','color',C{colSZ},'markersize',15);
                    ylim([0 2]);
                    xx = 0:1:coloniestoanalyze;
                    yy = ones(1,(coloniestoanalyze+1));
                    figure(4),plot(xx,yy,'-k','linewidth',2);
  
                end
                save('meanspluri.mat','bf_fin','window_fin');
            end
        end
            
        end
        
   % end
    
    
%end
if ~isempty(fr_stim)
                xx = 0:1:coloniestoanalyze;
                yy = ones(1,coloniestoanalyze+1);
                hold on,figure(4),plot(xx,xx,'-k','linewidth',2);
                ylim([0 2]);
                xlim([0 2]);
                hold on,figure(10),plot(xx,xx,'-k','linewidth',2);
                ylim([0 2]);
                xlim([0 2]);
                ylabel(['mean Nuc/Cyto smad4  ' num2str(p2) ' hours after stimulation']);
                xlabel('mean Nuc/Cyto smad4 Before stimulation');
                hold on,figure(5),subplot(1,2,2),plot(xx,yy,'-k','linewidth',2);
                ylabel(['mean Nuc/Cyto smad4  ' num2str(p2) ' hours after stimulation']);
                ylim([0 2]);
end

%%
% plot the mean values
jumptime = 5;% in frames
fr_stim = 22;%22 %38 %16
fldat = [2 3];
delta_t = 12; % 12% in minutes
p = fr_stim*delta_t/60;
global userParam;
userParam.colonygrouping = 120;
flag = 1;
test = cell(1,20);
resptime =75;% 15 50 36 in frames ( converted to hours later)
p2 = (resptime+jumptime)*delta_t/60;
coloniestoanalyze = 3;
cmap = parula;
C = {'b','g','r','m','m'};
%test = cell(1,coloniestoanalyze);
for jj = 1:coloniestoanalyze
    colSZ = jj;
    if isempty(fr_stim)
        load('meanspluri.mat');
        before(jj) = mean(nonzeros(bf_fin{jj}));
               
        errbf(jj) = std(nonzeros(bf_fin{jj}));
        
        figure(6), errorbar(jj,before(jj),errbf(jj),'*','markersize',15,'color',C{jj});hold on
        ylabel('Mean Nuc/Cyto smad4');
        ylim([0.3 1.9]);
        xlim([0 (coloniestoanalyze+1)]);
    end
    if ~isempty(fr_stim)
        load('meansdiff.mat');
       aft_fin{jj}(aft_fin{jj} == inf)=0;
       bf_fin{jj}(bf_fin{jj} == inf)=0;
        before(jj) = mean(nonzeros(bf_fin{jj})); %sum(bf_fin{jj}(:,1))/sum(bf_fin{jj}(:,2));
        after(jj) = mean(nonzeros(aft_fin{jj}));%sum(aft_fin{jj}(:,1))/sum(aft_fin{jj}(:,2));
            
        errbf(jj) = std(nonzeros(bf_fin{jj}));
        erraft(jj) = std(nonzeros(aft_fin{jj}));
   
        figure(7), errorbar(jj,before(jj),errbf(jj),'*','markersize',15,'color',C{jj});hold on;
        figure(7), errorbar(jj,after(jj),erraft(jj),'.','markersize',30,'color',C{jj});hold on;
        ylim([0.3 1.9]);
        xlim([0 (coloniestoanalyze+1)]);
        figure(8), errorbar(jj,after(jj),erraft(jj),'.','markersize',30,'color',C{jj});hold on;
         legend('before','after');
         title(['mean Nuc/Cyto smad4  ' num2str(p2) ' hours after stimulation']);
        ylabel('Mean Nuc/Cyto smad4');
        ylim([0.3 1.9]);
        xlim([0 (coloniestoanalyze+1)]);
        
    end
end

      
      
      
