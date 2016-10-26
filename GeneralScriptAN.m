%%
%ilastikfile = ('/Users/warmflashlab/Desktop/IlastikMasks_headless_PluriW0/NucMaskPluri_tg56.h5');
% ilastikfile = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/03-02-2016dataTrainingOutput/nuc_mask_{P}.h5');%('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/NucMsks3D/NucMasks3Djan8set1_z3.h5');
% ilastikfile2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/trainingset/cyto_mask_15_z1.h5');
% ilastikfile = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/trainingset/nuc_mask_15_z1.h5');

ilastikfile2 = ('/Users/warmflashlab/Desktop/July26Tiling2_CytoMasks/cytomask1_z0.h5');%/cyto_mask_z1.h5'
ilastikfile = ('/Users/warmflashlab/Desktop/July26Tiling2_NucMasks/nucmask1_z0.h5');%/nuc_mask_z1.h5'

nuc = h5read(ilastikfile,'/exported_data');
cyto = h5read(ilastikfile2,'/exported_data');
lblN = 2;
   k=47;% 41 15,14,16; ,42 % 43,44 , 52,53,54- no
    nuc1 = nuc(lblN,:,:,k);% for probabilities exported % for the nuc dataset (for the tiled expteimtne) the labels are reversed
    nuc1 = squeeze(nuc1);
    mask1 = nuc1;
    
    mask3 = imfill(mask1 > 0.82,'holes');
   % mask3 = bwareafilt(mask3,[850 850*10]);
    %imshow(mask3);
    
     cyto = cyto(lblN,:,:,k);% for probabilities exported
     cyto = squeeze(cyto);
     mask2 = cyto;
     Lcyto = imfill(mask2>0.85,'holes');
    % look at probability masks output
    figure(1), subplot(1,2,1),imshow(mask1');
    %figure(1), subplot(1,3,2),imshow(mask2');
    figure(1), subplot(1,2,2),imshow(mask3');
   
    Lnuc = mask3';%im2bw(mask1,0.5);
    %Lcyto = im2bw(mask2,0.85);
    Lcyto = imfill(mask2>0.85,'holes');
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

outfile = 'tiling2.mat';%3Dsegm_febdata

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
strdir = '0_testBacktocoord.mat';%[33 35 36 39 40] 
paramfiletrack = 'newTrackParamAN';
TrackGroupuCol(strdir,paramfiletrack);   %TrackGroupuCol60X
%%
k = 1;
dirlive = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/FinalOutfiles_completeTraces');
strdir = '_40X_imprBGandSegm.mat';
outfile = [ dirlive '/' num2str(k-1) strdir ];

trajmin =1;
plotcelltraces(outfile,trajmin)
%plotcelltracesandFixedData(outfile,trajmin)
hold on
%%
% correlations between the signaling cells in  2-cell colonies (cells within the same
% colony)
% sampling rate: 1/T (T seconds after each measitement is taken)

% fluctuations of signali of the cell should correlate with the
% fluctuations of the neighboring cell (in the same colony)
% these should not correlate, if the cells are within different colonies
%direc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-19-09-Correlations_analysis');
load('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-19-09-Correlations_analysis/2celltotestcorr.mat');
str = 'testcorr';
file = '2celltotestcorr.mat';
maxval2 = zeros(19,2);
for i=1:19 %( how many 2-cell colonies there are

var = [ str num2str(i) ] ;
h = load(file,var);
hh = struct2cell(h);
%varcor{i} = zeros(size(hh{1},1),2); % to store the correlations vector

delta_t = 17*60;% seconds
sr = 1/delta_t; % sampling rate
y1 = hh{1};
%correlations between the cells within the same colony
[X1,lag1] = xcorr(y1(:,1),y1(:,2));

varcor{i} = zeros(size(X1,1),2); % to store the correlations vector
varcor{i}(:,1) = X1;
varcor{i}(:,2) = lag1';
maxval = max(varcor{i}(:,1));
[r,~] = find(varcor{i}(:,1) == maxval);
maxval2(i,1:2) = varcor{i}(r,:);



figure(1), subplot(211),plot((0:numel(y1(:,1))-1)/sr,y1(:,1),'-*r'); hold on
plot((0:numel(y1(:,2))-1)/sr,y1(:,2),'-*b'); 

figure(1),subplot(212), plot(lag1/sr,X1,'.k');hold on
ylabel('Amplitude');
xlabel('Time,s');
title('Cross-correlation between two cells in same 2-cell colony');
end

figure(4), plot(maxval2(:,1),'b.','markersize',14);
mean(maxval2(:,1));
%%
% corr btw cells in different 2-cell colonies
load('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-19-09-Correlations_analysis/2celltotestcorr.mat');
str = 'testcorr';
file = '2celltotestcorr.mat';
maxval2diffcol = zeros(18,2);


for i=2:19 %( how many 2-cell colonies there are

var = [ str num2str(i) ] ;
h = load(file,var);
hh = struct2cell(h);

var2 = [ str num2str(1) ] ;% take the first cell 
h2 = load(file,var2);
hh2 = struct2cell(h2);
y2 = hh2{1};            % fixed cell within a different colony

delta_t = 17*60;% seconds
sr = 1/delta_t; % sampling rate
y1 = hh{1};
% here need to make sure that the y1 and y2 are the same length
s1 = size(y2,1);
s2 = size(y1,1);
if s1<=s2
    [X1,lag1] = xcorr(y2(:,2),y1(1:s1,1));
else
    [X1,lag1] = xcorr(y2(1:s1,2),y1(:,1));
end
varcor{i} = zeros(size(X1,1),2); % to store the correlations vector
varcor{i}(:,1) = X1;
varcor{i}(:,2) = lag1';
maxval = max(varcor{i}(:,1));
[r,~] = find(varcor{i}(:,1) == maxval);
maxval2diffcol(i,1:2) = varcor{i}(r,:);

% figure(3), subplot(211),plot((0:numel(y1(:,1))-1)/sr,y1(:,1),'-*r'); hold on
% plot((0:numel(y1(:,2))-1)/sr,y1(:,2),'-*b'); 

figure(3),subplot(212), plot(lag1/sr,X1,'k'); hold on
ylabel('Amplitude');
xlabel('Time,s');
title('Cross-correlation between two cells in different 2-cell colonies');
end

figure(2),hold on
plot(maxval2diffcol(:,1),'r.','markersize',14);
mean(maxval2diffcol(:,1));



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

fr_stim = 22 ;        %22(jan8data)  %38    %16 (feb16 and     july7 data)   july 26 data (fr_stim = 12)
delta_t = 17;%15       % 12 min                 15min           17min        (17min)
p = fr_stim*delta_t/60;
timecolSZ = fr_stim;%10
p2 = (timecolSZ)*delta_t/60;
cmap = colorcube;close all
cmap2 = hot;close all
trajmin = 30;%50
sz = 99;%81 99 83
strdir = '*_out.mat';%0_60X_testparam_allT
ff = dir(strdir);%'*_60X_testparam_allT.mat'ff = dir('*60Xjan8_R*.mat');
clear traces
clear traces_one
clear traces_two
clear traces_three
clear a
totalcol = zeros(6,1);
q = 1;
r = 1;
p = 1;
%goodpos = [0,10,3,4,5,7,8,9,10,12,14,18,21,22,24,26,27,28,29,31,32,33];
for k=1:length(ff)
    outfile = ff(k).name; %nms{k};
    %u = num2str(goodpos(k));
    %u2 = outfile(1:2);
    %if outfile(1) == u(1) ||  (u2(1)==u(1) && u2(2)==u(2) )
    
    load(outfile,'colonies');
    if ~exist('colonies','var');
        continue
    end
    %disp(['using outfile' num2str(goodpos(k)) ]);
    numcol = size(colonies,2); % how many colonies were grouped within the frame
    traces = cell(1,numcol);
    
    for j = 1:numcol
        if size(colonies(j).ncells_actual,1)>fr_stim                      % new segmentation
            colSZ =colonies(j).ncells_actual(fr_stim) ;                   % new segmentation
            %colSZ =colonies(j).numOfCells(timecolSZ-1) ;                   % for old mat files analysis
            if colSZ>0 && colSZ<6
                traces{j} = colonies(j).NucSmadRatio;                      %colonies(j).NucSmadRatio(:)
                % traces{j} = colonies(j).NucSmadRatioOld;                    % for old mat files analysis
                
                for h = 1:size(traces{j},2)
                    a = isfinite(traces{j}(:,h));                       %AN 
                    dat = zeros(size(traces{1},1),1);
                    dat(a == 1,1) = traces{j}(a==1,h);
                    
                    if length(nonzeros(dat))>trajmin       % FILTER OUT SHORT TRAJECTORIES(isnan(traces{j}(:,h))==0)
                        disp(['filter trajectories below' num2str(trajmin)]);
                        disp(['use' num2str(length(nonzeros(dat)))]);
                        totalcol(colSZ) = totalcol(colSZ)+1;
                        %traces{j}(traces{j}==0) = NaN;
                        figure(colSZ), plot(dat,'*','color',cmap(k,:));hold on% traces{j}(:,h)
                        ylim([0 2]);
                        xlim([5 sz]);
                        ylabel('mean Nuc/Cyto smad4 ');
                        xlabel('time, hours');
                        title(['All microColonies of size ' num2str(colSZ) ]);
                        
                        traces{j}(isnan(traces{j})==1) = 0;
                        d =  size(traces{j},2);
                        if colSZ == 1
                            traces_one{q}= dat;%dat
                        end
                        if colSZ == 2
                            traces_two{r}= dat;%traces{j}(:,h);
                        end
                        if colSZ >2% == 3 || colSZ == 4
                            traces_three{p}= dat;%traces{j}(:,h);
                        end
                        
                    end
                    q = q+1;
                    p = p+1;
                    r = r+1;
                end
            end
        end
    end                       % new segmentation
end
%end
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

%sz = 83;
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
s = size(traces_curr{k},2);% if the colony has less then 20 traces 
if s<=15
traces_one_new(:,q:q+s-1) = traces_curr{k}(:,:);%  smoothinditrace(traces_curr{k}(:,:),4) smoothtrace = smoothinditrace(vect,N);N = 4, smoothwindow
q = q+1;
end
end

fin_data = zeros(sz,2);
clear dat
for j =1:size(traces_one_new,1)
    for k=1:size(traces_one_new,2)
        dat =  nonzeros(traces_one_new(j,:));
        fin_data(j,1) = mean(dat(isfinite(dat)));
        fin_data(j,2) = std(dat(isfinite(dat)));
       
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
%load ('Feb_newSegm_NEWAnalysis30filt.mat','findatanew')%meanTraj_SDrmjunktraces
fr_stim = 16 ;        %22(jan8data) %38 %16(feb16 and july7 data)    12july26 data
delta_t = 15;%15       % 12 min              15min           17min   17 min
p = fr_stim*delta_t/60;
timecolSZ = fr_stim;%10
p2 = (timecolSZ)*delta_t/60;
%sz = 83;%81 

pfin = (sz)*delta_t/60;
cmap = colorcube;close all
cmap2 = hot;close all
 
 vect1 = (1:sz);
 startplot = 1;   %        START FROM TIME POITN5 SINCE THERE IS SOME JUNK IN THE FIRST TPT THAT MAKES THE SIGNAL JUMP UP
vect = (startplot:sz)*(delta_t)/60;
colN = {'b','g','r'};

for k=1:size(findatanew,2)
figure(10), errorbar(vect',findatanew{k}(startplot:end,1),findatanew{k}(startplot:end,2),'color',colN{k},'marker','*');hold on%findatanew{k}(:,2)
figure(10),title('All microColonies');
legend('1-cell colonies','2-cell colonies','3-cell colonies')
%text(40,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] ); 
xlim([0.5 pfin]);
ylim([0.5 2]);
ylabel('mean Nuc/Cyto smad4 ');
xlabel('time, hours');
figure(5), plot(vect',findatanew{k}(startplot:end,1),'color',colN{k},'marker','*'); hold on
figure(6), plot(vect',findatanew{k}(startplot:end,2),'color',colN{k},'marker','*');  hold on%smoothtrace smoothVAR(findatanew{k}(startplot:end,2),2)    power(findatanew{k}(startplot:end,2),2)
end

figure(5),title('Mean Trajectories');
legend('1-cell colonies','2-cell colonies','3- and 4-cell colonies')
text(40,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] ); 
xlim([0.5 pfin]);
ylim([0.5 2]);
ylabel('mean Nuc/Cyto smad4 ');
xlabel('time, hours');
figure(6),title('SD');
text(5,1,['colony size deremined at time  ' num2str(p2) ' hours'] );
legend('1-cell colonies','2-cell colonies','3- and 4-cell colonies');
xlim([0.5 pfin]);
ylim([0 0.7]);
ylabel('SD ');
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
