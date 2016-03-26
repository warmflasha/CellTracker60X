%%
%ilastikfile = ('/Users/warmflashlab/Desktop/IlastikMasks_headless_PluriW0/NucMaskPluri_tg56.h5');
ilastikfile = ('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/NucMsks3D/NucMasks3Djan8set1_z3.h5');
ilastikfile2 = ('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/CytoMasks3D/CytoMasks3Djan8set1_z3.h5');

nuc = h5read(ilastikfile,'/exported_data');
cyto = h5read(ilastikfile2,'/exported_data');
  k =20;% 41 15,14,16; ,42 % 43,44 , 52,53,54- no
    nuc = nuc(2,:,:,k);% for probabilities exported
    nuc = squeeze(nuc);
    mask1 = nuc;
    
    mask3 = imfill(mask1 > 0.9,'holes');
    
     cyto = cyto(2,:,:,k);% for probabilities exported
     cyto = squeeze(cyto);
     mask2 = cyto;
   % 
    figure(1), subplot(1,3,1),imshow(mask1);
    figure(1), subplot(1,3,2),imshow(mask2);
    figure(1), subplot(1,3,3),imshow(mask3);
   
    Lnuc = mask3;%im2bw(mask1,0.5);
    Lcyto = im2bw(mask2,0.9);
    figure(2),subplot(1,2,1), imshow(Lnuc);
    figure(2),subplot(1,2,2), imshow(Lcyto&~Lnuc);
    
   %%
   [MaskFin2] = Unmergetwonuclei(mask3);
   %[MaskFin3] = Unmergetwonuclei(MaskFin2);
   figure, imshow(MaskFin2);
   
%%
% obtain peaks for one frame 
zplane = [];
chan = [0 1];
paramfile = 'setUserParamLiveImagingAN';

ilastikDirec1  = ('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/Jan8IlastikMasks_newW0'); % 79=81 for position 26; 37-39 , pos 12
ilastikDirec2  = ('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/Jan8IlastikMasks_newW1');
imgDirec1 =('/Users/warmflashlab/Desktop/MaxProjections_Dif_Jan8run/W0');% already max projections
imgDirec2 =('/Users/warmflashlab/Desktop/MaxProjections_Dif_Jan8run/W1') ;% already max projections
dummy = ('/Users/warmflashlab/Desktop/Jan8setIlastikMasks_headless_DiffW1');
[nums, ilastikCytoAll]=folderFilesFromKeyword(dummy,'CytoMask');%make two ilastik directories
timegroups = 3;%
positions = length(nums)/timegroups;
 positions = 0:(positions-1);% vector with position numbers
       pos = 12;
    outfile = 'jan8set_test.mat';% basic name for all positions
peaks = nucCytoIlastik2peaksLoop(ilastikDirec1,ilastikDirec2,imgDirec1,imgDirec2,zplane,pos,chan,paramfile,outfile);% tsted
outfile = ([ num2str(pos) '_' num2str(outfile)]);

%%
 % run tracker on specific outfile
outfile = '8_3D_20hr_test_xyz.mat';
for k=1:length(peaks)
    if ~isempty(peaks{k})
       a = find(isnan(peaks{k}(:,1))) ;
       peaks{k}(a,:) = [];
    end
end
save(outfile,'imgfiles','imgfilescyto','peaks')
runTracker(outfile,'newTrackParamAN');
global userParam;
userParam.colonygrouping = 130;
cellsToDynColonies(outfile);

%%
% h5 track
ilstikfile3 = ('/Users/warmflashlab/Desktop/{test}_{track}.h5');
ilstikfile4 = ('/Users/warmflashlab/Desktop/{test}_{trackCyto}.h5');
k = 91;
trac_nuc = h5read(ilstikfile3,'/exported_data');
trac_nuc = squeeze(trac_nuc);
trac1 = trac_nuc(:,:,k);
bw2 = bwlabel(trac1);
figure,imshow(bw2,[]);

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
% run tracker and colony grouping an directory of mat files
%fr_stim = 16;
ff = dir('*jan8set_newparams*.mat');
for k=1:size(ff,1)
    outfile = ff(k).name;
%     outfile = '15_jan8set_optimization.mat';
runTracker(outfile,'newTrackParamAN');
global userParam;
%userParam.colonygrouping = 120;
% look at colonies around the stimulation frame (window of couple hours)?
cellsToDynColonies(outfile);
end
%%
% NEW: get new traces
% get time dependent data for a given .mat file ( plot each colony in
% different figure)
fr_stim = 22 ;%22 %38 %16
delta_t = 12; % 12 %15  in minutes
p = fr_stim*delta_t/60;
%global userParam;
colSZ = 3;
resptime =80;% 15 50 36 in frames ( converted to hours later)
jumptime = 5;% in frames
p2 = (resptime+jumptime)*delta_t/60;
coloniestoanalyze = 3;
%cmap = summer;
C = {'g','r','b','m','c'};
   % ff = dir('*12_jan8set_test*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs %Outfile
        %outfile = ff(k).name; %nms{k};
        %cellsToDynColonies(outfile);
        outfile = ('8_3D_20hr_test_xyz.mat');
        load(outfile,'colonies','peaks');
        tps = length(peaks);
        numcol = size(colonies,2); % how many colonies were grouped within the frame
        traces = cell(1,numcol);
        for j = 1:numcol
                 
                traces{j} = colonies(j).NucSmadRatio(:);
                traces{j}((traces{j} == 0)) = nan;
                figure(j+10), plot(traces{j},'-*','color',C{j});% cmap(j,:) 'r' traces
                ylim([0 2.5]);
                ylabel('mean Nuc/Cyto smad4  ');
                xlabel('frames');
            
        end
  
%%
% Get means over time but separately for different colony sizes
% new figure for each new colony size ( data drawn from multiple .mat
% files)

fr_stim = 22 ;        %22 %38 %16
delta_t = 12; 
p = fr_stim*delta_t/60;
timecolSZ = 20;
p2 = (timecolSZ)*delta_t/60;
cmap = colorcube;close all
cmap2 = hot;close all

%C = {'g','r','b','m'};
ff = dir('*_test*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs %Outfile
clear traces 
clear traces_one
clear traces_two
clear traces_three

q = 1;
r = 1;
p = 1;

for k=1:length(ff)
    outfile = ff(k).name; %nms{k};
    
    % outfile = ('12_3D_20hr_test_xyz.mat');
    load(outfile,'colonies','peaks');
    tps = length(peaks);
    numcol = size(colonies,2); % how many colonies were grouped within the frame
    traces = cell(1,numcol);
    
    for j = 1:numcol
        colSZ = colonies(j).numOfCells(timecolSZ); % colony size determined at the time of stimulation
        traces{j} = colonies(j).NucSmadRatio(:);
        traces{j}(traces{j}==0) = NaN;
        traces(cellfun(@isempty,traces)==1)=[];
        if colSZ>0 && colSZ<4
           
            figure(colSZ), plot(traces{j},'*','color',cmap(k,:));hold on% cmap(j,:) 'r' traces
            ylim([0 2.7]);
            ylabel('mean Nuc/Cyto smad4 ');
            xlabel('time, hours');
            title(['All microColonies of size ' num2str(colSZ) ]);
            text(60,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] );  
            
            traces{j}(isnan(traces{j})==1) = 0;
            d =  size(traces{j},2);
            if colSZ == 1
            traces_one{q}= traces{j};
            end
            if colSZ == 2
            traces_two{r}= traces{j};
            end
            if colSZ == 3
            traces_three{p}= traces{j};
            end
        end
         q = q+1; 
         p = p+1;
         r = r+1;
         
    end
end
traces_one = traces_one;
colSZ = 1;
traces_one(cellfun(@isempty,traces_one)==1)=[];

d = size(traces_one,2);
clear replace
clear sm

for k=1:d
a =size(traces_one{k},1);
if a < 99;
sm(k) = a;
end
end
a = find(sm>0);
sm1 = nonzeros(sm);

for jj=1:size(nonzeros(sm),1)
    
    replace{jj} = zeros(99,1);
    replace{jj}(1:sm1(jj),1) = traces_one{a(jj)}(:,1);
    traces_one{a(jj)} = replace{jj};
end
clear traces_one_new
for k=1:d
s = size(traces_one{k},2);
traces_one_new(:,q:q+s-1) = traces_one{k}(:,:);
q = q+1;
end

fin_data = zeros(99,2);

for j =1:size(traces_one_new,1)
        fin_data(j,1) = mean(nonzeros(traces_one_new(j,:)));
        fin_data(j,2) = std(nonzeros(traces_one_new(j,:)));
    
end
findatanew{colSZ} = fin_data;

vect1 = (1:99);
vect = (1:99)*12/60;



figure(1), hold on
plot(vect1',fin_data(:,1),'-*r','linewidth',3)

 
for k=1:size(findatanew,2)
figure(10+k), errorbar(vect',findatanew{k}(:,1),findatanew{k}(:,2),'color',cmap(k,:),'marker','*');hold on
figure(10+k),title('All microColonies');
legend('1-cell colonies','2-cell colonies','3-cell colonies')
text(40,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] ); 
xlim([0 20]);
ylim([0 2.5]);
ylabel('mean Nuc/Cyto smad4 ');
xlabel('time, hours');
figure(5), plot(vect',findatanew{k}(:,1),'color',cmap(k,:),'marker','*'); hold on
figure(6), plot(vect',power(findatanew{k}(:,2),2),'color',cmap(k,:),'marker','*');  hold on
end

figure(5),title('Mean Trajectories');
legend('1-cell colonies','2-cell colonies','3-cell colonies')
text(40,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] ); 
xlim([0 20]);
ylim([0.6 1.6]);
ylabel('mean Nuc/Cyto smad4 ');
xlabel('time, hours');
figure(6),title('Variance');
text(40,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] );
legend('1-cell colonies','2-cell colonies','3-cell colonies')
xlim([0 20]);
ylim([0 0.5]);
ylabel('Variance ');
xlabel('time, hours');
figure(10+k),title('All microColonies');
legend('1-cell colonies','2-cell colonies','3-cell colonies')
text(40,2.5,['colony size deremined at time  ' num2str(p2) ' hours'] ); 
xlim([0 20]);
ylim([0 2.5]);
ylabel('mean Nuc/Cyto smad4 ');
xlabel('time, hours');


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

      
      
      
