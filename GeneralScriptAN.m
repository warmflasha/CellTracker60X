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
 
outfile = '14_3D_20hr_test_xyz.mat';
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

% trac_cyto = h5read(ilstikfile4,'/exported_data');
% trac_cyto = squeeze(trac_cyto);
% trac2 = trac_cyto(:,:,k);
% bw3 = bwlabel(trac2);
% figure,imshow(bw3,[]);


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


% badinds = [bwn2.PixelIdxList] < userParam.areanuclow; %this area should also be a parameter
% badinds2 = [statsnucw0.Area] < userParam.areanuclow;
% statsnuc(badinds) = [];
% statsnucw0(badinds2) = [];
%%
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
fr_stim = 22 ;%22 %38 %16
delta_t = 12; % 12 %15  in minutes
p = fr_stim*delta_t/60;
%global userParam;
colSZ = 3;
resptime =80;% 15 50 36 in frames ( converted to hours later)
jumptime = 15;% in frames
p2 = (resptime+jumptime)*delta_t/60;
coloniestoanalyze = 3;
cmap = summer;
flag = 1;

C = {'b','r','g','m'};
   % ff = dir('*12_jan8set_test*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs %Outfile
  
    % for k=1:size(ff,1);
    %  k=1;
        
        %outfile = ff(k).name; %nms{k};
        %cellsToDynColonies(outfile);
        outfile = ('2_3D_20hr_test_xyz.mat');
        load(outfile,'colonies','peaks');
        tps = length(peaks);
        numcol = size(colonies,2); % how many colonies were grouped within the frame
        traces = cell(1,numcol);
        
        for j = 1:numcol
            
%             if colSZ == colonies(j).numOfCells(fr_stim);
%                 frames(k,1) =  k;% which outfile contains the traces for the colony size colSZ
%             end
            if flag ==1
                traces{j} = colonies(j).NucSmadRatio(:);
                traces{j}((traces{j} == 0)) = nan;
                figure(j+1), plot(traces{j},'-*','color',C{j});% cmap(j,:) 'r' traces
                ylim([0 2.5]);
                ylabel('mean Nuc/Cyto smad4  ');
                xlabel('frames');
            end
        end
    %end
%%
% NEW: get the before and after  plots
clear all;
fr_stim = 22;% 16 22 %38
fldat = [2 3];
delta_t = 12; % 12% in minutes
p = fr_stim*delta_t/60;
global userParam;
%userParam.colonygrouping = 120;
flag = 1;
resptime =75;% 15 50 36 in frames ( converted to hours later)
range = [26 20];
jumptime = 5;% in frames 5
p2 = (resptime+jumptime)*delta_t/60;
coloniestoanalyze = 3;
cmap = parula;
C = {'b','g','r','k ','m','y'};
for jj = 1:5
q(jj) = 1;
w(jj) = 1;
s(jj) = 1; 
end
%for jj = 1:coloniestoanalyze
    
    %colSZ = jj;
    ff = dir('*jan8set_test*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs %Outfile  %dec31_set_Diff
    for k=1:size(ff,1)
        outfile = ff(k).name ;
        %cellsToDynColonies(outfile);
        load(outfile,'colonies','peaks');
        tps = length(peaks);
        
        numcol = size(colonies,2); % how many colonies were grouped within the frame
        for j = 1:numcol
            if ~isempty(fr_stim) % && k~=22;
                
                %   if colSZ == colonies(j).numOfCells(fr_stim); % how many cells within colony at time fr_stim
                
                colSZ = colonies(j).numOfCells(fr_stim); % how many cells within colony at time fr_stim
                %colSZaft = colonies(j).numOfCells(fr_stim+jumptime);
                
                if colSZ >0 %&& (colSZ==colSZaft)
                   
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
                    disp([num2str(q(colSZ)) '   ' num2str(currlengthbf) '   ' num2str(bf_fin{colSZ}')]);
                    aft_fin{colSZ}((s(colSZ):(s(colSZ)+currlengthaft)-1),1) = nonzeros(aft);
                    %disp([num2str(s) '   ' num2str(currlengthaft) '   ' num2str(aft_fin{colSZ}')]);
                    window_fin{colSZ }((w(colSZ):(w(colSZ)+currlengthgap)-1),1) = nonzeros(window);
                    %disp([num2str(w) '   ' num2str(currlengthaft) '   ' num2str(window_fin{colSZ }')]);
                                        
                    q(colSZ) = q(colSZ)+currlengthbf;
                    w(colSZ) = w(colSZ)+currlengthgap;
                    s(colSZ) = s(colSZ)+currlengthaft;
                    
                    
                    hold on,figure(4), plot(bf,aft,'*','color',C{colSZ},'markersize',15);%C{colSZ}
                    title(['mean Nuc/Cyto smad4 ' num2str(p2) 'hrs after bmp4']);
                    ylim([0 2]);
                    hold on,figure(5),subplot(1,2,1), plot(colSZ,bf,'*','color',C{colSZ},'markersize',15);
                    ylim([0 2]);
                    hold on,figure(5),subplot(1,2,2), plot(colSZ,aft,'*','color',C{colSZ},'markersize',15);
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
            
            if isempty(fr_stim)
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
                    disp([num2str(q) '   ' num2str(currlengthbf) '   ' num2str(bf_fin{colSZ }')]);
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
coloniestoanalyze = 4;
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

      
      
      
