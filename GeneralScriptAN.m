%%
ilastikfile = ('/Users/warmflashlab/Desktop/IlastikMasks_headless_PluriW0/NucMaskPluri_tg56.h5');
nuc = h5read(ilastikfile,'/exported_data');
  k = 12;
    nuc = nuc(2,:,:,k);% for probabilities exported
    nuc = squeeze(nuc);
    mask1 = nuc;
    
    imshow(mask1);
    Lnuc = im2bw(mask1,0.5);
    figure, imshow(Lnuc);


%%
ff = dir('*jan8set*.mat');
k = 1;
outfile = ff(k).name;
load(outfile);

%%
N  =5;
n = uncompressBinaryImg(imgfiles(N).compressNucMask);
nc = uncompressBinaryImg(imgfilescyto(N).compressNucMask);
%close all
figure,subplot(1,2,1),imshow(n);

hold on
subplot(1,2,2),imshow(nc);
%%
fr_stim = 38;
ff = dir('**.mat');
for k=1:size(ff,1)
    outfile = ff(k).name;
runTracker(outfile,'newTrackParamAN');
global userParam;
userParam.colonygrouping = 120;
% look at colonies around the stimulation frame (window of couple hours)?
cellsToDynColonies(outfile);
end

%%
ff = dir('*10ngmlDifferentiated_22hrs*.mat');%jan8set % Pluri_42hrs % 22hr set: k = 3,4,6,16,22 traces k = 10 bad
k =22;
outfile = ff(k).name;
load(outfile);

fr_stim = 38;%22
fldat = [2 3];
delta_t = 5; % 12% in minutes
p = fr_stim*delta_t/60;
flag = 1;
resptime = 36;% 36in frames ( converted to hours later)
GetDynamicColonyTraces(outfile,fr_stim,fldat,delta_t);

%%
for j=1:size(ncells,1)
plot(ncells{j},'-*');hold on
ncells{j}(10)
end

%%
% plot the stats (before, after, amplitude...)
fr_stim = 22;%22 %38
fldat = [2 3];
delta_t = 12; % 12% in minutes
p = fr_stim*delta_t/60;
flag = 1;
resptime =150;% 15 50 36 in frames ( converted to hours later)
coloniestoanalyze = 3; % max size of the colonies that want to look at
for jj = 1:coloniestoanalyze
    
    colSZ = jj;
    
    %GetDynamicColonyTraces(outfile,fr_stim,fldat,delta_t);
    %ff = dir('*10ngmlDifferentiated*.mat');
    ff = dir('*jan8set*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs
    for k=1:size(ff,1)
        outfile = ff(k).name; %nms{k};
        load(ff(k).name);
        datafin = GetDynamicColonyStats(outfile,fr_stim,delta_t,flag,colSZ,resptime,coloniestoanalyze);
        
    end
end
if ~isempty(fr_stim)
figure(11)
xx = 0:0.1:2.5;
plot(xx,xx,'-k','linewidth',2);

end
if isempty(fr_stim) 
 figure(11) 
 xx = 0:1:coloniestoanalyze;
 yy = ones(1,coloniestoanalyze+1);
 plot(xx,yy,'--');
 figure(12)
 plot(xx,yy,'--');
end
%%
clear colonies
clear cells
clear ncells

% for k=1:length(peaks)
% a(k)= (~isempty(peaks{k}));
% end
% b = find(a == 1);
%  ncells = size(peaks{b(1)},1) ;


