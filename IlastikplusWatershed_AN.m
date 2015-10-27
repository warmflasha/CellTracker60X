function [statsnuc,statscyto, cellsinframe,Lnuc,Lcyto] = IlastikplusWatershed_AN(ilastikfile,pos,zplane,direc,flag)

ff=readAndorDirectory(direc);
chan = ff.w;

global userParam;
userParam.gaussRadius = 10;
userParam.gaussSigma = 3;
userParam.small_rad = 3;
userParam.presubNucBackground = 1;
userParam.backdiskrad = 300;
areanuclow = 700;
areanuchi = 9000;

info = h5info(ilastikfile);
info.Datasets;                                         % check these to make sure the dataset name is correct for the 'h5read' function

data = h5read(ilastikfile,'/exported_data');
data = squeeze(data);

time = size(data,3);                                   % how many time frames were expoted ( depends on the dataset, set in ilastik)
Lnuc = {};
Lcyto = {};
for k=3:time
    
Lnuc{k} = data(:,:,k) >1;  % need to leave only the nuclei masks and make the image binary
Lnuc{k} =  bwareafilt(Lnuc{k}',[areanuclow areanuchi]);

filename = getAndorFileName(ff,pos,ff.t(k),ff.z(zplane),chan(2)); % has to be channel 2 since all the masks should be applied to the gfp channel
I2 = imread(filename);

Lcyto{k} = WatershedsegmCytoplasm_AW(Lnuc{k},I2,flag);% get the cyto mask using watershed

% at this point should have an array of nuc and cyto masks(from ilastik
% and watershed respectively)
%I2 = imread(filename);
I2 = imopen(I2,strel('disk',userParam.small_rad)); % remove small bright stuff
I2 = smoothImage(I2,userParam.gaussRadius,userParam.gaussSigma); %smooth
I2 =presubBackground_self(I2);

%get the NUCLEAR mean intensity for each labeled object
cc_nuc{k} = bwconncomp(Lnuc{k},8);
statsnuc{k} = regionprops(cc_nuc{k},I2,'Area','Centroid','PixelIdxList','MeanIntensity');
aa = [statsnuc{k}.Centroid];
Inucchan  = [statsnuc{k}.MeanIntensity];
anuc{k} = [statsnuc{k}.Area]; % nuclear area in cyto channel
xnuc = aa(1:2:end);
ynuc = aa(2:2:end);

%get the cytoplasmic mean intensity for each labeled object
cc_cyto{k} = bwconncomp(Lcyto{k});
statscyto{k} = regionprops(cc_cyto{k},I2,'Area','Centroid','PixelIdxList','MeanIntensity');
acyto{k} = [statscyto{k}.Area];
xcyto = aa(1:2:end);
ycyto = aa(2:2:end);

nucmeanInt{k}  = [statsnuc{k}.MeanIntensity]; %k is the index over images!!!!! statsnuc is a cell array
cytomeanInt{k}  = [statscyto{k}.MeanIntensity];


% To show which objects are labeled in the final watershed segmentation in nuc
% lear and cyto channels; also show the final foreground markers
if flag == 1
    Lrgb = label2rgb(Lnuc{k}, 'jet', 'k', 'shuffle');
    figure,subplot(1,2,1),imshow(I2,[]);hold on
    h = imshow(Lrgb);
    h.AlphaData = 0.3;
    Lrgbcyto = label2rgb(Lcyto{k}, 'jet', 'k', 'shuffle');
    subplot(1,2,2),imshow(I2,[]);hold on
    h = imshow(Lrgbcyto);
    h.AlphaData = 0.3;
end

%plot the results

% ratios{k}(:,1) = nucmeanInt{k}./cytomeanInt{k};
% ratios{k}(:,2) = k;
% ratios{k}(:,3) = xnuc';
% ratios{k}(:,4) = ynuc';
% cellsinframe{k} = size(statsnuc{k},1);% need to put this as a label next to each datapoint
% meanRa{k} = mean(ratios{k}(:,1));
% 
% %vect = (1:time);
% figure(4),plot(k,meanRa{k},'r--*','Markersize',15);hold on
% ylim([0.4 1.4]);
% xlim([0 (time+1)]);
% xlabel('time, frames');
% ylabel('nuc/cyto mean for cells in the frame');
% title('GFPsmad4RFPh2b cells, 10ng/ml bmp4 added after frame 16, ~ 9 hours imaging time');
% legend('Mean over all cells in frame') ;
% 
% 
% %  twocellcol1{k} = mean(ratios{k}(2:3,1));
% %  twocellcol2{k} = mean(ratios{k}(2:3,1));
%  
%   onecellcol1{k} = ratios{k}(1,1);
%   onecellcol2{k} = ratios{k}(2,1); % 
%   %onecellcol3{k} = ratios{k}(3,1); % 
% 
% %figure(5),plot(k,twocellcol{k},'b--*');hold on
% figure(5),plot(k,meanRa{k},'r--*','Markersize',15);hold on
% %plot(k,onecellcol3{k},'b--*');
% plot(k,onecellcol1{k}, 'm--*');
% plot(k,onecellcol2{k}, 'g--*');
% ylim([0.2 1.5]);
% xlabel('time, frames');
% ylabel('nuc/cyto mean for cells in the frame');
% title('GFPsmad4RFPh2b cells, 10ng/ml bmp4 added after frame 16,~ 9 hours imaging time');
% legend('mean over all cells in frame','one-cell colony(1)','one-cell colony(1)') ;



end

%  figure(4), savefig('MeanDynamicsframe10.fig');
%  figure(5),savefig('CellTracesframe10.fig');
% 
  save('Frame13_Analysis');

end



