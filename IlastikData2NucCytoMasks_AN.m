function [ratios, cellsinframe] = IlastikData2NucCytoMasks_AN(ilastikfile,ilastikfilecyto,pos,zplane,direc,flag)
% direc shoule be the one from where you get the images to apply the masks
% ilastikfile = h5 file containing the masks of all nuclei( expoted from
% the ilastik )
% ilastikfilecyto=h5 file containing the masks of all cytoplasms( expoted from
% the ilastik )
% pos= the position for which the segmentation was done
% zplane= the relevant plane to use out of z stack
% flag= if 1, then will show the objects segmented in color  
% 
% 

%direc = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/SingleCellSignalingAN_20150805_123245 PM';
%I2 = imread('SingleCellSignalingAN_t0000_f0019_z0003_w0001.tif');% gfp channel (gfp-smad4 cells)
% I = imread('SingleCellSignalingAN_t0000_f0019_z0003_w0000.tif');% nuc chan
ff=readAndorDirectory(direc);
chan = ff.w;

areanuclow = 700;
areanuchi = 9000;
areacytolow = 2000;
areacytohi = 20000;
global userParam;
userParam.gaussRadius = 10;
userParam.gaussSigma = 3;
userParam.small_rad = 3;
userParam.presubNucBackground = 1;
userParam.backdiskrad = 300;


info = h5info(ilastikfile);
infocyto = h5info(ilastikfilecyto);
info.Datasets;                                         % check these to make sure the dataset name is correct for the 'h5read' function
infocyto.Datasets;

data = h5read(ilastikfile,'/exported_data');
datacyto = h5read(ilastikfilecyto,'/exported_datacyto');

data = squeeze(data);
datacyto = squeeze(datacyto);
time = size(data,3);                                   % how many time frames were expoted ( depends on the dataset, set in ilastik)
Lnuc = {};
Lcyto = {};
for k=3:time
    
Lnuc{k} = data(:,:,k) >1;                             % need to leave only the nuclei masks and make the image binary
Lnuc{k} =  bwareafilt(Lnuc{k}',[areanuclow areanuchi]);%need to transpose the image matrix since for some reason the masks come out rotated from ilastik

Lcyto{k} = datacyto(:,:,k) >1; 
Lcyto{k} =  bwareafilt(Lcyto{k}',[areacytolow areacytohi]);%transposed for the same reason
Lcytofin{k} = Lcyto{k} & ~ Lnuc{k};                        % cyto masks initially include both nuclei+cyto, so need to eliminate nuc
Lcytofin{k} =  bwareafilt(Lcytofin{k},[areacytolow areacytohi]);%

Lnuc{k} = bwlabel(Lnuc{k},8);                             % convert back to the label matrix in order to match the correct nuc to the correct cyto
Lcytofin{k} = bwlabel(Lcytofin{k},8); 

% subtract the background from the original grayscale image
% preprocess
filename = getAndorFileName(ff,pos,ff.t(k),ff.z(zplane),chan(2)); % has to be channel 2 since all the masks should be applied to the gfp channel
I2 = imread(filename);
% I2 = imopen(I2,strel('disk',userParam.small_rad)); % remove small bright stuff
% I2 = smoothImage(I2,userParam.gaussRadius,userParam.gaussSigma); %smooth
% I2 =presubBackground_self(I2);


cc_nuc{k} = bwconncomp(Lnuc{k},8);
statsnuc{k} = regionprops(cc_nuc{k},I2,'Area','Centroid','PixelIdxList','MeanIntensity'); % statst for the nuc 'area' in cyto channel
%statsnuc_nuc = regionprops(cc_nuc,I,'Area','Centroid','PixelIdxList','MeanIntensity'); % stats for the nuclear channel
aa = [statsnuc{k}.Centroid];
Inucchan  = [statsnuc{k}.MeanIntensity];
anuc{k} = [statsnuc{k}.Area]; % nuclear area in cyto channel
xnuc = aa(1:2:end);
ynuc = aa(2:2:end);

%get the cytoplasmic mean intensity for each labeled object
cc_cyto{k} = bwconncomp(Lcytofin{k});
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
    Lrgbcyto = label2rgb(Lcytofin{k}, 'jet', 'k', 'shuffle');
    subplot(1,2,2),imshow(I2,[]);hold on
    h = imshow(Lrgbcyto);
    h.AlphaData = 0.3;
end

% plot the results

ratios{k}(:,1) = nucmeanInt{k}./cytomeanInt{k};
ratios{k}(:,2) = k;
cellsinframe{k} = size(statsnuc{k},1);% need to put this as a label next to each datapoint
meanRa{k} = mean(ratios{k}(:,1));

%vect = (1:time);
figure(4),plot(k,meanRa{k},'r--*','Markersize',15);hold on
%ylim([0 1.7]);
xlim([0 time+1]);
xlabel('time, frames');
ylabel('nuc/cyto mean for cells in the frame');
title('GFPsmad4RFPh2b cells, 10ng/ml bmp4 added after frame 16');


 %twocellcol(k) = mean(ratios{k}(2:3,1));
 onecellcol1{k} = ratios{k}(1,1);
 onecellcol2{k} = ratios{k}(2,1); 

hold on
plot(k,onecellcol1{k},'b--*');
plot(k,onecellcol2{k}, 'm--*');
ylim([0.6 1.2]);
legend('all two cells','1-cell colony(1)','1-cell colony(2)') ;

%  save('Frame0_nuc2cyto_ratios','ratios','meanRa','onecellcol1','onecellcol2');
%  savefig('Cellraces.fig');

end



end


%--------------chunk of code from Sapna
% immask = h5read(strcat(segfilepath ,dirinfo(file).name), '/exported_data');
% %nchannel = size(immask,1);
% %ntime = size(immask,4);
% %Case 1: nchannel = 1
% 
%  ntime = 2;
%  time = 1:ntime
%     mask = immask(segchannel,:,:,time);
%     mask_s = squeeze(mask);
%     thresh = graythresh(mask_s); 
%     mask_b = im2bw(mask_s, thresh); % cells:1, background:2
%     %mask_bc = imcomplement(mask_b);% we need background as '0'
%     mask_bct = mask_b'; % ilastik default settings have axes exchanged somehow.
    
%--------------chunk of code from Sapna