% only the nuclear segmentation of 40X data of the live cell chips that
% were fixed and stained and then images again (back to coordinates)

%ilastikdircyto = ('/Users/warmflashlab/Desktop/July26Tiling2_CytoMasks');
ilastikdirnuc =('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-26-Tiling2_24hrtotalBMP410ngml/nucmasksFixedChip');
% imagedir1 =('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-26-Tiling2_24hrtotalBMP410ngml/uColTiling2_20hrinBMP4_20160726_24717PM/nucrawdata_z2z3');
% imagedir2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-26-Tiling2_24hrtotalBMP410ngml/uColTiling2_20hrinBMP4_20160726_24717PM/cytorawdata_z2z3');
imagedir=('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-26-Tiling2_24hrtotalBMP410ngml/backtoCdx2pos_tiling2_20160730_43330 PM');

timegroup = [];%1

%chan = [1 2];
chanal = [1 2 3];
paramfile = 'setUserParamLiveImagingAN_40X';       %'setUserParamLiveImagingAN_40X   setUserParamLiveImagingAN';
paramfile3D = 'setUserParam3DsegmentationAN_40X';  %'setUserParam3DsegmentationAN    setUserParam3DsegmentationAN_40X'

outfile = 'testBacktocoord.mat';%3Dsegm_febdata

  positions = [11]; %july 26 tiling 2 position (40X)
  %positions = [0 2 3 4 5 8 11 15 16 17 18 20 21 22 25 26 29 30 32 33 10]; %  february 3 dataset
    %positions = [5 7 8 9 11 13 15 16 17 18 20 21 23 25 27 28 29 33 34 35 36 37 38 39 40];
 % July 7 dataset 40 poszitions therelastpart

pl =1;%2,3 5
strnuc = 'nucmask';   %NucMasks3Djan8set   %nucmask3DFebset %nuc_mask_cytomask3DJuly7set
k = 1;
%for k = 1:size(positions,2)
    pos = positions(k);
    rundataset3Donlynucsegm(ilastikdirnuc,imagedir,pos,paramfile,timegroup,outfile,paramfile3D,pl,strnuc,chanal(1))
%end
matfile = '12_testBacktocoord.mat';
load(matfile);
Lnuc = imgfiles(1).NucMask;
figure(10), imshow(Lnuc);
figure(10),hold on,plot(peaks{1}(:,1),peaks{1}(:,2),'r*');
hold on,text(peaks{1}(:,1)+10,peaks{1}(:,2)+10,num2str(peaks{1}(:,6)./peaks{1}(:,5)),'color','g'); hold on
hold on,text(peaks{1}(:,1)+20,peaks{1}(:,2)+20,num2str(peaks{1}(:,7)./peaks{1}(:,5)),'color','m'); hold on
