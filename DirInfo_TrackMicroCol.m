% info on where the images are and whidh position / timepoint to run
 
 
% tiling data set July 7

% ilastikdircyto = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/trainingset/cyto_mask_15_z1.h5');%/cyto_mask_z1.h5'
% ilastikdirnuc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/trainingset/nuc_mask_15_z1.h5');%/nuc_mask_z1.h5'
ilastikdircyto = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/July7TilingLCellilastik_CytoMasks');%/cyto_mask_z1.h5'
ilastikdirnuc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/July7TilingLCellilastik_NucMasks');%/nuc_mask_z1.h5'
imagedir1 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/nuc_raw');%rawimages_nuc%
imagedir2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/cyto_raw');%rawimages_cyto

% Feb2016 data
 
% ilastikdircyto = ('/Users/warmflashlab/Desktop/Feb2016ilastik_CytoMasks');
% ilastikdirnuc = ('/Users/warmflashlab/Desktop/Feb2016ilastik_NucMasks');%/nuc_mask_z1.h5'
% imagedir1 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/03-02-2016-uCol_diff_AF(83tptsusable)/Nuc_raw_data');
% imagedir2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/03-02-2016-uCol_diff_AF(83tptsusable)/Cyto_raw_data');

% location of stuff for the January8 data 3D
% ilastikdircyto = ('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/CytoMasks3D');
% ilastikdirnuc =('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/NucMsks3D');
% imagedir1 =('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/3DanalysisRAWimg_W0');
% imagedir2 = ('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/3DanalysisRAWimg_W1');

% location of stuff for the January8 data 2D ( September rerun )

%  imagedir1 =('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/2D_ToProcess_Jan8Data/W0');
%  imagedir2 =('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/2D_ToProcess_Jan8Data/W1');
%  ilastikdirnuc = ('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/2016-09-08-nucmask');
%  ilastikdircyto =('/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/2016-09-08-cytomask');

%imagedir =('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/testPos1_rawImgs');
%pos = 12;
%tpt =1;
timegroup = [];%1 []

chan = [1 2];
chanal = 1;
paramfile = 'setUserParamLiveImagingAN_40X';%'setUserParamLiveImagingAN_40X setUserParamLiveImagingAN';
paramfile3D = 'setUserParam3DsegmentationAN_40X ';%'setUserParam3DsegmentationAN    setUserParam3DsegmentationAN_40X'

outfile = '40X_imprBGandSegm.mat';%3Dsegm_febdata
%positions = [3 11 13 15 17 19 20 24 26 31];  jan8 2D
  %positions = [12 14 15 18 19 20 24 25 27 29 30]; %janyary8 dataset positions that were processed
  %positions = [0 2 3 4 5 8 11 15 16 17 18 20 21 22 25 26 29 30 32 33 10]; %  february 3 dataset
    %positions = [5 7 8 9 11 13 15 16 17 18 20 21 23 25 27 28 29 33 34 35 36 37 38 39 40];
 % July 7 dataset 40 positions therelastpart
positions = [33 35 36 39 40];% 33: 2z, 35: 2z 36: 2z 39 z1 40: 1z
pl =2;%2,3 5

strcyto = 'cytomask3DJuly7set'; % jan8datacytomask CytoMasks3Djan8set  %cytomask3DFebset % nucmask3DJuly7setcyto_mask_
strnuc = 'nucmask3DJuly7set';   % jan8datanucmask  NucMasks3Djan8set   %nucmask3DFebset %nuc_mask_  cytomask3DJuly7set
%k = 5;
for k = 1:size(positions,2)
    pos = positions(k);
    if pos == 39 || pos == 40 
        pl = 1;
    end
    rundataset3D(ilastikdirnuc,ilastikdircyto,imagedir1,imagedir2,pos,paramfile,timegroup,outfile,paramfile3D,pl,strnuc,strcyto,chanal);
end


