% info on where the images are and whidh position / timepoint to run

ilastikdircyto = ('/Users/warmflashlab/Desktop/Feb2016ilastik_CytoMasks');
ilastikdirnuc = ('/Users/warmflashlab/Desktop/Feb2016ilastik_NucMasks');
imagedir1 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/03-02-2016-uCol_diff_AF(83tptsusable)/Nuc_raw_data');%rawimages_nuc
imagedir2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/03-02-2016-uCol_diff_AF(83tptsusable)/Cyto_raw_data');%rawimages_cyto

%imagedir =('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/testPos1_rawImgs');
%pos = 12;
%tpt =1;
timegroup = 1;
chan = [1 2];
chanal = 1;
paramfile = 'setUserParamLiveImagingAN';
paramfile3D = 'setUserParam3DsegmentationAN';

outfile = 'frame_test.mat';


% janyary8 dataset : positions = [1 2 5 7 8 12 14 15 18 19 20 24 25 27 29 30]; %positions that were processed
  positions = [0 1 2 3 4 5 8 9 11 12 13 14 15 16 17 18 20 21 22 25 26 27 32 33]; % 29 30 february 3 dataset

pl = 3;
strcyto = 'cytomask3DFebset'; %CytoMasks3Djan8set
strnuc = 'nucmask3DFebset';   %NucMasks3Djan8set
%for k=1:length(positions)
k = 13;
    pos = positions(k);
rundataset3D(ilastikdirnuc,ilastikdircyto,imagedir1,imagedir2,pos,paramfile,timegroup,outfile,paramfile3D,pl,strnuc,strcyto,chanal);
%end