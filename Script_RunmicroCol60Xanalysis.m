% script to run the live cell imaging microcolonies analysis( 60X data)
% need to be in the directory with ilastik-supplied masks 

%necessary parameters :

dir = ('~/Dropbox/Warmflash Lab/Anastasiia/AN_liveCellOutfiles/TestTraces');
%direc =  ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/ANmicrocoloniesNov12(4)_20151116_102334AM');
direc = ('/Volumes/data/Anastasia/LiveCellImagingGFPs4RFPh2b/ANmicrocoloniesNov12(4)_20151116_102334 AM');
%direc = ('E:\RICE_Research_databackup\ANmicrocoloniesPluri42hr_20151129_95830 AM');

pos = 12;
%zplane = 2; % for dataset from november 12 zplne = 4; for dataset from nov 25(pluri_42hr) zplane = 2;
dt = 0;
tg = [1 2 3 4]; % for dataset from november 12 4 time groups; for dataset from nov 25(pluri_42hr) 3 time groups;
fr_stim = [];
isshift = [];
col = [2 3];
N = 1;
delta_t = 5; % imaging was done every 5 minutes for the nov 12 dataset(diff condition); im pluri, this number is 11 mins
%%
% after running, this function will save the outfile_pos_tps.mat (tps =
% all time points for this position)
[peaks,dims,imgfilescyto,imgfiles] = RunFullTimeSerias60X_AN(dir,zplane,direc,pos,dt,tg);
%%
[colonies] = runMicroColonyGrouping(dir,fr_stim,pos,isshift);
% need to be with the folder with the final outfiles for each position (
% outfile_pos_tps.mat (tps = time points)
[datcell] = AnalyzeCellTraces_AN(dir,col,2,fr_stim,delta_t,flag,pos);