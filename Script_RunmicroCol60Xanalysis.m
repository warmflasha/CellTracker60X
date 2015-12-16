% script to run the live cell imaging microcolonies analysis( 60X data)
% need to be in the directory with ilastik-supplied masks 

%necessary parameters :
dir = '.';
direc =  ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/ANmicrocoloniesNov12(4)_20151116_102334AM');
%direc2 = ('/Volumes/Seagate Backup Plus Drive/RICE_Research_databackup/ANmicrocoloniesPluri42hr_20151129_95830 AM');
pos = 11;
zplane = 4; % for dataset from november 12 zplne = 4; for dataset from nov 25(pluri_42hr) zplane = 3;
dt = 0;
tg = [1 2 3 4]; % for dataset from november 12 4 time groups; for dataset from nov 25(pluri_42hr) 3 time groups;
fr_stim = 38;
isshift = 1;
col = [6 7];
N = 1;
delta_t = 5; % imaging was done every 5 minutes for the nov 12 dataset(diff condition); im pluri, this number is 11 mins


[peaks,dims,imgfilescyto,imgfiles] = RunFullTimeSerias60X_AN(dir,zplane,direc2,pos,dt,tg);

[colonies] = runMicroColonyGrouping(dir,fr_stim,pos,isshift);
% need to be with the folder with the final outfiles for each position (
% outfile_pos_tps.mat (tps = time points)
[datcell] = AnalyzeCellTraces_AN(dir,col,2,fr_stim,delta_t,flag);