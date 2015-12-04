function [peaks,dims,imgfilescyto,imgfiles] = RunTimeSeries60XuColoniesAN(ilastikfile,ilastikfilecyto,pos,zplane,direc,dt,tg)

info = h5info(ilastikfile);
info.Datasets;   

data = h5read(ilastikfile,'/exported_data');
data = squeeze(data);

time = size(data,4);   

peaks = cell(1,time); 
peakscyto = cell(1,time);

imgfiles = struct([]);
imgfilescyto = struct([]);
%lblobjects = struct([]);
if dt == 0
ff=readAndorDirectory(direc);
timegroups = size(ff.t,2);


filename = getAndorFileName(ff,pos,ff.t(tg),ff.z(zplane),ff.w(2)); % has to be channel 2 since all the masks should be applied to the gfp channel
filename2 = getAndorFileName(ff,pos,ff.t(tg),ff.z(zplane),ff.w(1)); % to get info from the nuc channel
imgs = bfopen(filename);
imgs_nuc = bfopen(filename2);

nframes = size(imgs{1},1);
time = nframes;
end
for k=1:time
    
[datacell,Lnuc,Lcytofin] = IlastikplusWatershed_AW(ilastikfile,ilastikfilecyto,pos,zplane,direc,k,dt,tg,imgs,imgs_nuc);

peaks{k} = datacell;
if isempty(datacell)
continue
end

imgfiles(k).compressNucMask = compressBinaryImg(Lnuc, size(Lnuc) );
imgfilescyto(k).compressNucMask = compressBinaryImg(Lcytofin, size(Lcytofin) );
%lblobjects(k).lblobj = a;
end
dims = size(peaks);
 %save('Outfile_test2','peaks','dims','imgfiles','imgfilescyto','colonies');
 end

