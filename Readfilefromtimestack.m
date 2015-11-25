
direc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/ANmicrocoloniesNov12(4)_20151116_102334AM');
ff=readAndorDirectory(direc);
size(ff.t,2)%number of groups of time points
%filename = getAndorFileName(files,pos,time,z,w);
filename = getAndorFileName(ff,pos,ff.t(1),ff.z(zplane),ff.w(2));
% filename is just one string
imgs = bfopen(filename);
nframes = size(imgs{1},1);% number of images within this time chunk
for k=1:nframes
fnm = imgs{1}{k,1};% fnm s an actual image , specified by pos, zplane, ff.w(2) and the real frame number is k
end
%imgs{1}(:,2): the first cell of imgs has all the filenames 
% so the 


% img80 = imgs{1}{80,1};
% open splitHyperStackByTime
% reader = bfGetReader(filename);
% reader.getSizeC
% reader.getSizeT
% reader.getSizeZ
% iPlane = reader.getIndex(0,0,79);
% img80b=bfGetPlane(reader,iPlane);
% figure; imshow(img80b,[]);

% other useful functions:
% saveOneZstack; splitHyperStackByTime