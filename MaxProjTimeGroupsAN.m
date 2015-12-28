

function [max_img] = MaxProjTimeGroupsAN(direc)
%if dt == 0
ff=readAndorDirectory(direc);
timegroups = size(ff.t,2);
tg = 1;
pos = 1;
nz = size(ff.z,2);
filename = cell(1,nz);

for j=1:nz
 
filename{j} = getAndorFileName(ff,pos,ff.t(tg),ff.z(j),ff.w(1)); % has to be channel 2 since all the masks should be applied to the gfp channel
end
% filename{1} = plane z0000; now if open it with bfopen will uncover all
% the timepoints within first time group taken at z0000

% for x = 1:nz
% imgs{x} = bfopen(filename{x}); open each timegroup and each z slice
% end

% then find the max image within those and save
% code below will do this, but need to replace filename{ii} with respective
% imgs{x}{1}{ii,1}
% saving is not coded too

%nframes = size(imgs{1},1);
%time = nframes;

for ii=1:nz
    %filename = getAndorFileName(files,pos,time,files.z(ii),chan);
    img_now = imread(filename{ii});
    if ii==1
        max_img=img_now;
    else
        max_img=max(img_now,max_img);
       end
end
% imgs = bfopen(filename);
% imgs_nuc = bfopen(filename2);
% nframes = size(imgs{1},1);
% time = nframes;
% end

%reader = bfGetReader(filename);