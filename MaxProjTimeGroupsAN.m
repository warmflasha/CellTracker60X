
function  MaxProjTimeGroupsAN(direc,direc2,pos,tg,chan,Nchoose)
% direc = directory with raw images
% pos = the position to be processed
% tg = time group, vector
% chan = corresponds to ff.w(chan), e.g. ff.w(1) = nuclear channel
% direc2 = where to save the output projections
%

ff=readAndorDirectory(direc);
timegroups = size(ff.t,2);
nz = size(ff.z,2);
filename = cell(1,nz);
imgs = cell(1,nz);
for xx = 1:size(tg,2)
%xx = 0;
for j=1:nz
    if isempty(ff.t)
        filename{j} = getAndorFileName(ff,pos,[],ff.z(j),ff.w(chan));
    end
    if size(ff.t,2) > 0
    filename{j} = getAndorFileName(ff,pos,ff.t(1),ff.z(j),ff.w(chan));
    end
end

% filename{1} = plane z0000; now if open it with bfopen will uncover all
% the timepoints within first time group taken at z0000

for x = 1:nz
    imgs{x} = bfopen(filename{x});% open each timegroup and each z slice
    
end
% then find the max image within those and save
% code below will do this
nframes = size(imgs{1}{1},1);% how many time frames are there within each chunk same for all imgs{j}
if ~isempty(Nchoose)
nframes = Nchoose;% if onlt want selected number of time poins to be processed (e.g. there was a move after certain time point tand the data is not usable)
end
max_img = cell(1,nframes);
for k=1:nframes
    for ii=1:nz
        img_now = imgs{ii}{1}{k,1};
        if ii==1
            max_img=img_now;
        else
            max_img=max(img_now,max_img);
        end
    end
    %imwrite(max_img,[direc2 'W' num2str(ff.w(chan)) '_maxprojection_' num2str(pos) '_t' num2str(k) '.tif'],'Compression','none');%'writemode','append',
   if pos < 10 && isempty(ff.t)
    imwrite(max_img,[direc2 'Maxprojection_f000' num2str(pos) '_w000' num2str(ff.w(chan)) '.tif'],'writemode','append','Compression','none');%
   end
   if pos >= 10 && isempty(ff.t)
    imwrite(max_img,[direc2 'Maxprojection_f00' num2str(pos)  '_w000' num2str(ff.w(chan)) '.tif'],'writemode','append','Compression','none');%
   end
   
   if pos < 10 && ~isempty(ff.t)
    imwrite(max_img,[direc2 'Maxprojection_f000' num2str(pos) '_t000' num2str(tg(xx)) '_w000' num2str(ff.w(chan)) '.tif'],'writemode','append','Compression','none');%
   end
   if pos >= 10 && ~isempty(ff.t)
    imwrite(max_img,[direc2 'Maxprojection_f00' num2str(pos) '_t000' num2str(tg(xx)) '_w000' num2str(ff.w(chan)) '.tif'],'writemode','append','Compression','none');%
   end
   
   
   % the line above saves all the projections into a miltitif (suitable
    % for ilastik training )
    % Maxprojection_f0000_t0000_w000
end
end
end

