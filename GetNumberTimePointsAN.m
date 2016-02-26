function [nT] = GetNumberTimePointsAN(imagedir,pos,timegroup)%

% reading raw data
ff=readAndorDirectory(imagedir);
nz = 1;
filename = cell(1,nz);
imgs = cell(1,nz);

%dirinfo1(start1).name
for j=1
    filename{j} = getAndorFileName(ff,pos,ff.t(timegroup),ff.z(j),ff.w(1));
end
for m = 1
    imgs{m} = bfopen(filename{m});  %
    img_now = imgs{m}{1}{1,1};      % get the size infor from the first time point
    nT = size(imgs{m}{1},1);
end







end