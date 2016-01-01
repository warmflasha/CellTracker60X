function  GetZProjections_allPositionsAN(direc,direc2,Npos,tg,chan)

% direc = directory with raw images
% Npos = number of positions, need to know from experiment (or get from the
% directory , line 15, need to debug
% tg = time group, vector
% chan = corresponds to ff.w(chan), e.g. ff.w(1) = nuclear channel
% direc2 = where to save the output projections
% tg = a vector with the number of separate time groups, need to know it
% from experiment

%get the number of positions for the live cell run 
% or input it as a parameter

%[Npos,~]=folderFilesFromKeyword(direc,'f00_','_t0000_w0000_z0000');% all cyto masks for the frame 0 ( four time groups)'Cyto','{0002}'['Outfile_000' num2str(pos) '_t']


for j = 0:Npos;%length(nums)
   MaxProjTimeGroupsAN(direc,direc2,j,tg,chan)
end

end