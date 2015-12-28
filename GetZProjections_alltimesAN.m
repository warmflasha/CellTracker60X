function  GetZProjections_alltimesAN(direc,direc2,pos,tg,chan)

% direc = directory with raw images
% pos = the position to be processed
% tg = time group
% chan = corresponds to ff.w(chan), e.g. ff.w(1) = nuclear channel
% direc2 = where to save the output projections
% tg = a vector with the number of separate time groups


for j = 1:length(tg)
    for k=1:length(chan)
        MaxProjTimeGroupsAN(direc,direc2,pos,tg(j),chan(k))
        
    end
end

end