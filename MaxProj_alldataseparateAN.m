function  MaxProj_alldataseparateAN(direc,direc2,tg,chan)
% direc = directory with raw images

% tg = number fo time points to be used from the dataset
% chan = corresponds to ff.w(chan), e.g. ff.w(1) = nuclear channel
% need to run this function separately for each channel
% direc2 = where to save the output projections

%filename = getAndorFileName(files,pos,time,z,w)

ff=readAndorDirectory(direc);
nz = size(ff.z,2);
if isempty(tg)
    tg = size(ff.t,2);
end
Npos = length(ff.p);

for j= 1:(Npos)
    for xx = 1:(tg)
        %     %-------- the very first live cell imaging  dataset does not have all the z slices , random ones are missing, so the max
        %     %projection screws up
        %    imgfiles  = dir(direc);
        %    filesInDir = imgfiles(~([imgfiles.isdir]));
        %    for h = 1:nz
        %        if j <10 && xx <10
        %          stringToBeFound = (['SingleCellSignalingAN_t000'  num2str(xx) '_f000'  num2str(j) '_z000'  num2str(h-1) '_w000'  num2str(ff.w(chan)) '.tif']);
        %        end
        %        if j >10 && xx >10
        %          stringToBeFound = (['SingleCellSignalingAN_t00'  num2str(xx) '_f00'  num2str(j) '_z000'  num2str(h-1) '_w000'  num2str(ff.w(chan)) '.tif']);
        %        end
        %    for jj = 1:size(filesInDir,1)
        %    found = strfind(filesInDir(jj).name,stringToBeFound);
        %    if isempty(found)
        %        not_there(j,1:4) = [xx j h ff.w(chan)];
        %    end
        %    end
        %    end
        
        %------------
%         try
            if isempty(ff.z)
                filename = getAndorFileName(ff,j,xx,0,ff.w(chan));
                max_img = imread(filename);
            end
            
            
            for ii=1:length(ff.z)
                filename = getAndorFileName(ff,j,xx,ff.z(ii),ff.w(chan));
                
                img_now = imread(filename);
                
                if ii==1
                    max_img=img_now;
                else
                    max_img=max(img_now,max_img);
                end
            end
            
%         catch exception
%             rethrow(exception)
        end
    end
    if j < 10
        imwrite(max_img,[direc2 'Maxprojection_f000' num2str(j) '_t000' num2str(1) '_w000' num2str(ff.w(chan)) '.tif'],'writemode','append','Compression','none');%
    end
    if j >= 10
        imwrite(max_img,[direc2 'Maxprojection_f00' num2str(j) '_t000' num2str(1) '_w000' num2str(ff.w(chan)) '.tif'],'writemode','append','Compression','none');%
    end
end
end

