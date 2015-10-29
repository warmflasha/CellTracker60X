function [outdatnuc,outdatcyto,Lnuc,Lcytofin] = IlastikplusWatershed_AN(ilastikfile,ilastikfilecyto,pos,zplane,direc,img,flag)

ff=readAndorDirectory(direc);
chan = ff.w;

global userParam;
userParam.gaussRadius = 10;
userParam.gaussSigma = 3;
userParam.small_rad = 3;
userParam.presubNucBackground = 1;
userParam.backdiskrad = 300;
areanuclow = 1300;
areanuchi = 9000;

info = h5info(ilastikfile);
info.Datasets;   
infocyto = h5info(ilastikfilecyto);
infocyto.Datasets ;

% check these to make sure the dataset name is correct for the 'h5read' function

data = h5read(ilastikfile,'/exported_data');
data = squeeze(data);

datacyto = h5read(ilastikfilecyto,'/exported_data');
datacyto = squeeze(datacyto);
%time = size(data,3);   
Lnuc = [];
Lcytofin = [];
    
Lnuc = data(:,:,img) <2;  % need to leave only the nuclei masks and make the image binary
Lnuc =  bwareafilt(Lnuc',[areanuclow areanuchi]);

filename = getAndorFileName(ff,pos,ff.t(img),ff.z(zplane),chan(2)); % has to be channel 2 since all the masks should be applied to the gfp channel
I2 = imread(filename);

Lcyto = WatershedsegmCytoplasm_AW(Lnuc,I2);% get the cyto mask using watershed

LcytoIl = datacyto(:,:,img) >1; 
LcytoIl = (LcytoIl');
LcytoIl = LcytoIl & ~ Lnuc;                        % cyto masks initially include both nuclei+cyto, so need to eliminate nuc
Lcytofin = LcytoIl & Lcyto; 

% at this point should have an array of nuc and cyto masks(from ilastik
% and watershed respectively)

I2proc = imopen(I2,strel('disk',userParam.small_rad)); % remove small bright stuff
I2proc = smoothImage(I2proc,userParam.gaussRadius,userParam.gaussSigma); %smooth
I2proc = presubBackground_self(I2proc);

%get the NUCLEAR mean intensity for each labeled object
cc_nuc = bwconncomp(Lnuc,8);
statsnuc = regionprops(cc_nuc,I2proc,'Area','Centroid','PixelIdxList','MeanIntensity');
badinds = [statsnuc.Area] < 200; 
statsnuc(badinds) = [];
anuc = [statsnuc.Area]';
aa = [statsnuc.Centroid];
xnuc = aa(1:2:end)';
ynuc = aa(2:2:end)';
nucmeanInt = [statsnuc.MeanIntensity]';

%get the cytoplasmic mean intensity for each labeled object
cc_cyto = bwconncomp(Lcytofin);
statscyto = regionprops(cc_cyto,I2proc,'Area','Centroid','PixelIdxList','MeanIntensity');
badinds = [statscyto.Area] < 200; 
statscyto(badinds) = [];
acyto = [statscyto.Area]';
aa = [statscyto.Centroid];
xcyto = aa(1:2:end)';
ycyto = aa(2:2:end)';
cytomeanInt = [statscyto.MeanIntensity]';


outdatnuc = [xnuc ynuc anuc nucmeanInt ];
outdatcyto = [xcyto ycyto acyto cytomeanInt];
% ind = find(outdatnuc(:,5)<1700);
% if ~ isempty(ind); 
% outdatnuc(ind,:)= [];
% end

  if flag == 1
      figure, subplot(1,3,1),imshow(I2proc,[]);hold on
      subplot(1,3,2),imshow(Lcytofin);hold on
      subplot(1,3,3),imshow(Lnuc);hold on
  end
 
  
end



