%% read file
% ff=readAndorDirectory('.');
% %nuc = andorMaxIntensity(ff,0,0,1);
% fn=getAndorFileName(ff,9,0,0,0);
% nuc = imread(fn);
nuc = imread('SingleCellSignalingAN_t0000_f0000_z0003_w0000.tif');
%nuc = I;
nuc_o = nuc;
%% preprocess
global userParam;
userParam.gaussRadius = 10;
userParam.gaussSigma = 3;
userParam.small_rad = 3;
userParam.presubNucBackground = 1;
userParam.backdiskrad = 300;

nuc = imopen(nuc,strel('disk',userParam.small_rad)); % remove small bright stuff
nuc = smoothImage(nuc,userParam.gaussRadius,userParam.gaussSigma); %smooth
nuc =presubBackground_self(nuc);
%%  Normalize image
diskrad = 100;
low_thresh = 500;

nuc(nuc < low_thresh)=0;
norm = imdilate(nuc,strel('disk',diskrad));
normed_img = im2double(nuc)./im2double(norm);
normed_img(isnan(normed_img))=0;
%% gradient image
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(normed_img), hy, 'replicate');
Ix = imfilter(double(normed_img), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
%% circle find and display
%[cc, rr, met]=imfindcircles(gradmag,[20 40],'Method','TwoStage','Sensitivity',0.95);
[cc, rr, met]=imfindcircles(gradmag,[20 40],'Method','TwoStage','Sensitivity',0.999);
figure; imshow(nuc,[]); hold;
for ii=1:length(cc)
    drawcircle(cc(ii,:),rr(ii),'m');
end

%throw out circles with nothing inside
cavg = zeros(length(rr),1);
for ii=1:length(rr)
[cavg(ii), mm]=averageImageInCircle(nuc,floor(cc(ii,:)),rr(ii));
end
badinds = cavg < 1000;
cc(badinds,:)=[]; rr(badinds,:)=[];

% convert circlees to cells (will merge close circles) 
cen = circles2cells(cc,rr);
%% display results
figure; 
 
imshow(nuc_o,[]);hold on; plot(cen(:,1),cen(:,2),'r*');
title('Original Image with cells identified');
% 

% figure; 
% imshow(gradmag,[]);
% hold on;
% for ii=1:length(rr)
%     drawcircle(cc(ii,:),rr(ii),'m');
% end
% title('gradient image with circles');