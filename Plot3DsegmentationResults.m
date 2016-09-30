

%%
% plot the layers of the variables, returned during 3D segmentation
% useful for troubleshooting

for k=1:size(newmask_lbl,3)
    figure(10), subplot(2,size(newmask_lbl,3),k),imshow(newmask_lbl(:,:,k),[]);
    figure(11), subplot(2,size(newmask_lbl,3),k),imshow(maskzcyto_lbl(:,:,k),[]);
end
%%
for k=1:size(icyto_new,3)
    figure(3), subplot(2,size(icyto_new,3),k),imshow(icyto_new(:,:,k),[]);
    figure(4), subplot(2,size(icyto_new,3),k),imshow(inuc_new(:,:,k),[]);
end
%%
for k=1:size(icyto,3)
    figure(5), subplot(2,size(icyto,3),k),imshow(icyto(:,:,k),[]);
    figure(6), subplot(2,size(icyto,3),k),imshow(inuc(:,:,k),[]);
end

%%

for k=1:size(pmasks,3)
    figure(22), subplot(2,size(pmasks,3),k),imshow(pmasks(:,:,k),[]);
    
end
%%
for k=1:size(smasks,3)
    figure(9), subplot(2,size(pmasks,3),k),imshow(smasks(:,:,k),[]);%smaskssmasks(:,:,k)
    
end

%%
for k=1:size(masktmp,2)
    figure(22), subplot(2,size(masktmp,2),k),imshow(masktmp{k},[]);
    
end
%%
for k=1:size(MaskFin2,2)
    figure(22), subplot(2,size(MaskFin2,2),k),imshow(MaskFin2{k},[]);
    
end

%%
j =1; 
f = 10;
Lnuc = imgfiles(j).NucMask;
%Lnuc = imgfilescyto(j).Cyto;
%j = 2; % peaks{j}
for k=1:size(Lnuc,3)
 figure(f),subplot(1,size(Lnuc,3),k),imshow(Lnuc(:,:,k),[]);%(7),subplot(2,size(Lnuc,3),k),
end
figure(f),hold on
plot(peaks{j}(:,1),peaks{j}(:,2),'*r','markersize',15);
figure(f),  hold on
text(peaks{j}(:,1)+10,peaks{j}(:,2),num2str(peaks{j}(:,8)./peaks{j}(:,5)),'Color','m'); % cdx2:dapi
title('CDX2/DAPI')



%%
for k=1:size(Lnuc,3)
 figure(9),subplot(2,size(Lnuc,3),k),imshow(Lnuc(:,:,k),[]);
end
%%
%Lcytofin;
for k=1:size(Lcytofin,3)
 figure(11),subplot(2,size(Lcytofin,3),k),imshow(Lcytofin(:,:,k),[]);
end

%%

for k=1:size(I2proc,3)
 figure(18),subplot(2,size(I2proc,3),k),imshow(I2proc(:,:,k),[]);
end
%%


for k=1:size(pnuc,3)
 figure(9),subplot(2,3,k),imshow(pnuc(:,:,k),[]);
end
%%
for k=1:size(pnuc,3)
 figure(11),subplot(2,3,k),imshow(pcyto(:,:,k),[]);
end
%%
for j=1
    plane =1;                 
Lnuc = imgfiles(j).NucMask(:,:,plane);%uncompressBxinaryImg(imgfiles(j).NucMask);
Lcytofin = imgfilescyto(j).Cyto(:,:,plane);%uncompressBinaryImg(imgfilescyto(j).Cyto);
for k=1:2
figure(j),subplot(1,2,1),showMaskWithNumber(Lnuc)
figure(j),subplot(1,2,2),imshow(Lcytofin,[]);
end
end

