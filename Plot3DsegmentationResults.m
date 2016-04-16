

%%
% plot the layers of the variables, returned during 3D segmentation
% useful for troubleshooting

for k=1:size(newmask_lbl,3)
    figure(10), subplot(2,size(newmask_lbl,3),k),imshow(newmask_lbl(:,:,k),[]);
    figure(11), subplot(2,size(newmask_lbl,3),k),imshow(maskzcyto_lbl(:,:,k),[]);
end

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
    figure(22), subplot(2,3,k),imshow(pmasks(:,:,k),[]);
    
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
j = 10;
Lnuc = imgfiles(j).NucMask;
%j = 2; % peaks{j}
for k=1:size(Lnuc,3)
 figure(7),imshow(Lnuc(:,:,k),[]);%(7),subplot(2,size(Lnuc,3),k),
end
figure(7),hold on
plot(peaks{j}(:,1),peaks{j}(:,2),'*r','markersize',15);
figure(7),  hold on
text(peaks{j}(:,1)+10,peaks{j}(:,2),num2str(peaks{j}(:,6)./peaks{j}(:,7)),'Color','m');


%%
for k=1:size(smasks,3)
    figure(8), subplot(2,3,k),imshow(smasks(:,:,k),[]);%smasks
    
end

%%
for k=1:size(Lnuc,3)
 figure(9),subplot(2,size(Lnuc,3),k),imshow(Lnuc(:,:,k),[]);
end
%%
Lcytofin;
for k=1:size(Lcytofin,3)
 figure(11),subplot(2,size(Lcytofin,3),k),imshow(Lcytofin(:,:,k),[]);
end

%%

for k=1:size(pnuc,3)
 figure(9),subplot(2,3,k),imshow(pnuc(:,:,k),[]);
end

for k=1:size(pnuc,3)
 figure(11),subplot(2,3,k),imshow(pcyto(:,:,k),[]);
end
%%
for j=1:10;
                    
Lnuc = imgfiles(j).NucMask;%uncompressBinaryImg(imgfiles(j).NucMask);
Lcytofin = imgfilescyto(j).Cyto;%uncompressBinaryImg(imgfilescyto(j).Cyto);
for k=1:2
figure(j),subplot(1,2,1),showMaskWithNumber(Lnuc)
figure(j),subplot(1,2,2),imshow(Lcytofin,[]);
end
end
