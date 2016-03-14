

%%
% plot the layers of the variables, returned during 3D segmentation
% useful for troubleshooting

for k=1:size(newmask_lbl,3)
    figure(1), subplot(2,size(newmask_lbl,3),k),imshow(newmask_lbl(:,:,k),[]);
    figure(2), subplot(2,size(newmask_lbl,3),k),imshow(maskzcyto_lbl(:,:,k),[]);
end

for k=1:size(icyto_new,3)
    figure(3), subplot(2,size(icyto_new,3),k),imshow(icyto_new(:,:,k),[]);
    figure(4), subplot(2,size(icyto_new,3),k),imshow(inuc_new(:,:,k),[]);
end
%%
for k=1:5
    figure(5), subplot(2,3,k),imshow(icyto(:,:,k),[]);
    figure(6), subplot(2,3,k),imshow(inuc(:,:,k),[]);
end

%%

for k=1:size(pmasks,3)
    figure(21), subplot(2,3,k),imshow(pmasks(:,:,k),[]);
    
end


%%

Lnuc = uncompressBinaryImg(imgfiles(j).NucMask);
%j = 2; % peaks{j}
for k=1:size(Lnuc,3)
 figure(7+j),subplot(1,size(Lnuc,3),k),imshow(Lnuc(:,:,k),[]);%(7),subplot(2,size(Lnuc,3),k),
end
figure(7+j), subplot(1,size(Lnuc,3),2), hold on
plot(peaks{j}(:,1),peaks{j}(:,2),'*r','markersize',15);
figure(7+j), subplot(1,size(Lnuc,3),2), hold on
text(peaks{j}(:,1)+10,peaks{j}(:,2),num2str(peaks{j}(:,6)./peaks{j}(:,7)),'Color','m');

%%
for k=1:size(pmasks,3)
    figure(5), subplot(2,3,k),imshow(pmasks(:,:,k),[]);%pmasks
    
end


%%
for k=1:size(Lnuc,3)
 figure(9),subplot(2,size(Lnuc,3),k),imshow(Lnuc(:,:,k),[]);
end
%%
Lcytofin = LcytoIl;
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
j = 10;
Lnuc = imgfiles(j).NucMask;%uncompressBinaryImg(imgfiles(j).NucMask);
Lcytofin = imgfilescyto(j).Cyto;%uncompressBinaryImg(imgfilescyto(j).Cyto);
for k=1:2
figure(16),subplot(1,2,1),showMaskWithNumber(Lnuc)
figure(16),subplot(1,2,2),imshow(Lcytofin,[]);
end
