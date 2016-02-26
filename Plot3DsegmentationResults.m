

%%
% plot the layers of the variables, returned during 3D segmentation
% useful for troubleshooting

for k=1:size(newmask_lbl,3)
    figure(1), subplot(2,size(newmask_lbl,3),imshow(newmask_lbl(:,:,k),[]);
    figure(2), subplot(2,size(newmask_lbl,3),imshow(maskzcyto_lbl(:,:,k),[]);
end

for k=1:size(icyto_new,3)
    figure(3), subplot(2,size(icyto_new,3),k),imshow(icyto_new(:,:,k),[]);
    figure(4), subplot(2,size(icyto_new,3),k),imshow(inuc_new(:,:,k),[]);
end

for k=1:5
    figure(5), subplot(2,3,k),imshow(icyto(:,:,k),[]);
    figure(6), subplot(2,3,k),imshow(inuc(:,:,k),[]);
end

%%
j = 2;
Lnuc = imgfiles(j).NucMask;
for k=1:size(Lnuc,3)
 figure(7),subplot(1,size(Lnuc,3),k),imshow(Lnuc(:,:,k),[]);%(7),subplot(2,size(Lnuc,3),k),
end
figure(7), subplot(1,size(Lnuc,3),2), hold on
plot(peaks{j}(:,1),peaks{j}(:,2),'*r','markersize',15);
figure(7), subplot(1,size(Lnuc,3),2), hold on
text(peaks{j}(:,1)+10,peaks{j}(:,2),num2str(peaks{j}(:,6)./peaks{j}(:,7)),'Color','m');

%%

Lcytofin = imgfilescyto(j).NucMask;

for k=1:size(Lcytofin,3)
 figure(8),subplot(2,size(Lcytofin,3),k),imshow(Lcytofin(:,:,k),[]);
end
%%

for k=1:5
 figure(9),subplot(1,5,k),imshow(inuc(:,:,k),[]);
end

for k=1:5
 figure(11),subplot(1,5,k),imshow(pcyto(:,:,k),[]);
end
%%
j = 10;
Lnuc = imgfiles(j).NucMask;
for k=1:size(Lnuc,3)
figure(15),subplot(1,size(Lnuc,3),k),showMaskWithNumber(Lnuc(:,:,k))
end