

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

for k=1:size(pmaskscyto,3)
    figure(20), subplot(2,3,k),imshow(pmaskscyto(:,:,k),[]);
    
end


%%

Lnuc = imgfiles(j).NucMask;
j = 2;
for k=1:size(Lnuc,3)
 figure(7+j),subplot(1,size(Lnuc,3),k),imshow(Lnuc(:,:,k),[]);%(7),subplot(2,size(Lnuc,3),k),
end
figure(7+j), subplot(1,size(Lnuc,3),2), hold on
plot(outdat(:,1),outdat(:,2),'*r','markersize',15);
figure(7+j), subplot(1,size(Lnuc,3),2), hold on
text(outdat(:,1)+10,outdat(:,2),num2str(outdat(:,6)./outdat(:,7)),'Color','m');

%%
for k=1:size(pmasks,3)
    figure(5), subplot(2,3,k),imshow(pmasks(:,:,k),[]);%pmasks
    
end

%%
j = 40;
Lcytofin = imgfilescyto(j).NucMask;

for k=1:size(Lcytofin,3)
 figure(8),subplot(2,size(Lcytofin,3),k),imshow(Lcytofin(:,:,k),[]);
end
%%

for k=1:5
 figure(9),subplot(2,3,k),imshow(pnuc(:,:,k),[]);
end

for k=1:5
 figure(11),subplot(2,3,k),imshow(pcyto(:,:,k),[]);
end
%%
j = 40;
Lnuc = imgfiles(j).NucMask;
for k=1:size(Lnuc,3)
figure(16),subplot(1,size(Lnuc,3),k),showMaskWithNumber(Lnuc(:,:,k))
end