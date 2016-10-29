
function plotselecteddata(dir,matfile1,matfile2,index1,param1,chanmax)

% dir = directory with the outfile
% matfile1 = string name of the mat file with the plate1 data
% matfile2 = string name of the mat file with the selected image  number
% saved in the variable imN3
% index1 = columns of colonies.data to plot, second should be dapi, if onlt
% one is supplied, that col will be plotted unnormalized
% param1 = string name to characterize the column you are plotting (e.g.
% Sox2)

dapimax =5000;%now used as the area thresh in colony analysis; dapimax is set to max 60000 within the generalized mean function

usemeandapi =[];
C= {'c','r'};
flag = 0;
clear rawdata1
clear totalcells


if isempty(matfile2);
    flag1 = 1; % to plot all images
    imN3 = [];
else
    load(matfile2,'imN3');
    flag1 = 0;
end

[rawdata1,totalcells] =  Intensity_vs_ColSizeSelectImages(matfile1,[],dir,index1,param1,dapimax,chanmax,usemeandapi,flag,imN3,flag1);

hold on,figure(1), plot(rawdata1{1}(1:8),'-*','color',C{1},'markersize',18,'linewidth',2);
h1 = figure(1);
h1.CurrentAxes.FontSize = 18;
legend('Selected images');
xlim([0 8]);
ylim([0 7]);
xlabel('Colony Size');
ylabel(['Expression of ' num2str(param1) ' marker']);

hold on,figure(2), plot(totalcells,'-*','color',C{2},'markersize',18,'linewidth',2);
h = figure(2);
h.CurrentAxes.FontSize = 18;
legend('Selected images');
xlim([0 8]);
ylim([0 1000])
ylabel('Total Cells');
xlabel('Colony Size');

end