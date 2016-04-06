function spotinfo(inputdirec, channels, samplenum)


%%
ch = channels;
samples = samplenum;

newfolder = strcat(inputdirec, '/', 'mrna_information');
mkdir (newfolder)

for i = 1:length(ch)
    filename = strcat(inputdirec, sprintf('/spotsquantify_ch%d/data/FISH_spots_data_new.mat', ch(i)));
    load(filename);
    clear spotinfomat
    
    for j = 1:samples
        spotinfo = spotlist_new{j};
        spotinforeq = [spotinfo(:,14), spotinfo(:,1), spotinfo(:,7), spotinfo(:,8),spotinfo(:,5)/One_mRNA];
        if j == 1
           spotinfomat = spotinforeq;
        else
            spotinfomat = [spotinfomat; spotinforeq];
        end
         
    end
    
    filename = strcat(inputdirec, sprintf('/ch%dallspots.mat', ch(i)));
    save(filename, 'spotinfomat');
end

