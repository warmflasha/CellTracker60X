clear all;
dir1 = '.';
nch = 1; % channel no. to be analysed
sn = 3;  % no. of different channels
pos = [1 1 1];
z1 = [16 16 16];
negperc = 65;
negsamp = 2;

for i = 1:sn
    sname{i} = sprintf('sample%d', i);
end

%% Quantifying mRNA 
% Note: each section below can be run only after the previous one has been
% run.
% 
%Spatzcell code begins!
tic;
 
RunSpotRecognitiontest(dir1, z1, pos, sn, nch, sname);

GroupSpotsAndPeakHistsTest(dir1, z1, pos, sn, nch, negsamp, sname, negperc);

GetSingleMrnaIntTest(dir1, z1, pos, sn, nch, sname);

GroupCellSpotsTest(dir1, z1, pos, sn, nch, sname);

toc;

%%
n_ch = [1]; % Channels that are analysed and need to be tabulated. List out all the channels that need to be tabulated.
sn = 3;
%dir1 = pwd;

tabulatemRNAfishnotile(dir1, sn, n_ch);

%%
