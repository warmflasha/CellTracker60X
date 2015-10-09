#!/bin/bash

i=0
num=1

ndir=/Users/warmflash2/Desktop/train_ilastik/pxlseg/

prjtpath=/Users/warmflash2/Desktop/train_ilastik/multitiffpxlclass_2.ilp

savepath=/Users/warmflash2/Desktop/train_ilastik/seglabel.h5

segchannel=1

mkdir $ndir

samplepath=/Users/warmflash2/Desktop/train_ilastik/multtif/;

files=/Users/warmflash2/Desktop/train_ilastik/multtif/*

for f in $files
do
./ilastik-1.1.6-OSX.app/Contents/MacOS/python ilastik-1.1.6-OSX.app/Contents/Resources/ilastik.py --headless  --project=$prjtpath $f

i=$((i+num))

scp $savepath ${ndir}fileseg${i}.h5

echo $f $i
done



export PATH=/Applications/MATLAB_R2015a.app/bin:$PATH


matlab -nodesktop -nosplash -r ilastikout2peaks\(\'$ndir\',\'$samplepath\',$segchannel\); 

exit;



