#!/bin/bash

i=0
num=1
mkdir /Users/warmflash2/Desktop/train_ilastik/pxlseg/

files=/Users/warmflash2/Desktop/train_ilastik/multtif/*

for f in $files
do
./ilastik-1.1.6-OSX.app/Contents/MacOS/python ilastik-1.1.6-OSX.app/Contents/Resources/ilastik.py --headless  --project=/Users/warmflash2/Desktop/train_ilastik/multitiffpxlclass_2.ilp $f

i=$((i+num))

scp /Users/warmflash2/Desktop/train_ilastik/seglabel.h5 /Users/warmflash2/Desktop/train_ilastik/pxlseg/seg$i.h5

echo $f $i
done



export PATH=/Applications/MATLAB_R2015a.app/bin:$PATH

matlab -nodesktop -nosplash -r ilastikout2peaks\(\'/Users/warmflash2/Desktop/train_ilastik/pxlseg\'\); 

exit;



