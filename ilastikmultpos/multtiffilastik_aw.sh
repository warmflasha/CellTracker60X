#!/bin/bash

i=0
num=1
mkdir /Users/aryeh/Desktop/ilastik_tracking_trial/pxlseg/

files=/Users/aryeh/Desktop/ilastik_tracking_trial/rawdata/*

for f in $files
do
/Applications/ilastik-1.1.6-OSX.app/Contents/MacOS/python /Applications/ilastik-1.1.6-OSX.app/Contents/Resources/ilastik.py --headless  --project=/Users/aryeh/Desktop/coCulture3_f0053.ilp $f

i=$((i+num))

scp /Users/aryeh/Desktop/ilastik_tracking_trial/trial1.h5 ~/Desktop/ilastik_tracking_trial/pxlseg/seg$i.h5

echo $f $i
done



export PATH=/Applications/MATLAB_R2015a.app/bin:$PATH

matlab -nodesktop -nosplash -r ilastikout2peaks\(\'/Users/aryeh/Desktop/ilastik_tracking_trial/pxelseg/\'\); 

exit;



