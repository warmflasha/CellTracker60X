#!/bin/bash
# ilastik segmentation output to peaks: input- multiple time points, multichannel, multitiff files.
#ndir: new directory where segmentation results are saved
#projectpath: path corresponding to ilastik project that was usecd to train a sample dataset
#savepath: path corresponding to ilastik output export directory
#segchannel: channel used for segmenting images
#samplepath: path corresponding to sample images that are yet to be segmented - serves as an input to the matlab function
#files: image files in the samplepath that are to be segmented 


i=0
num=1

ndir=/Users/warmflashlab/Desktop/Dec31setIlastikMasks_headless_DiffW0/

projectpath=/Users/warmflashlab/Desktop/Dec31imging_TrainingSet/Nuc_Training.ilp

savepath=/Users/warmflashlab/Desktop/Dec31imging_TrainingSet/\{NucMasks\}_\{P\}.h5

#segchannel=1

mkdir $ndir

samplepath=/Users/warmflashlab/Desktop/MaxProjections_december31Diff/W0/;

files=/Users/warmflashlab/Desktop/MaxProjections_december31Diff/W0/*;

for f in $files
do
/Applications/ilastik-1.1.8-OSX.app/Contents/MacOS/python /Applications/ilastik-1.1.8-OSX.app/Contents/Resources/ilastik.py --headless  --project=$projectpath $f

i=$((i+num))

scp $savepath ${ndir}NucMaskDiffDec31set_tg${i}.h5

echo $f $i
done

exit;



