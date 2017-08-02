#!/bin/bash

if [ -z "$1" ]; then
  echo "BrainSuite v17a cortical surface extraction script"
  echo "Authored by David W Shattuck (http://shattuck.bmap.ucla.edu)"
  echo "UCLA Brain Mapping Center"
  echo "for more information, please see: http://brainsuite.org"
  echo "usage: $0 input_image"
  echo "	where input_image is a raw MRI file, in .nii, .img, .nii.gz, or .img.gz format."
  exit
fi

# Override this if you want to call the script in a different directory than the binaries reside
# This may have to be hardcoded to work on an SGE grid. You can also set an environment variable to do this.
if [ -z $BrainSuiteDir ]; then
	BrainSuiteDir=$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )")
fi;
BrainSuiteBin="${BrainSuiteDir}/bin/"
echo "Using BrainSuite binaries and data in $BrainSuiteDir"

filename=$1
basename=${filename%.gz}
basename=${basename%.hdr}
basename=${basename%.img}
basename=${basename%.nii}

# change EXT to produce different outputs, e.g., uncompressed nifti
EXT=nii.gz

echo working with basename: $basename

# use the VERBOSE flag to change the output option for all of the programs
# VERBOSE="-v 1 --timer"
# These values are set to match the defaults used in the GUI version of BrainSuite16a1
BSEOPTIONS="-n 3 -d 25 -s 0.64 -p --trim --auto"
BFCOPTIONS="-L 0.5 -U 1.5"
# PVCOPTIONS="-l 0.1"
ATLAS="${BrainSuiteDir}/atlas/brainsuite.icbm452.lpi.v08a.img"
ATLASLABELS="${BrainSuiteDir}/atlas/brainsuite.icbm452.v15a.label.img"
ATLASES="--atlas ${ATLAS} --atlaslabels ${ATLASLABELS}"
# using the centroids option can improve the initial alignment during the cerebrum labeling step
# CEREBROOPTIONS="--centroids"

# you can specify an options file in the directory from which you launch the script
if [ -f BrainSuiteOptions ]; then source BrainSuiteOptions; fi;

echo ---- Running BrainSuite Cortical Surface Extraction Sequence ----
#${BrainSuiteBin}bse $VERBOSE -i ${filename} -o ${basename}.bse.$EXT --mask ${basename}.mask.$EXT $BSEOPTIONS
if [ $? -ne 0 ]; then echo "cortical extraction halted because bse failed to run or return an error."; exit 1; fi;
${BrainSuiteBin}bfc $VERBOSE -i ${basename}.bse.$EXT -o ${basename}.bfc.$EXT $BFCOPTIONS $BFCFILES
if [ $? -ne 0 ]; then echo "cortical extraction halted because bfc failed to run or return an error."; exit 1; fi;
${BrainSuiteBin}pvc $VERBOSE -i ${basename}.bfc.$EXT -o ${basename}.pvc.label.$EXT -f ${basename}.pvc.frac.$EXT $PVCOPTIONS
if [ $? -ne 0 ]; then echo "cortical extraction halted because pvc failed to run or return an error."; exit 1; fi;
${BrainSuiteBin}cerebro $VERBOSE $ATLASES -i ${basename}.bfc.$EXT \
	-l ${basename}.hemi.label.$EXT -m ${basename}.mask.$EXT -o ${basename}.cerebrum.mask.$EXT \
	-c 2 --air ${basename}.air --warp ${basename}.warp ${CEREBROOPTIONS}
if [ $? -ne 0 ]; then echo "cortical extraction halted because cerebro failed to run or return an error."; exit 1; fi;
${BrainSuiteBin}cortex $VERBOSE -f ${basename}.pvc.frac.$EXT -h ${basename}.hemi.label.$EXT -o ${basename}.init.cortex.mask.$EXT -a -w -p 50
if [ $? -ne 0 ]; then echo "cortical extraction halted because cortex failed to run or return an error."; exit 1; fi;
${BrainSuiteBin}scrubmask $VERBOSE -i ${basename}.init.cortex.mask.$EXT -o ${basename}.cortex.scrubbed.mask.$EXT -f 0 -b 2
if [ $? -ne 0 ]; then echo "cortical extraction halted because scrubmask failed to run or return an error."; exit 1; fi;
${BrainSuiteBin}tca $VERBOSE -i ${basename}.cortex.scrubbed.mask.$EXT -o ${basename}.cortex.tca.mask.$EXT -m 2500  --delta 20 $TCASETTINGS
if [ $? -ne 0 ]; then echo "cortical extraction halted because tca failed to run or return an error."; exit 1; fi;
${BrainSuiteBin}dewisp $VERBOSE -i ${basename}.cortex.tca.mask.$EXT -o ${basename}.cortex.dewisp.mask.$EXT $DEWISPSETTINGS
if [ $? -ne 0 ]; then echo "cortical extraction halted because dewisp failed to run or return an error."; exit 1; fi;
${BrainSuiteBin}dfs $VERBOSE -i ${basename}.cortex.dewisp.mask.$EXT -o ${basename}.inner.cortex.dfs -n 10 -a 0.5 -w 5.0
if [ $? -ne 0 ]; then echo "cortical extraction halted because dfs failed to run or return an error."; exit 1; fi;
${BrainSuiteBin}pialmesh $VERBOSE -f ${basename}.pvc.frac.$EXT -i ${basename}.inner.cortex.dfs -m ${basename}.cerebrum.mask.$EXT -o ${basename}.pial.cortex.dfs \
	-n 100 -r 1 -s 0.40 -t 1.05 --max 20 --smooth 0.025 --interval 10 --nc 0.2
if [ $? -ne 0 ]; then echo "cortical extraction halted because pialmesh failed to run or return an error."; exit 1; fi;
${BrainSuiteBin}hemisplit $VERBOSE -i ${basename}.inner.cortex.dfs -l ${basename}.hemi.label.$EXT \
	--left ${basename}.left.inner.cortex.dfs --right ${basename}.right.inner.cortex.dfs \
	-p ${basename}.pial.cortex.dfs -pl ${basename}.left.pial.cortex.dfs -pr ${basename}.right.pial.cortex.dfs
if [ $? -ne 0 ]; then echo "cortical extraction halted because hemisplit failed to run or return an error."; exit 1; fi;

exit 0

