#!/bin/bash

exe_name=$0
exe_dir=`dirname "$0"`

# If MCR R2023a is installed in a non-default location, define correct path 
# on next line and uncomment it (remove the leading "#")
#BrainSuiteMCR="/path/to/your/MCR";

if [ -z "$BrainSuiteMCR" ]; then
  if [ -e /usr/local/MATLAB/MATLAB_Runtime/R2023a ]; then
    BrainSuiteMCR="/usr/local/MATLAB/MATLAB_Runtime/R2023a";
  elif [ -e /usr/local/MATLAB/R2023a/runtime ]; then
    BrainSuiteMCR="/usr/local/MATLAB/R2023a";
  else
    echo
    echo "Could not find Matlab 2023a with Matlab Compiler or Matlab 2023a (9.14)."
    echo "Please install the Matlab 2023a MCR from MathWorks at:"
    echo
    echo "https://www.mathworks.com/products/compiler/matlab-runtime.html"
    echo 
    echo "If you already have Matlab 2023a with the Matlab Compiler or Matlab 2023a"
    echo "installed, please edit ${exe_name} by uncommenting and editing the line:"
    echo "#BrainSuiteMCR=\"/path/to/your/MCR\";"
    echo "(replacing /path/to/your/MCR with the path to your Matlab or MCR installation)"
    echo "near the top of the file"
    echo
    exit 78
  fi
fi

read -d '' usage <<EOF

  bfp : BrainSuite fMRI pipeline (bfp)
  Authored by Anand A. Joshi, Signal and Image Processing Institute,
              Soyoung Choi, USC
  Department of Electrical and Computer Engineering, Viterbi School of Engineering, USC
  website: http://brainsuite.org/bfp
  
  Usage: bfp configfile t1 fmri studydir subid sessionid TR 

  configfile: Configuration file that you edited during installation, (<bfpdir>/supp_data/config.ini)
  t1: T1 image in NIFTI-1 format e.g. 'sub01-T1w.nii.gz'
  fmri: string or cell array of fMRI NIFTI-1 formatted files e.g. 'sub_01_rest_bold.nii.gz'
  studydir, subid: Outputs will be saved in the studydir/subid for each subject e.g. studydir='mystudy', subid='sub-01'
  sessionid: string or cell array of session id for the input fmris e.g. 'rest' for resting state

  e.g. bfp /home/ajoshi/bfp/supp_data/config.ini /home/ajoshi/sub-01_T1w.nii.gz /home/ajoshi/sub-01_rest.nii.gz /home/ajoshi/mystudy sub-01 rest 2
    TR is in sec. If you don't know the TR, set it to 0 and bfp will try to read it from the header of fmri nifti file

EOF

# Parse inputs
if [ $# -lt 7 ]; then
  echo
  echo "$usage";
  echo
  exit;
fi

configfile=$1
t1=$2
fmri=$3
studydir=$4
subid=$5
sessionid=$6
TR=$7



# Set up path for MCR applications.
PATH=${exe_dir}:${PATH} ;
LD_LIBRARY_PATH=.:${BrainSuiteMCR}/runtime/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BrainSuiteMCR}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BrainSuiteMCR}/sys/os/glnxa64;
MCRJRE=${BrainSuiteMCR}/sys/java/jre/glnxa64/jre/lib/amd64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ; 
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;  
XAPPLRESDIR=${BrainSuiteMCR}/X11/app-defaults ;
export PATH;
export LD_LIBRARY_PATH;
export XAPPLRESDIR;


# Run the BFP Sequence
${exe_dir}/bfp "${configfile}" "${t1}" "${fmri}" "${studydir}" "${subid}" "${sessionid}" "${TR}" 

exit
