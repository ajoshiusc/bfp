#!/bin/bash

exe_name=$0
exe_dir=`dirname "$0"`

# If MCR R2015b is installed in a non-default location, define correct path 
# on next line and uncomment it (remove the leading "#")
#BrainSuiteMCR="/path/to/your/MCR";

if [ -z "$BrainSuiteMCR" ]; then
  if [ -e /usr/local/MATLAB/MATLAB_Runtime/v90 ]; then
    BrainSuiteMCR="/usr/local/MATLAB/MATLAB_Runtime/v90";
  elif [ -e /usr/local/MATLAB/R2015b/runtime ]; then
    BrainSuiteMCR="/usr/local/MATLAB/R2015b";
  else
    echo
    echo "Could not find Matlab 2015b with Matlab Compiler or MCR 2015b (v7.17)."
    echo "Please install the Matlab 2015b MCR from MathWorks at:"
    echo
    echo "http://www.mathworks.com/products/compiler/mcr/"
    echo 
    echo "If you already have Matlab 2015b with the Matlab Compiler or MCR 2015b"
    echo "installed, please edit ${exe_name} by uncommenting and editing the line:"
    echo "#BrainSuiteMCR=\"/path/to/your/MCR\";"
    echo "(replacing /path/to/your/MCR with the path to your Matlab or MCR installation)"
    echo "near the top of the file"
    echo
    exit 78
  fi
fi

read -d '' usage <<EOF

  generateVolGOrdfMRI : Generate Volume Grayordinates.
  Authored by Anand A. Joshi, Signal and Image Processing Institute
  Department of Electrical Engineering, Viterbi School of Engineering, USC
  
  Usage: generateVolGOrdfMRI.sh GOrdVolIndFile, subbasename, fmri2standard, GOrdVolFile

    GOrdVolIndFile: the indices of volume grayordinates. It is a .mat file.
    subbasename: subject base name
    fmri2standard: fMRI registered to the surface
    GOrdVolFile = Grayordinate volume file (.mat)
EOF

# Parse inputs
if [ $# -lt 4 ]; then
  echo
  echo "$usage";
  echo
  exit;
fi

GOrdVolIndFile=$1
subbasename=$2 
fmri2standard=$3
GOrdVolFile = $4



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
${exe_dir}/generateVolGOrdfMRI.sh "${GOrdVolIndFile}" "${subbasename}" "${fmri2standard}" "${GOrdVolFile}" 

exit
