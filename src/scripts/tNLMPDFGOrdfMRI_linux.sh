#!/bin/bash

exe_name=$0
exe_dir=`dirname "$0"`

# If MCR R2019b is installed in a non-default location, define correct path 
# on next line and uncomment it (remove the leading "#")
#BrainSuiteMCR="/path/to/your/MCR";

if [ -z "$BrainSuiteMCR" ]; then
  if [ -e /usr/local/MATLAB/MATLAB_Runtime/v97 ]; then
    BrainSuiteMCR="/usr/local/MATLAB/MATLAB_Runtime/v97";
  elif [ -e /usr/local/MATLAB/R2019b/runtime ]; then
    BrainSuiteMCR="/usr/local/MATLAB/R2019b";
  else
    echo
    echo "Could not find Matlab 2019b with Matlab Compiler or MCR 2019b (v7.17)."
    echo "Please install the Matlab 2019b MCR from MathWorks at:"
    echo
    echo "http://www.mathworks.com/products/compiler/mcr/"
    echo 
    echo "If you already have Matlab 2019b with the Matlab Compiler or MCR 2019b"
    echo "installed, please edit ${exe_name} by uncommenting and editing the line:"
    echo "#BrainSuiteMCR=\"/path/to/your/MCR\";"
    echo "(replacing /path/to/your/MCR with the path to your Matlab or MCR installation)"
    echo "near the top of the file"
    echo
    exit 78
  fi
fi

read -d '' usage <<EOF

  tNLMPDFGOrdfMRI.sh : Performs tNLM PDF filtering
  Authored by Anand A. Joshi, Signal and Image Processing Institute
  Department of Electrical Engineering, Viterbi School of Engineering, USC
  
  Usage: tNLMPDFGOrdfMRI.sh GOrdInFile GOrdOutFile fpr memory MultiThreading

    GOrdInFile: Input grayordinate file. It is a .mat file.
    GOrdOutFile: Grayordinate output file (.mat)
    fpr: False positive rate parameter of the filter (.mat)
    memory: 'auto' or number in GB
    MutiThreading: 0 for off or 1 for off.
    scbpath: path of scb file for tNLMPDF filtering
EOF

# Parse inputs
if [ $# -lt 6 ]; then
  echo
  echo "$usage";
  echo
  exit;
fi

GOrdInFile=$1
GOrdOutFile=$2 
fpr=$3
memory=$4
MultiThreading=$5
scbpath=$6

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


# Run the executable
${exe_dir}/tNLMPDFGOrdfMRI "${GOrdInFile}" "${GOrdOutFile}" "${fpr}" "${memory}" "${MultiThreading}" "${scbpath}"

exit
