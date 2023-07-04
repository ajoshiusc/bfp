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

  usc_rigid_reg : performs rigid registration
  Authored by Anand A. Joshi, Signal and Image Processing Institute
  Department of Electrical Engineering, Viterbi School of Engineering, USC
  website: https://bitbucket.org/brainsuite/bfp/src/a851ada6d9006ca2b01fca38447e31caebbc1a7b/docs/?at=master
  
  Usage: transform_data_affine data_file data_coord output_file reg_moving_file reg_static_file reg_mat_file method

  e.g. transform_data_affine 'fmri_res.nii.gz' 'm' 'fmri_res2standard.nii.gz' 'example_func.nii.gz' 'standard.nii.gz' 'fmri_example_func2standard.rigid_registration_result.mat' 'linear'

EOF

# Parse inputs
if [ $# -lt 7 ]; then
  echo
  echo "$usage";
  echo
  exit;
fi

data_file=$1
data_coord=$2
output_file=$3
reg_moving_file=$4
reg_static_file=$5
reg_mat_file=$6
method=$7

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
${exe_dir}/transform_data_affine "${data_file}" "${data_coord}" "${output_file}" "${reg_moving_file}" "${reg_static_file}" "${reg_mat_file}" "${method}" 

exit
