#!/usr/bin/env bash
echo $@

##########################################################################################################################
## SCRIPT TO PREPROCESS THE FUNCTIONAL SCAN
## parameters are passed from 0_preprocess.sh
##
## Written by the Underpants Gnomes (a.k.a Clare Kelly, Zarrar Shehzad, Maarten Mennes & Michael Milham)
## for more information see www.nitrc.org/projects/fcon_1000
##
##########################################################################################################################

## anatomical image
#t1=$1
t1=$1
#/home/ajoshi/coding_ground/bfp/data/sub-01-run1/anat/t1
## name of the fMRI scan
fmri=$2
#/home/ajoshi/coding_ground/bfp/data/sub-01-run1/func-1/fmri
## fmri scan directory
func_dir=$3
#/home/ajoshi/coding_ground/bfp/data/sub-01-run1/func-1
## first timepoint (remember timepoint numbering starts from 0)
#TRstart=$4
#0
## last timepoint
#TRend=$5
#450
## TR
TR=$4
#2
nuisance_template=$5
#/home/ajoshi/coding_ground/bfp/src/nuisance.fsf
## number of time points
#n_vols=$8
#451
n_vols=$(3dinfo -nv ${fmri}.nii.gz)
TRstart=0
TRend=$((n_vols-1))
echo $TR 
## set your desired spatial smoothing FWHM - we use 6 (acquisition voxel size is 3x3x4mm)
FWHM=$9
#6
sigma=`echo "scale=10 ; ${FWHM}/2.3548" | bc`

## Set high pass and low pass cutoffs for temporal filtering
hp=$10
#0.005
lp=$11
#0.1

## Example func image

example=$(basename "$fmri")_example

## directory setup

##########################################################################################################################
##---START OF SCRIPT----------------------------------------------------------------------------------------------------##
##########################################################################################################################

echo ---------------------------------------
echo !!!! PREPROCESSING FUNCTIONAL SCAN !!!!
echo ---------------------------------------

cwd=$( pwd )
cd ${func_dir}

## 1. Dropping first # TRS
echo "Dropping first TRs"
3dcalc -a ${fmri}.nii.gz[${TRstart}..${TRend}] -expr 'a' -prefix ${fmri}_dr.nii.gz

##2. Deoblique
echo "Deobliquing"
3drefit -deoblique ${fmri}_dr.nii.gz

##3. Reorient into fsl friendly space (what AFNI calls RPI)
echo "Reorienting"
3dresample -orient RPI -inset ${fmri}_dr.nii.gz -prefix ${fmri}_ro.nii.gz

##4. Motion correct to average of timeseries
echo "Motion correcting"
3dTstat -mean -prefix ${fmri}_ro_mean.nii.gz ${fmri}_ro.nii.gz 
3dvolreg -Fourier -twopass -base ${fmri}_ro_mean.nii.gz -zpad 4 -prefix ${fmri}_mc.nii.gz -1Dfile ${fmri}_mc.1D ${fmri}_ro.nii.gz

##5. Remove skull/edge detect
echo "Skull stripping"
3dAutomask -prefix ${fmri}_mask.nii.gz -dilate 1 ${fmri}_mc.nii.gz
3dcalc -a ${fmri}_mc.nii.gz -b ${fmri}_mask.nii.gz -expr 'a*b' -prefix ${fmri}_ss.nii.gz

##6. Get eighth image for use in registration
echo "Getting example_func for registration"
3dcalc -a ${fmri}_ss.nii.gz[7] -expr 'a' -prefix ${example}_func.nii.gz

##7. Spatial smoothing
echo "Smoothing"
fslmaths ${fmri}_ss.nii.gz -kernel gauss ${sigma} -fmean -mas ${fmri}_mask.nii.gz ${fmri}_sm.nii.gz

##8. Grandmean scaling
echo "Grand-mean scaling"
fslmaths ${fmri}_sm.nii.gz -ing 10000 ${fmri}_gms.nii.gz -odt float

##9. Temporal filtering
echo "Band-pass filtering"
3dFourier -lowpass ${lp} -highpass ${hp} -retrend -prefix ${fmri}_filt.nii.gz ${fmri}_gms.nii.gz

##10.Detrending
echo "Removing linear and quadratic trends"
3dTstat -mean -prefix ${fmri}_filt_mean.nii.gz ${fmri}_filt.nii.gz
3dDetrend -polort 2 -prefix ${fmri}_dt.nii.gz ${fmri}_filt.nii.gz
3dcalc -a ${fmri}_filt_mean.nii.gz -b ${fmri}_dt.nii.gz -expr 'a+b' -prefix ${fmri}_pp.nii.gz

##11.Create Mask
echo "Generating mask of preprocessed data"
fslmaths ${fmri}_pp.nii.gz -Tmin -bin ${fmri}_pp_mask.nii.gz -odt char


## 12.FUNC->T1
## You may want to change some of the options
flirt -ref ${t1}.bfc.nii.gz -in ${example}_func.nii.gz -out ${example}_func2t1.nii.gz -omat ${example}_func2t1.mat -cost corratio -dof 12 -interp trilinear
# Create mat file for conversion from subject's anatomical to functional
convert_xfm -inverse -omat t12${example}_func.mat ${example}_func2t1.mat
echo t12${example}_func.mat 

## TBD

## 12.FUNC->standard (3mm)
## You may want to change some of the options
flirt -ref standard.nii.gz -in ${example}_func.nii.gz -out ${example}_func2standard.nii.gz -omat ${example}_func2standard.mat -cost corratio -dof 12 -interp trilinear
# Create mat file for conversion from subject's anatomical to functional
convert_xfm -inverse -omat standard2${example}_func.mat ${example}_func2standard.mat
## TBD

## apply registration
#flirt -ref standard -in ${example}_func -out ${example}_func2standard -applyxfm -init ${example}_func2standard.mat -interp trilinear



## 13. 
nuisance_dir=${func_dir}/nuisance_$(basename "$fmri")

echo --------------------------------------------
echo !!!! RUNNING NUISANCE SIGNAL REGRESSION !!!!
echo --------------------------------------------


## 14. make nuisance directory
mkdir -p ${nuisance_dir}

# 15. Seperate motion parameters into seperate files
echo "Splitting up ${subject} motion parameters"
awk '{print $1}' ${fmri}_mc.1D > ${nuisance_dir}/mc1.1D
awk '{print $2}' ${fmri}_mc.1D > ${nuisance_dir}/mc2.1D
awk '{print $3}' ${fmri}_mc.1D > ${nuisance_dir}/mc3.1D
awk '{print $4}' ${fmri}_mc.1D > ${nuisance_dir}/mc4.1D
awk '{print $5}' ${fmri}_mc.1D > ${nuisance_dir}/mc5.1D
awk '{print $6}' ${fmri}_mc.1D > ${nuisance_dir}/mc6.1D

# Extract signal for global, csf, and wm
## 16. Global
echo "Extracting global signal for ${subject}"
3dmaskave -mask ${fmri}_pp_mask.nii.gz -quiet ${fmri}_pp.nii.gz > ${nuisance_dir}/global.1D

## 17. csf
flirt -ref ${example}_func.nii.gz -in ${t1}.pvc.label.nii.gz -out ${t1}.func.pvc.label.nii.gz -applyxfm -init t12${example}_func.mat -interp nearestneighbour

fslmaths ${t1}.func.pvc.label.nii.gz -thr 5.5 -bin ${t1}.func.csf.mask.nii.gz
fslmaths ${t1}.func.pvc.label.nii.gz -thr 2.5 -uthr 3.5 -bin ${t1}.func.wm.mask.nii.gz

echo "Extracting signal from csf"
3dmaskave -mask ${t1}.func.csf.mask.nii.gz -quiet ${fmri}_pp.nii.gz > ${nuisance_dir}/csf.1D

## 18. wm
echo "Extracting signal from white matter for ${subject}"
3dmaskave -mask ${t1}.func.wm.mask.nii.gz -quiet ${fmri}_pp.nii.gz > ${nuisance_dir}/wm.1D

## 6. Generate mat file (for use later)
## create fsf file

echo "Modifying model file"
sed -e s:nuisance_dir:"${nuisance_dir}":g <${nuisance_template} >${nuisance_dir}/temp1
sed -e s:nuisance_model_outputdir:"${nuisance_dir}/residuals.feat":g <${nuisance_dir}/temp1 >${nuisance_dir}/temp2
sed -e s:nuisance_model_TR:"${TR}":g <${nuisance_dir}/temp2 >${nuisance_dir}/temp3
sed -e s:nuisance_model_numTRs:"${n_vols}":g <${nuisance_dir}/temp3 >${nuisance_dir}/temp4
sed -e s:nuisance_model_input_data:"${func_dir}/${fmri}_pp.nii.gz":g <${nuisance_dir}/temp4 >${nuisance_dir}/nuisance.fsf 

#rm ${nuisance_dir}/temp*

echo "Running feat model"
feat_model ${nuisance_dir}/nuisance

minVal=`3dBrickStat -min -mask ${fmri}_pp_mask.nii.gz ${fmri}_pp.nii.gz`

## 7. Get residuals
echo "Running film to get residuals"
film_gls --rn=${nuisance_dir}/stats --noest --sa --ms=5 --in=${fmri}_pp.nii.gz --pd=${nuisance_dir}/nuisance.mat --thr=${minVal}

## 8. Demeaning residuals and ADDING 100
3dTstat -mean -prefix ${nuisance_dir}/stats/res4d_mean.nii.gz ${nuisance_dir}/stats/res4d.nii.gz
3dcalc -a ${nuisance_dir}/stats/res4d.nii.gz -b ${nuisance_dir}/stats/res4d_mean.nii.gz -expr '(a-b)+100' -prefix ${fmri}_res.nii.gz

## 9. Resampling residuals to MNI space
flirt -ref ${func_dir}/standard -in ${fmri}_res -out ${fmri}_res2standard -applyxfm -init ${func_dir}/${example}_func2standard.mat -interp trilinear

