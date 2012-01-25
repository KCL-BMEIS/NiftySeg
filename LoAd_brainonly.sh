#!/bin/sh

# This is set as the minimum number of arguments for the script to run
ndefargs=2

# Parameters go here, where optionally, you can specify a default value
max_iter=100
bc_order=4

# Check number of args
image=$1
mask=$2

# Check mandatory args
if [ ! -f $image ]; then
    echo "The image '$image' to be segmented does not exist!"
    exit
fi

if [ ! -f $mask ]; then
    echo "The mask image '$mask' does not exist!" 
    exit
fi

# Parse remaining command line options after the mandatory ones.
shift $ndefargs
while [ "$#" -gt 0 ]
do
  case $1 in
    -max_iter)
      max_iter=$2
      shift 1
      ;;
    -bc_order)
      bc_order=$2
      shift 1
      ;;
    *)
      echo "Error: option $1 not recognised"
      exit
      ;;
  esac
  shift 1
done

# Create command line args for binary.
if [ "$max_iter" != "" ]; then
  max_iter_arg=" -max_iter $max_iter "
fi
if [ "$bc_order" != "" ]; then
  bc_order_arg=" -bc_order $bc_order "
fi

name=`echo $image | awk -F '.nii' {'print $1'}`
pathatlas=~/Source/nifty-seg/priors
atlas=${pathatlas}/T1.nii.gz
grey_prior=${pathatlas}/CGM_prior.nii.gz 
white_prior=${pathatlas}/WM_prior.nii.gz
csf_prior=${pathatlas}/ECSF_prior.nii.gz
deep_grey_prior=${pathatlas}/DGM_prior.nii.gz
internal_csf_prior=${pathatlas}/ICSF_prior.nii.gz

seg_maths ${mask} -bin -dil 3 -ero 3 -fill -dill ${name}__mask_filled.nii

if [ ! -f ${name}__affine.mat ]; then
  echo Starting reg_aladin
  linearArgs="-target ${image} -tmask ${name}__mask_filled.nii -source ${atlas} -aff ${name}__affine.mat -ln 3 -lp 2 -result ${name}__affine_result.nii.gz"
  LinearExecutable=`which reg_aladin`
  $LinearExecutable $linearArgs
fi

if [ ! -f ${name}__cpp.nii ]; then
  echo Starting reg_f3d
  nonLinearExecutable=`which reg_f3d`
  nonLinearArgs=" -ln 3 -lp 2 -sx 5 -target ${image} -tmask ${name}__mask_filled.nii -source ${atlas} -aff ${name}__affine.mat -cpp ${name}__cpp.nii"
  $nonLinearExecutable $nonLinearArgs
  echo Finished reg_f3d

  ###############################################################
  # Transform images from atlas space to native T1 space
  ################################################################
  reg_resample -target ${image} -source $grey_prior -cpp ${name}__cpp.nii -result ${name}__registered_grey.nii -TRI
  reg_resample -target ${image} -source $white_prior -cpp ${name}__cpp.nii -result ${name}__registered_white.nii -TRI
  reg_resample -target ${image} -source $csf_prior -cpp ${name}__cpp.nii -result ${name}__registered_csf.nii -TRI
  reg_resample -target ${image} -source $deep_grey_prior -cpp ${name}__cpp.nii -result ${name}__registered_deep_grey.nii -TRI
  reg_resample -target ${image} -source $internal_csf_prior -cpp ${name}__cpp.nii -result ${name}__registered_internal_csf.nii -TRI
  echo Finishing Resampling
fi

################################################################
# Run LoAd
################################################################

cmd="seg_LoAd -in ${image} -mask ${name}__mask_filled.nii -priors ${name}__registered_white.nii ${name}__registered_grey.nii ${name}__registered_csf.nii ${name}__registered_deep_grey.nii ${name}__registered_internal_csf.nii -out ${name}_segmentation.nii -v 1 $max_iter_arg $bc_order_arg"
echo "Running command:$cmd"
eval $cmd 

rm ${name}__*
