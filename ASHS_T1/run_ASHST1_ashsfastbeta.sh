#/bin/bash
#$ -S /bin/bash
set -e -x

#set -x -e

##############################################
# Setup environment
# Software PATH
MTLATLAS=/home/longxie/ASHS_T1/ASHSexp/exp201/atlas/final_fast_ashs_beta
ICVATLAS=/data/jet/vpiskin/ASHS/build_ICV_atlas_trimmed/ICV_atlas/final/
export ASHS_ROOT=/data/picsl/longxie/pkg/ashs/ashs-fast-beta
export PATH=$PATH:$ASHS_ROOT/bin
C3DPATH=$ASHS_ROOT/ext/Linux/bin
#MATLAB_BIN=/share/apps/matlab/R2017a/bin/matlab
#SRMATLABCODEDIR=/home/longxie/ASHS_PHC/SuperResolution/SRToolBox/matlabfunction
#SRPATH=/home/longxie/pkg/PatchSuperResolution/build_release
#CODEDIR=/home/longxie/ASHS_T1/pipeline_package


# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

infile=$1
OUTDIR=$2
DELETETMP=$3

#############################################
function main()
{

  filenameall=$(basename "$infile")
  extension="${filenameall##*.}"
  filename="${filenameall%.*}"

  if [[ $extension == "gz" ]]; then
    extension="${filename##*.}"
    if [[ $extension == "nii" ]]; then
      filename="${filename%.*}"
    fi
  fi

  # make output directory
  mkdir -p $OUTDIR
  if [[ $infile != $OUTDIR/${filename}.nii.gz ]]; then
    cp -f $infile $OUTDIR/${filename}.nii.gz
  fi


  # perform MTL segmentation
  echo "Step 1/3: Performing medial temporal lobe segmentation"
  MTLSeg
  echo "Step 1/3: Done!"
  
  # perform ICV segmentation
  echo "Step 2/3: Performing intracranial volume segmentation"
  ICVSeg
  echo "Step 2/3: Done!"

  # reorganize and summarize the result
  echo "Step 3/3: Reorganize output and summarize the result"
  Summarize
  echo "Step 3/3: Done!"
 

}

#############################################
function MTLSeg()
{
  $ASHS_ROOT/bin/ashs_main.sh \
    -a $MTLATLAS -d -T -I $filename -g $infile \
    -f $infile \
    -s 1-7 \
    -m $ASHS_ROOT/bin/identity.mat -M \
    -w $OUTDIR/MTLSeg
}

#############################################
function ICVSeg()
{
  $ASHS_ROOT/bin/ashs_main.sh \
    -a $ICVATLAS -d -T -I $filename -g $infile \
    -f $infile \
    -s 1-7 \
    -B \
    -w $OUTDIR/ICVSeg
}

#############################################
function Summarize()
{
  cp $OUTDIR/MTLSeg/tse.nii.gz \
     $OUTDIR/${filename}_denoised_SR.nii.gz

  cp $OUTDIR/MTLSeg/final/${filename}_left_lfseg_heur.nii.gz \
     $OUTDIR/${filename}_MTLSeg_left.nii.gz
  cp $OUTDIR/MTLSeg/final/${filename}_right_lfseg_heur.nii.gz \
     $OUTDIR/${filename}_MTLSeg_right.nii.gz

  cp $OUTDIR/MTLSeg/final/${filename}_left_heur_volumes.txt \
     $OUTDIR/${filename}_MTLSeg_left_volumes.txt
  cp $OUTDIR/MTLSeg/final/${filename}_right_heur_volumes.txt \
     $OUTDIR/${filename}_MTLSeg_right_volumes.txt

  cp $OUTDIR/ICVSeg/final/${filename}_left_lfseg_corr_nogray.nii.gz \
     $OUTDIR/${filename}_ICVSeg.nii.gz

  cp $OUTDIR/ICVSeg/final/${filename}_left_corr_nogray_volumes.txt \
     $OUTDIR/${filename}_ICVSeg_volumes.txt

  mkdir $OUTDIR/qc
  cp -r $OUTDIR/MTLSeg/qc $OUTDIR/qc/ASHST1_qc
  cp -r $OUTDIR/ICVSeg/qc $OUTDIR/qc/ICV_qc


  if [[ $DELETETMP != "0" ]]; then  
    rm -rf $OUTDIR/MTLSeg $OUTDIR/ICVSeg
  fi
}

######################################################

main


