#!/bin/bash
#$ -S /bin/bash
set -e  -x


export ASHS_ROOT=/data/picsl/longxie/pkg/ashs-fast-modified-nobl
export PATH=$PATH:$ASHS_ROOT/bin
MATLAB_BIN=/share/apps/matlab/R2016a/bin/matlab

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

# directories
ROOT=/data/picsl/longxie/DeepLabelFusion/PMC_Prisma_3T_withflip
CODEDIR=$ROOT/matlabcode
ANALYSISDIR=$ROOT/analysis_input
EXPSETUPDIR=$ROOT/expsetup
EXPROOT=$ROOT/exp_ASHST1/

# fold info
fold=0
expid=$(printf %03d $fold)
traintxt=$EXPSETUPDIR/train${expid}.txt
testtxt=$EXPSETUPDIR/test${expid}.txt
EXPDIR=$EXPROOT/exp${expid}
EXPTRAINDIR=$EXPDIR/atlas
EXPTESTDIR=$EXPDIR/testing
EXPRESULTDIR=$EXPDIR/result
DeepLFDIR=$EXPDIR/DeepLFdataset
DeepLF7TT2DIR=$EXPDIR/DeepLF7TT2dataset
DUMPDIR=$EXPROOT/dump
mkdir -p $DUMPDIR
SUBJTXT7TT2=/home/longxie/ASHS_7TT2/analysis_input/7TAtlasIDs.csv
EXPTEST7TT2DIR=$EXPDIR/testing_7TT2
EXPEVAL7TT2DIR=$EXPDIR/evaluate_7TT2
EXPRESULT7TT2DIR=$EXPDIR/result_7TT2


##################################################
# parameters
pr=7
pr_sample=13

##################################################
function main()
{
  reset_dir

  # link cross validation data
  #LinkData

  #######################
  # prepare for UNET and new pipeline
  PrepareDeepLFUnet

  #######################
  # prepare for the 7T T2 dataset
  #PrepareDeepLF7TT2

}

##################################################
function LinkData()
{
  # process training
  IDS=($(cat $traintxt))
  N=${#IDS[*]}

  # go through all the subjects
  PREFIX=LD
  for ((i=0;i<${N};i++)); do

    # submit jobs for evaluation
    id=${IDS[i]}
    qsub -cwd -o $DUMPDIR -j y -V \
      -q all.q \
      -l h_vmem=4.1G,s_vmem=4G \
      -N "${PREFIX}_${i}_${id}" \
      $0 LinkData_sub $i $id training
    sleep 0.1

  done

  # process test data
  IDS=($(cat $testtxt))
  N=${#IDS[*]}

  # go through all the subjects
  PREFIX=LD
  for ((i=0;i<${N};i++)); do

    # submit jobs for evaluation
    id=${IDS[i]}
    qsub -cwd -o $DUMPDIR -j y -V \
      -q all.q \
      -l h_vmem=4.1G,s_vmem=4G \
      -N "${PREFIX}_${i}_${id}" \
      $0 LinkData_sub $i $id test
    sleep 0.1

  done


  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function LinkData_sub()
{
  i=$1
  ii=$(printf %04d $((i+1)))
  id=$2
  type=$3

  # directories
  if [[ $type == "training" ]]; then
    SUBJTESTDIR=$EXPTRAINDIR/xval/xval${ii}/test/xval${ii}_test_${id}
  else
    SUBJTESTDIR=$EXPTESTDIR/$id/
  fi

  # reslice ground truth segmentation
  mkdir -p $SUBJTESTDIR/multiatlas/DeepLF
  mkdir -p $SUBJTESTDIR/bootstrap/DeepLF
  for side in left right; do
    c3d $SUBJTESTDIR/tse_native_chunk_${side}.nii.gz \
      $SUBJTESTDIR/refseg/refseg_${side}.nii.gz \
      -int 0 -reslice-identity \
      -o $SUBJTESTDIR/refseg/refseg_${side}_chunk.nii.gz

    # derive a mask
    c3d $SUBJTESTDIR/multiatlas/tseg_${side}_train*/atlas_to_native_segvote.nii.gz \
      -accum -add -endaccum \
      -binarize -dilate 1 5x5x2 \
      -o $SUBJTESTDIR/multiatlas/DeepLF/mask_${side}.nii.gz

    c3d $SUBJTESTDIR/bootstrap/tseg_${side}_train*/atlas_to_native_segvote.nii.gz \
      -accum -add -endaccum \
      -binarize -dilate 1 5x5x2 \
      -o $SUBJTESTDIR/bootstrap/DeepLF/mask_${side}.nii.gz

    # compute auto manual error map
    c3d $SUBJTESTDIR/refseg/refseg_${side}_chunk.nii.gz \
      -replace 6 0 7 0 13 0 14 0 \
      $SUBJTESTDIR/multiatlas/fusion/lfseg_heur_${side}.nii.gz \
      -replace 6 0 7 0 13 0 14 0 \
      -scale -1 -add -binarize \
      -o $SUBJTESTDIR/multiatlas/DeepLF/AutoSegErrorMap_${side}.nii.gz

    c3d $SUBJTESTDIR/refseg/refseg_${side}_chunk.nii.gz \
      -replace 6 0 7 0 13 0 \
      $SUBJTESTDIR/bootstrap/fusion/lfseg_heur_${side}.nii.gz \
      -replace 6 0 7 0 13 0 \
      -scale -1 -add -binarize \
      -o $SUBJTESTDIR/bootstrap/DeepLF/AutoSegErrorMap_${side}.nii.gz

  done
}

##################################################
function PrepareDeepLFUnet()
{
  mkdir -p $DeepLFDIR

  # go through all the subject
  PREFIX=DLFDS
  for type in training test; do

  if [[ $type == "training" ]]; then
    IDS=($(cat $traintxt | awk -F, '{print $1}'))
  else
    IDS=($(cat $testtxt | awk -F, '{print $1}'))
  fi
  N=${#IDS[*]}
  #N=1

  for ((i=0;i<${N};i++)); do
  #for ((i=0;i<1;i++)); do

  ii=$((i/2))
  ii=$((ii+1))
  ii=$(printf %04d $ii)
  id=${IDS[i]}

  # directories
  if [[ $type == "training" ]]; then
    SUBJTESTDIR=$EXPTRAINDIR/xval/xval${ii}/test/xval${ii}_test_${id}
    Natlas=$((N-2))
  else
    SUBJTESTDIR=$EXPTESTDIR/$id/
    Natlas=$N
  fi
  OUTDIR=$DeepLFDIR/$type/$id
  rm -rf $OUTDIR
  mkdir -p $OUTDIR
  cp -r $SUBJTESTDIR/multiatlas $OUTDIR/

    for side in left right; do

    # submit jobs for evaluation
    id=${IDS[i]}
    qsub -cwd -o $DUMPDIR -j y -V \
      -q all.q \
      -l h_vmem=4.1G,s_vmem=4G \
      -N "${PREFIX}_${type}_${i}_${id}_${side}" \
      $0 PrepareDeepLFUnet_sub $i $id $side $type
    sleep 0.1
 
  done
  done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function PrepareDeepLFUnet_sub()
{
  i=$1
  id=$2
  side=$3
  ii=$((i/2))
  ii=$((ii+1))
  ii=$(printf %04d $ii)
  type=$4
  IDS=($(cat $traintxt))
  N=${#IDS[*]}

  # directories
  if [[ $type == "training" ]]; then
    SUBJTESTDIR=$EXPTRAINDIR/xval/xval${ii}/test/xval${ii}_test_${id}
  else
    SUBJTESTDIR=$EXPTESTDIR/$id/
  fi
  REFSEGDIR=$SUBJTESTDIR/refseg
  Timg1FN=$SUBJTESTDIR/tse_native_chunk_${side}.nii.gz
  Timg2FN=$SUBJTESTDIR/mprage_to_tse_native_chunk_${side}.nii.gz
  TsegFN=$REFSEGDIR/refseg_${side}_chunk.nii.gz

  if [[ ! -f $TsegFN ]]; then
    c3d $Timg1FN $REFSEGDIR/refseg_${side}.nii.gz \
      -int 0 -reslice-identity \
      -o $TsegFN
  fi

  # output dir
  OUTDIR=$DeepLFDIR/$type/$id
  mkdir -p $OUTDIR

  # copy files
  cp $Timg1FN $OUTDIR/
  cp $Timg2FN $OUTDIR/
  cp $TsegFN $OUTDIR/

  # resample
  fslreorient2std \
    $Timg1FN \
    $TMPDIR/tse_native_chunk_${side}.nii.gz
  #c3d $TMPDIR/tse_native_chunk_${side}.nii.gz \
  #  -int 0 -resample-mm 0.4x0.4x0.4mm \
  #  -o $OUTDIR/tse_native_chunk_${side}_resampled.nii.gz
  # resample T1 to T2 native chunk space
  #greedy -d 3 \
  #  -rf $OUTDIR/tse_native_chunk_${side}_resampled.nii.gz \
  #  -ri NN \
  #  -rm $SUBJTESTDIR/mprage.nii.gz \
  #      $OUTDIR/mprage_to_tse_native_chunk_${side}_resampled.nii.gz \
  #  -r $SUBJTESTDIR/flirt_t2_to_t1/greedy_t2_to_t1_chunk_${side}.mat \
  #     $SUBJTESTDIR/flirt_t2_to_t1/flirt_t2_to_t1.mat
  #c3d $OUTDIR/tse_native_chunk_${side}_resampled.nii.gz \
  #  $Timg2FN \
  #  -int 0 -reslice-identity \
  #  -o $OUTDIR/mprage_to_tse_native_chunk_${side}_resampled.nii.gz
  #c3d $OUTDIR/tse_native_chunk_${side}_resampled.nii.gz \
  #  $TsegFN \
  #  -int 0 -reslice-identity \
  #  -o $OUTDIR/refseg_${side}_chunk_resampled.nii.gz

  echo "${id}_${side}_target,$OUTDIR/tse_native_chunk_${side}_resampled.nii.gz,$OUTDIR/refseg_${side}_chunk_resampled.nii.gz,Control" > $OUTDIR/$(printf %03d $i)_${type}_${id}_${side}_info.csv

  NATLASES=($(ls $OUTDIR/multiatlas/ | grep tseg_${side}_train | wc -l))

  for ((j=0;j<$NATLASES;j++)); do
    ATLASDIR=$OUTDIR/multiatlas/tseg_${side}_train$(printf %03d $j)
    #c3d $OUTDIR/tse_native_chunk_${side}_resampled.nii.gz -as TSE \
    #  $ATLASDIR/atlas_to_native.nii.gz \
    #  -int 0 -reslice-identity \
    #  -o $ATLASDIR/atlas_to_native_resampled.nii.gz \
    #  -push TSE \
    #  $ATLASDIR/atlas_to_native_mprage.nii.gz \
    #  -int 0 -reslice-identity \
    #  -o $ATLASDIR/atlas_to_native_mprage_resampled.nii.gz \
    #  -push TSE \
    #  $ATLASDIR/atlas_to_native_segvote.nii.gz \
    #  -int 0 -reslice-identity \
    #  -o $ATLASDIR/atlas_to_native_segvote_resampled.nii.gz

      echo "${id}_${side}_atlas$(printf %03d $j),$ATLASDIR/atlas_to_native_resampled.nii.gz,$OUTDIR/atlas_to_native_segvote_resampled.nii.gz,Control" >> $OUTDIR/$(printf %03d $i)_${type}_${id}_${side}_info.csv
  done 
}

##################################################
function PrepareDeepLF7TT2()
{
  mkdir -p $DeepLF7TT2DIR

  # go through all the subject
  PREFIX=DLFDS7TT2
  N=$(cat $SUBJTXT7TT2 | wc -l)

  for ((i=2;i<${N};i++)); do
 
  aid=$(cat -A $SUBJTXT7TT2 | head -n $i | tail -n 1 | cut -d, -f 1)
  id=$(echo $aid | cut -d _ -f 1)
  scandate=$(echo $aid | cut -d _ -f 2)
  OUTDIR=$DeepLF7TT2DIR/test/$id/
  rm -rf $OUTDIR
  mkdir -p $OUTDIR
  cp -r $EXPTEST7TT2DIR/$id/multiatlas $OUTDIR

  for side in left right; do

    # submit jobs for evaluation
    qsub -cwd -o $DUMPDIR -j y -V \
      -q all.q \
      -l h_vmem=4.1G,s_vmem=4G \
      -N "${PREFIX}_test_${id}_${side}" \
      $0 PrepareDeepLF7TT2_sub $id $side $EXPTEST7TT2DIR/$id $OUTDIR
    sleep 0.1

  done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function PrepareDeepLF7TT2_sub()
{
  id=$1
  side=$2
  SUBJTESTDIR=$3
  OUTDIR=$4
  REFSEGDIR=$SUBJTESTDIR/refseg
  Timg1FN=$SUBJTESTDIR/tse_native_chunk_${side}.nii.gz
  Timg2FN=$SUBJTESTDIR/mprage_to_tse_native_chunk_${side}.nii.gz
  TsegFN=$REFSEGDIR/refseg_${side}_chunk.nii.gz

  if [[ ! -f $TsegFN ]]; then
    c3d $Timg1FN $REFSEGDIR/refseg_${side}.nii.gz \
      -int 0 -reslice-identity \
      -o $TsegFN
  fi

  # copy files
  cp $Timg1FN $OUTDIR/
  cp $Timg2FN $OUTDIR/
  cp $TsegFN $OUTDIR/

  # resample
  #fslreorient2std \
  #  $Timg1FN \
  #  $TMPDIR/tse_native_chunk_${side}.nii.gz
  c3d $Timg1FN -swapdim RSA \
    -int 0 -resample-mm 0.4x0.4x0.4mm \
    -o $OUTDIR/tse_native_chunk_${side}_resampled.nii.gz
  c3d $OUTDIR/tse_native_chunk_${side}_resampled.nii.gz \
    $Timg2FN \
    -int 0 -reslice-identity \
    -o $OUTDIR/mprage_to_tse_native_chunk_${side}_resampled.nii.gz
  c3d $OUTDIR/tse_native_chunk_${side}_resampled.nii.gz \
    $TsegFN \
    -int 0 -reslice-identity \
    -o $OUTDIR/refseg_${side}_chunk_resampled.nii.gz

  NATLASES=($(ls $OUTDIR/multiatlas/ | grep tseg_${side}_train | wc -l))

  for ((j=0;j<$NATLASES;j++)); do

    ATLASDIR=$OUTDIR/multiatlas/tseg_${side}_train$(printf %03d $j)
    c3d $OUTDIR/tse_native_chunk_${side}_resampled.nii.gz -as TSE \
      $ATLASDIR/atlas_to_native.nii.gz \
      -int 0 -reslice-identity \
      -o $ATLASDIR/atlas_to_native_resampled.nii.gz \
      -push TSE \
      $ATLASDIR/atlas_to_native_mprage.nii.gz \
      -int 0 -reslice-identity \
      -o $ATLASDIR/atlas_to_native_mprage_resampled.nii.gz \
      -push TSE \
      $ATLASDIR/atlas_to_native_segvote.nii.gz \
      -int 0 -reslice-identity \
      -o $ATLASDIR/atlas_to_native_segvote_resampled.nii.gz

  done
}


##################################################
function reset_dir()
{
  rm -rf $DUMPDIR/*
  mkdir -p $DUMPDIR
}

##################################################
# Main entrypoint
cmd=$0
if [[ $# -lt 2 ]]; then

  main $@

else
  cmd=$1
  shift
  $cmd $@
fi

