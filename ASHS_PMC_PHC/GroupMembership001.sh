#!/bin/bash
#$ -S /bin/bash
set -x -e

########################################################
# INIT
#ANTSPATH=/home/avants/bin/ants
ANTSPATH=/data/picsl/pauly/bin/ants/
NEWANTSDIR=/share/apps/ANTs/2014-06-23/build/bin/
C3DPATH=/data/picsl/longxie/pkg/c3d_tool/bin
PKGDIR=/data/picsl/longxie/pkg/
MATLAB_BIN=/share/apps/matlab/R2016a/bin/matlab
export PATH=$ANTSPATH:$C3DPATH:$PATH

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

# Directories
ROOT=/data/jet/longxie/ASHS_PHC
WORKDIR=$ROOT/thickness_newlabel
MATLABCODEDIR=$WORKDIR/matlabcode/
GTGROUPDIR=$WORKDIR/group/
SUBJ_TXT=$WORKDIR/analysis_input/subj.txt
SUBJ_WITHTEMP_TXT=$WORKDIR/analysis_input/subj_withtemp.txt

########################################################
# Parameters needs to be specify
# Experiment number
expid=001
expdir=$WORKDIR/exp/exp_member${expid}
DUMPDIR=${expdir}/dump
mkdir -p ${expdir}/dump

########################################
# 1.1 copy data
LABEL_IDS_ALL=(BKG                ERC  BA35 BA36 PHC  aCS  pCS)
LABEL_MRG_ALL=("0 1 2 3 4 7 8 16" "10" "11" "12" "13" "15" "14")
# Labels to get segmentation
LABEL_IDS=(BKG                ERC  BA35 BA36 PHC  aCS  pCS)
LABEL_MRG=("0 1 2 3 4 7 8 16" "10" "11" "12" "13" "15" "14")
LABEL_NEW=(0                  1    2    3    4    5    6)
MESH_LABEL=(ERC  BA35 BA36 PHC)
KINDS="tse ${LABEL_IDS[*]}"
DATADIR=${expdir}/data
ASHSRUNDIR=$ROOT/OriginASHS/fullset/ashs

# 1.2 group membership 
WGT_1=1
REG_LABELS_1=(${LABEL_IDS[*]})
ANTs_t_1="SyN[0.25]"
ANTs_r_1="Gauss[1.5,0.5]"
ANTs_i_1="15x8x0"
ANTs_x_1="Y"
ANTs_G_1=""
ANTs_all_metrics_1="--use-all-metrics-for-convergence"
groups=(1 2 3)
GROUPDIR=$expdir/GroupMembership
MTDIR=/home/longxie/ASHS_PHC/thickness_newlabel/exp/exp556/MultiTemps

#######################################################
function main()
{
  reset_dir

  ################################
  # preparation
  # 1.1 copy data
  #copy_data

  # 1.2 decide which group membersip
  GroupMembership

}

############################################################
function copy_data()
{
  # Get the data
  mkdir -p $DATADIR

  # Load id number
  IDs=($(cat $SUBJ_TXT))

  # Submit job to copy data
  PREFIX=CP${expid}
  for side in left right; do
    for ((i=0;i<${#IDs[*]};i++)); do

      id=${IDs[i]}
      fn=$(ls $ASHSRUNDIR | grep $id)

      qsub -cwd -o $DUMPDIR -j y \
           -l h_vmem=2.1G,s_vmem=2G \
           -N "${PREFIX}_${id}_${side}" \
           $0 copy_subject $fn $side
      sleep 0.1

     done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function copy_subject()
{
  fn=$1
  side=$2

  # ASHS segmentation
  SEG=$ASHSRUNDIR/${fn}/final/${fn}_${side}_lfseg_corr_usegray_cleanup_dividedCS.nii.gz

  # Link the subfield images
  if [[ -f $SEG ]]; then
    ln -sf \
       $ASHSRUNDIR/$fn/tse_native_chunk_${side}.nii.gz \
       $DATADIR/${fn}_${side}_tse.nii.gz
  fi

  # Generate a binary image for the label
  for ((i=0;i<${#LABEL_IDS_ALL[*]};i++)); do

    c3d $DATADIR/${fn}_${side}_tse.nii.gz \
      $SEG -replace $(for k in ${LABEL_MRG_ALL[i]}; do echo $k 999; done) -thresh 999 999 1 0 \
      -int 0 -reslice-identity \
      -o $DATADIR/${fn}_${side}_${LABEL_IDS_ALL[i]}.nii.gz

  done

  # Vote in the ASHS template subject space to create segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $DATADIR/${fn}_${side}_${sub}.nii.gz; done) \
    -vote -type ushort -o $DATADIR/${fn}_${side}_seg.nii.gz
}

############################################################
function GroupMembership()
{
  # Get the data
  rm -rf $GROUPDIR/tmp
  mkdir -p $GROUPDIR/tmp

  # Load id number
  IDs=($(cat $SUBJ_TXT))
  #IDs="1030_T00"

  # Submit job to copy data
  PREFIX=GM${expid}
  for side in left right; do
    for ((i=0;i<${#IDs[*]};i++)); do

      id=${IDs[i]}
      fn=$(ls $ASHSRUNDIR | grep $id)

      qsub -cwd -o $DUMPDIR -j y \
           -l h_vmem=2.1G,s_vmem=2G \
           -N "${PREFIX}_${id}_${side}" \
           $0 GroupMembership_sub $fn $side $GROUPDIR/tmp
      sleep 0.1

     done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # all group and similarity info
  rm -rf $GROUPDIR/group_${side}.txt \
         $GROUPDIR/OVL_${side}.txt \
         $GROUPDIR/IDs_*.txt
  for side in left right; do
    for ((i=0;i<${#IDs[*]};i++)); do

      id=${IDs[i]}
      fn=$(ls $ASHSRUNDIR | grep $id)

      cat $GROUPDIR/tmp/${fn}_${side}_group.txt \
        >> $GROUPDIR/group_${side}.txt
      cat $GROUPDIR/tmp/${fn}_${side}_OVL.txt \
       >> $GROUPDIR/OVL_${side}.txt

      # save IDs
      grp=$(cat $GROUPDIR/tmp/${fn}_${side}_group.txt)
      echo $id >> $GROUPDIR/IDs_${side}_${grp}.txt

    done
  done

  # remove tmp dir
  rm -rf $GROUPDIR/tmp
}

function GroupMembership_sub()
{
  fn=$1
  side=$2
  GPTMPDIR=$3
  OVL_max=0
  OVL_all=""

  for ((i=0;i<${#groups[*]};i++)); do

    grp=${groups[i]}

    ########################
    # registration
    # Use ml_affine for nice affine alignment
    /data/picsl/pauly/wolk/ashs/ext/Linux/bin/ml_affine \
      $DATADIR/${fn}_${side}_seg.nii.gz \
      $MTDIR/group${grp}/work/iter_${side}_07/template_${side}_${grp}_seg.nii.gz \
      $GPTMPDIR/template${grp}_to_${fn}_${side}_mlaffine.txt

    # Convert that to ITK format
    c3d_affine_tool \
      $GPTMPDIR/template${grp}_to_${fn}_${side}_mlaffine.txt \
      -oitk $GPTMPDIR/template${grp}_to_${fn}_${side}_mlaffine_itk.txt

    CMD=""
    for sub in ${REG_LABELS_1[*]}; do
      CMD="$CMD -m MSQ[$DATADIR/${fn}_${side}_${sub}.nii.gz,$MTDIR/group${grp}/work/iter_${side}_07/template_${side}_${grp}_${sub}.nii.gz,$WGT_1]"
    done

    if [[ $ANTs_x_1 == "Y" ]]; then
      c3d $DATADIR/${fn}_${side}_seg.nii.gz \
        -binarize -dilate 1 15x15x15vox \
        -o $GPTMPDIR/template${grp}_to_${fn}_${side}_mask.nii.gz
      ANTs_mask_1="-x $GPTMPDIR/template${grp}_to_${fn}_${side}_mask.nii.gz"
    else
      ANTs_mask_1=""
    fi

    # Perform ANTs registration
    ANTS 3 $CMD \
      -t $ANTs_t_1 \
      -r $ANTs_r_1 \
      -i $ANTs_i_1 \
      $ANTs_G_1 \
      $ANTs_mask_1 \
      -a $GPTMPDIR/template${grp}_to_${fn}_${side}_mlaffine_itk.txt \
      --continue-affine 0 \
      $ANTs_all_metrics_1 \
      -o $GPTMPDIR/template${grp}_to_${fn}_${side}_group.nii.gz \
      | tee $GPTMPDIR/template${grp}_to_${fn}_${side}_antsoutput.txt

    for sub in ${KINDS[*]}; do

      WarpImageMultiTransform 3 \
        $MTDIR/group${grp}/work/iter_${side}_07/template_${side}_${grp}_${sub}.nii.gz \
        $GPTMPDIR/template${grp}_to_${fn}_${side}_reslice_${sub}.nii.gz \
        -R $DATADIR/${fn}_${side}_tse.nii.gz \
        $GPTMPDIR/template${grp}_to_${fn}_${side}_groupWarp.nii.gz \
        $GPTMPDIR/template${grp}_to_${fn}_${side}_groupAffine.txt

    done

    # Create seg
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $GPTMPDIR/template${grp}_to_${fn}_${side}_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $GPTMPDIR/template${grp}_to_${fn}_${side}_reslice_seg.nii.gz

    for sub in ${KINDS[*]}; do

      rm -rf $GPTMPDIR/template${grp}_to_${fn}_${side}_reslice_${sub}.nii.gz

    done

    ###################################
    # measure similarity
    OVL_tmp=$(c3d \
      $DATADIR/${fn}_${side}_seg.nii.gz \
      -replace 1 0 4 0 6 0 \
      $GPTMPDIR/template${grp}_to_${fn}_${side}_reslice_seg.nii.gz \
      -replace 1 0 4 0 6 0 \
      -label-overlap \
      | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}')

    if [[ $OVL_tmp > $OVL_max ]]; then

      group_final=$grp
      OVL_max=$OVL_tmp

    fi

    OVL_all="$OVL_all $OVL_tmp"

  done

  ####################################
  # write result to disk
  echo $OVL_all > $GPTMPDIR/${fn}_${side}_OVL.txt
  echo $group_final > $GPTMPDIR/${fn}_${side}_group.txt
}

############################################################
function reset_dir()
{
  rm -rf $DUMPDIR/*
}

############################################################
if [[ $# -lt 2 ]]; then

  main

else

  cmd=$1
  shift
  $cmd $@

fi
