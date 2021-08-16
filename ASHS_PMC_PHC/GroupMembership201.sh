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
ANALYSISDIR=$WORKDIR/analysis_input
MATLABCODEDIR=$WORKDIR/matlabcode/
GTGROUPDIR=$WORKDIR/group/
SUBJ_TXT=$WORKDIR/analysis_input/subj.txt
SUBJ_WITHTEMP_TXT=$WORKDIR/analysis_input/subj_withtemp.txt
FLIPXFN=$WORKDIR/analysis_input/flipx_itk.txt
IDENTITYFN=$WORKDIR/analysis_input/identity_itk.txt

########################################################
# Parameters needs to be specify
# Experiment number
expid=201
expdir=$WORKDIR/exp/exp_member${expid}
WEIGHTDIR=$expdir/weights
DUMPDIR=${expdir}/dump
mkdir -p ${expdir}/dump
GROUPDIR=$expdir/GroupMembership

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

#######################################################
function main()
{
  reset_dir

  EvalWeight

  GroupMembership 
}

##################################################
function EvalWeight()
{
  # Get Target IDS
  IDS=$(cat $SUBJ_TXT)
  #IDS="DW104"
  mkdir -p $WEIGHTDIR

  # Loop through each subjects and do label fusion
  PREFIX=EW${expid}
  for id in $IDS; do
    for side in left right; do

      # Submit job to do label fusion
      qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -N "${PREFIX}_${id}_${side}" \
         $0 EvalWeight_sub $id $side
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # organize file for PRC
  for side in left right; do
    rm -rf $WEIGHTDIR/weights_PRC_${side}.txt
    rm -rf $WEIGHTDIR/weights_ERC_${side}.txt

    for allid in $IDS; do
      OVERLAP=""
      OVERLAPERC=""
      idx=1
      allfn=$(ls $ASHSRUNDIR | grep $allid)

      for xvalid in $(cat $ANALYSISDIR/subj_xval.txt); do

        #xvalfn=$(ls $FULLSETDIR/ashs | grep $xvalid)

        if [[ $allid == $xvalid ]]; then
          OVERLAP="$OVERLAP -1"
          OVERLAPERC="$OVERLAPERC -1"
        else
          SUBJDIR=$WEIGHTDIR/$allfn/$side
          OVERLAP="$OVERLAP $(cat $SUBJDIR/weights_PRC_${side}.txt | head -n $idx | tail -n 1)"
          OVERLAPERC="$OVERLAPERC $(cat $SUBJDIR/weights_ERC_${side}.txt | head -n $idx | tail -n 1)"
          idx=$((idx+1))
        fi

      done

      # output weight
      echo $OVERLAP >> $WEIGHTDIR/weights_PRC_${side}.txt
      echo $OVERLAPERC >> $WEIGHTDIR/weights_ERC_${side}.txt

    done
  done
}

function EvalWeight_sub()
{
  id=$1
  id=$(ls $ASHSRUNDIR | grep $id)
  side=$2
  SUBJDIR=$ASHSRUNDIR/$id

  # outputdir
  OUTDIR=$WEIGHTDIR/$id/$side
  mkdir -p $OUTDIR

  # perform label fusion
  /data/picsl-build/pauly/malf/gcc64rel/label_fusion 3 \
    -g $SUBJDIR/bootstrap/tseg_${side}_train*/atlas_to_native.nii.gz \
    -l $SUBJDIR/bootstrap/tseg_${side}_train*/atlas_to_native_segvote.nii.gz \
    -m Joint[0.1,2] \
    -rp 3x3x1 \
    -rs 0x0x0 \
    -w $OUTDIR/weight_${side}_%03d.nii.gz \
    $SUBJDIR/tse_native_chunk_${side}.nii.gz \
    $OUTDIR/lfseg_raw_${side}.nii.gz

  # weights
  SEG=$ASHSRUNDIR/$id/final/${id}_${side}_lfseg_corr_usegray_cleanup_dividedCS.nii.gz
  c3d $OUTDIR/lfseg_raw_${side}.nii.gz \
    $SEG -int 0 -reslice-identity \
    -o $OUTDIR/lfseg_corr_usegray_cleanup_dividedCS_${side}.nii.gz
  rm -rf $OUTDIR/weights_${side}.txt
  rm -rf $OUTDIR/weights_PRC_${side}.txt
  for w in $(ls $OUTDIR | grep weight | grep -v weights); do

    OVERLAP=$(c3d $OUTDIR/$w $OUTDIR/lfseg_corr_usegray_cleanup_dividedCS_${side}.nii.gz -lstat | awk '{print $2}' | awk '{printf("%s ",$1)}' |awk '{printf("%.05f %.05f %.05f %.05f %.05f %.05f %.05f %.05f %.05f %.05f %.05f %.05f",$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14)}'| awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}')

    echo $OVERLAP >> $OUTDIR/weights_${side}.txt

    OVERLAP_PRC=$(c3d $OUTDIR/$w $OUTDIR/lfseg_corr_usegray_cleanup_dividedCS_${side}.nii.gz -replace 11 999 12 999 15 999 -thresh 999 999 1 0 -lstat | awk '{print $2}' | awk '{printf("%s ",$1)}' | awk '{printf("%.05f", $3)}' | awk '{print $1}')

    echo $OVERLAP_PRC >> $OUTDIR/weights_PRC_${side}.txt

    OVERLAP_ERC=$(c3d $OUTDIR/$w $OUTDIR/lfseg_corr_usegray_cleanup_dividedCS_${side}.nii.gz -replace 10 999 -thresh 999 999 1 0 -lstat | awk '{print $2}' | awk '{printf("%s ",$1)}' | awk '{printf("%.05f", $3)}' | awk '{print $1}')

    echo $OVERLAP_ERC >> $OUTDIR/weights_ERC_${side}.txt

  done
}

############################################################
function GroupMembership()
{
  # Get the data
  mkdir -p $GROUPDIR
  rm -rf $GROUPDIR/tmp
  mkdir -p $GROUPDIR/tmp

  # Load id number
  IDs=($(cat $SUBJ_TXT))

  # Submit job to copy data
  PREFIX=GM${expid}
  for side in left right; do

    qsub -cwd -o $DUMPDIR -j y \
         -l h_vmem=2.1G,s_vmem=2G \
         -N "${PREFIX}_${side}" \
         $0 GroupMembership_sub $side
    sleep 0.1

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function GroupMembership_sub()
{
  side=$1

  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('/data/picsl/longxie/ASHS/thickness_newlabel/matlabcode');
    GroupMembershipLabelFusion('$WEIGHTDIR/weights_PRC_${side}.txt','$WORKDIR/group/group_xval_${side}.txt','$GROUPDIR/group_${side}.txt');
MATCODE
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







