#!/bin/bash
#$ -S /bin/bash
set -x -e

##############################################################################
# INIT
#ANTSPATH=/home/avants/bin/ants
ANTSPATH=/data/picsl/pauly/bin/ants/
NEWANTSDIR=/share/apps/ANTs/2014-06-23/build/bin/
C3DPATH=/data/picsl/longxie/pkg/c3d_tool/bin
PKGDIR=/data/picsl/longxie/pkg/
MATLAB_BIN=/share/apps/matlab/R2013a/bin/matlab
export PATH=$ANTSPATH:$C3DPATH:$PATH

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

# Directories
ROOT=/data/jet/longxie/ASHS_PHC
WORKDIR=$ROOT/thickness_newlabel
ASHSRUNDIR=$ROOT/SuperResolution/ASHSexp/exp003/fullset/ashs
MATLABCODEDIR=$WORKDIR/matlabcode/
GTGROUPDIR=$WORKDIR/group/
SUBJ_TXT=$WORKDIR/analysis_input/subj.txt
SUBJ_WITHTEMP_TXT=$WORKDIR/analysis_input/subj_withtemp.txt

###############################################################
# parameters
LABEL_IDS=(CA      DG  SUB HIPPO       ERC  BA35 BA36 PHC )
LABEL_MRG=("1 2 4" "3" "8" "1 2 3 4 8" "10" "11" "12" "13")

ROOT=/home/longxie/ASHS_PHC
WORKDIR=$ROOT/thickness_newlabel
SREXP="exp003"
ASHSRUNDIR=$ROOT/SuperResolution/ASHSexp/$SREXP/fullset/ashs
VOLUMEDIR=$WORKDIR/volume
mkdir -p $VOLUMEDIR

SUBJ_TXT=$WORKDIR/analysis_input/subj.txt

# Load id number
IDs=($(cat $SUBJ_TXT))

# generate header
STRING="ID"
for side in left right; do
  for label in ${LABEL_IDS[*]}; do
    STRING="$STRING,${label}_${side}"
  done
done
echo $STRING > $VOLUMEDIR/SR${SREXP}_volume.txt

# Submit job to copy data
for ((i=0;i<${#IDs[*]};i++)); do

  id=${IDs[i]}
  fn=$(ls $ASHSRUNDIR | grep $id)
  volume_all=${id}

  for side in left right; do

    SEG=$ASHSRUNDIR/${fn}/final/${fn}_${side}_lfseg_corr_usegray_dividedCS.nii.gz

    for ((j=0;j<${#LABEL_IDS[*]};j++)); do
      label=${LABEL_IDS[j]}
      volume=$(c3d $SEG -replace $(for k in ${LABEL_MRG[j]}; do echo $k 999; done) -thresh 999 999 1 0 -dup -lstat | head -n 3 | tail -n 1 | awk '{print $7}')
      volume_all="$volume_all,$volume"
    done

  done

  echo $volume_all >> $VOLUMEDIR/SR${SREXP}_volume.txt
done






