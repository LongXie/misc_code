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
export PATH=$ANTSPATH:$C3DPATH:$PATH

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

# Directories
ROOT=/data/jet/longxie/ASHS_PHC
WORKDIR=$ROOT/thickness_newlabel
ASHSRUNDIR=$ROOT/SuperResolution/ASHSexp/exp002/fullset/ashs
MATLABCODEDIR=$WORKDIR/matlabcode/
GTGROUPDIR=$WORKDIR/group/
SUBJ_TXT=$WORKDIR/analysis_input/subj.txt


########################################
# 0. copy data
LABEL_IDS_ALL=(BKG      CA      DG  SUB ERC  BA35 BA36 PHC  CS   OTS)
LABEL_MRG_ALL=("0 7 16" "1 2 4" "3" "8" "10" "11" "12" "13" "14" "16")
# Labels to get segmentation
LABEL_IDS=(BKG      CA      DG  SUB ERC  BA35 BA36 PHC  CS)
LABEL_MRG=("0 7 16" "1 2 4" "3" "8" "10" "11" "12" "13" "14")
LABEL_NEW=(0        1       2   3   4    5    6    7    8)
KINDS="tse mprage ${LABEL_IDS[*]}"

expid=004
expdir=$WORKDIR/SCexp/exp${expid}
FINALGRAPHDIR=$expdir/clustering/finalgraph/
FINALTEMPDIR=$expdir/FinalTemp
# 4.1
ALLPATHSDIR=$FINALTEMPDIR/paths/allpaths
# 4.3
BESTPATHSDIR=$FINALTEMPDIR/paths/bestpaths
DATADIR=${expdir}/data
PWDIR=${expdir}/pairwise/

groups=(1 2 3)

ITER=8
ANTs_start=3
WGT=1
REG_LABELS=(${LABEL_IDS[*]})
ANTs_t="SyN[0.25]"
ANTs_r="Gauss[0.5,0]"
ANTs_i="80x80x20"
ANTs_x="Y"
ANTs_all_metrics_2="--use-all-metrics-for-convergence"

TESTDIR=/home/longxie/ASHS_PHC/thickness_newlabel/test/
mkdir -p $TESTDIR
rm -rf $TESTDIR/*
#id=DW209_T00_20091022
side=left
#grp=1
idx_mov=1
idx_fix=84


IDS=($(cat $SUBJ_TXT))
id_mov=${IDS[$idx_mov]}

# get the path
paths=($(cat $BESTPATHSDIR/$side/${id_mov}_${side}.txt))
N=${#IDS[*]}
N_grps=${#groups[*]}
N_all=${#paths[*]}
#idx=$((N_all-N_grps-1+grp_fix))
cur_path=${paths[$idx_fix]}

# generate IDS array with templates
for ((i=0;i<$N;i++)); do
  id_tmp=${IDS[i]}
  IDS[$i]=$(ls $ASHSRUNDIR | grep $id_tmp)
done

for ((i=0;i<$N_grps;i++)); do
  k=$((i+1))
  IDS[$((N+i))]="template${k}"
done

id_mov=${IDS[$idx_mov]}
id_fix=${IDS[$idx_fix]}

# loop
i=1
warps=""
invwarps=""
while true; do

  # from subject
  idx_from=$(echo $cur_path | cut -d , -f $i)
  id_from=${IDS[$((idx_from-1))]}

   # to subject
  set +e
  i=$((i+1))
  idx_to=$(echo $cur_path | cut -d , -f $i)
  id_to=${IDS[$((idx_to-1))]}
  set -e

  if [[ $idx_to == "" ]]; then
    break
  fi

  # get warps
  warps="$PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_pairwiseWarp.nii.gz $PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_mlaffine_itk.txt $warps"
  
  invwarps="$invwarps -i $PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_mlaffine_itk.txt $PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_pairwiseInverseWarp.nii.gz"

  done

# compose transform
$ANTSPATH/ComposeMultiTransform 3 \
  $TESTDIR/${id_mov}_${side}_InitWarp.nii.gz \
  -R $DATADIR/${id_fix}_${side}_tse.nii.gz \
  $warps

$ANTSPATH/ComposeMultiTransform 3 \
  $TESTDIR/${id_mov}_${side}_InitInverseWarp.nii.gz \
  -R $DATADIR/${id_mov}_${side}_tse.nii.gz \
  $invwarps

# apply init warps
for sub in $KINDS; do

    WarpImageMultiTransform 3 \
      $DATADIR/${id_mov}_${side}_${sub}.nii.gz \
      $TESTDIR/${id_mov}_${side}_inittotemp_reslice_${sub}.nii.gz \
      -R $DATADIR/${id_fix}_${side}_tse.nii.gz \
      $TESTDIR/${id_mov}_${side}_InitWarp.nii.gz

done

# Create the segmentation for the template
c3d $(for sub in ${LABEL_IDS[*]}; do echo $TESTDIR/${id_mov}_${side}_inittotemp_reslice_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $TESTDIR/${id_mov}_${side}_inittotemp_reslice_seg.nii.gz





cp $DATADIR/${id_mov}_${side}_*.nii.gz $TESTDIR/
cp $DATADIR/${id_fix}_${side}_*.nii.gz $TESTDIR/




