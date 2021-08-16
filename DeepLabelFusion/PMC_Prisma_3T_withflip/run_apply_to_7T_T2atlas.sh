#!/bin/bash
#$ -S /bin/bash
set -e  -x

RUNCMD1="./run_crossval_update03062020.sh"
RUNCMD2="run_crossval_update03062020.sh"

export ASHS_ROOT=/data/picsl/longxie/pkg/ashs/ashs-fast-beta
export PATH=$PATH:$ASHS_ROOT/bin
MATLAB_BIN=/share/apps/matlab/R2013a/bin/matlab

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

# directories
ROOT=/home/longxie/DeepLabelFusion/PMC_Prisma_3T_withflip
ANALYSISDIR=$ROOT/analysis_input
EXPSETUPDIR=$ROOT/expsetup
EXPROOT=$ROOT/exp_ASHST1/
CODEDIR=$ROOT/scripts
#RAWDATADIR=/data/picsl/lwisse/NACC_Prisma_Atlas_NotFinal
RAWDATADIR=/data/picsl/pauly/wolk/abc_prisma/data
expid=exp_all
EXPDIR=$EXPROOT/exp${expid}
EXPTRAINDIR=$EXPDIR/atlas
EXPTESTDIR=$EXPDIR/testing_7TT2
EXPEVALDIR=$EXPDIR/evaluate_7TT2

##################################################
# parameters
#EVALLABELS=( SUB ERC BA35 BA36 PRC     PHC CS HIPPGDSC    HIPPDSC     EXPHIPP        ALLGDSC                ALLDSC)
#RANGES=(    1   2   4   "1 2 4" 3  7    8   10  11   12   "11 12" 13  14 "1 2 3 4 8" "1 2 3 4 8" "10 11 12 13" "1 2 3 4 8 10 11 12 13" "1 2 3 4 8 10 11 12 13")

EVALLABELS=(Hippo        CA1 CA2 CA3 CA      DG SUB Tail ERC BA35 BA36 PHC HippoSulcus CS MISC1 MISC2)
RANGES=(    "1 2 3 4 8"  1   2   4   "1 2 4" 3  8   5    9   10   11   12  6           14 7     13)

##################################################
function main()
{
  # 1. perform training
  #PrepareTrain
  #TrainASHS

  # 2. testing
  #Testing

  # 3. remove the intermediate files
  # only keep the segmentations and posteriors
  #CleanDir

  # 4. compute overlap
  # don`t run evaluation for now
  #Evaluate
}

##################################################
function PrepareTrain()
{
    mkdir -p $EXPTRAINDIR

    # prepare manifest file
    rm -rf $EXPTRAINDIR/manifest.txt
    for id in $(cat $ANALYSISDIR/subjlist_withflip.txt); do

      #for side in left right; do

      #  RAWSEG=$ROOT/ABC_3T_Atlas_final/${id}_${side}_finalLW.nii.gz
      #  if [[ ! -f $RAWSEG ]]; then
      #    echo "$id $side GT seg does not exist"
      #  fi
      #OUTSEG=$RAWDATADIR/$id/${id}_${side}_finalLW.nii.gz
      #  cp $RAWSEG $OUTSEG
      #done

      id="$id"
      MPRAGE=$(ls $RAWDATADIR/$id/*mprage.nii.gz)
      TSE=$(ls $RAWDATADIR/$id/*tse.nii.gz)
      LSEG=$(ls $RAWDATADIR/$id/*_left_finalLW.nii.gz)
      RSEG=$(ls $RAWDATADIR/$id/*_right_finalLW.nii.gz)
      echo "$id \
            $MPRAGE \
            $TSE \
            $LSEG \
            $RSEG" \
            >> $EXPTRAINDIR/manifest.txt

    done
}

##################################################
function TrainASHS()
{
  # run training
  if [[ 1 == 1 ]]; then
    time $ASHS_ROOT/bin/ashs_train.sh \
      -D $EXPTRAINDIR/manifest.txt \
      -L $ANALYSISDIR/Labels_3T_Atlas.txt \
      -w $EXPTRAINDIR \
      -C $ANALYSISDIR/ashs_config_vacind_multimodal.sh \
      -N \
      -s 1-7 \
      -Q -z $CODEDIR/ashs-fast-z.sh \
      | tee -a $EXPTRAINDIR/ashs_train.out.txt

      #-x $EXPTRAINDIR/xval.txt  \

  fi
}


