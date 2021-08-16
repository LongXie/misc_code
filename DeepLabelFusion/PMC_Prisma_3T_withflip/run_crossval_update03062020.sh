#!/bin/bash
#$ -S /bin/bash
#set -e
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
RAW7TDATADIR=/home/longxie/ASHS_7TT2/rawdata/ABC_7T_Atlas_final
SUBJTXT7TT2=/home/longxie/ASHS_7TT2/analysis_input/7TAtlasIDs.csv

##################################################
# parameters
#EVALLABELS=( SUB ERC BA35 BA36 PRC     PHC CS HIPPGDSC    HIPPDSC     EXPHIPP        ALLGDSC                ALLDSC)
#RANGES=(    1   2   4   "1 2 4" 3  7    8   10  11   12   "11 12" 13  14 "1 2 3 4 8" "1 2 3 4 8" "10 11 12 13" "1 2 3 4 8 10 11 12 13" "1 2 3 4 8 10 11 12 13")

EVALLABELS=(Hippo          CA1 CA2 CA3 CA      DG SUB Tail ERC BA35 BA36 PHC HippoSulcus CS MISC1 MISC2)
RANGES=(    "1 2 3 4 5 8"  1   2   4   "1 2 4" 3  8   5    9   10   11   12  6           14 7     13)

##################################################
function main()
{
  TRAINFNS=($(ls $EXPSETUPDIR | grep train | grep .txt))

  # loop over all the experiments
  #for ((i=0;i<${#TRAINFNS[*]};i++)); do
  #for ((i=1;i<59;i++)); do
  i=$1   

    trainfn=${TRAINFNS[i]}
    expid=$(echo $trainfn | sed -e "s/train//" | sed -e "s/.txt//")
    testfn="test${expid}.txt"
    traintxt=$EXPSETUPDIR/$trainfn
    testtxt=$EXPSETUPDIR/$testfn
    EXPDIR=$EXPROOT/exp${expid}
    EXPTRAINDIR=$EXPDIR/atlas
    EXPTESTDIR=$EXPDIR/testing
    EXPRESULTDIR=$EXPDIR/result
    EXPEVALDIR=$EXPDIR/evaluate
    EXPTEST7TT2DIR=$EXPDIR/testing_7TT2
    EXPEVAL7TT2DIR=$EXPDIR/evaluate_7TT2
    EXPRESULT7TT2DIR=$EXPDIR/result_7TT2

    EXPTEST7TT2SAMEDIR=$EXPDIR/testing_7TT2_same_subject
    EXPEVAL7TT2SAMEDIR=$EXPDIR/evaluate_7TT2_same_subject
    EXPRESULT7TT2SAMEDIR=$EXPDIR/result_7TT2_same_subject

    EXPVARYATLASDIR=$EXPDIR/VaryNumAtlas

    EXPTRUEVARYATLASDIR=$EXPDIR/TrueVaryNumAtlas

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
 
    # 5. compute entrophy
    #Entropy
 
    # 6. apply to 7T T2 atlas
    #Testing_7TT2

    # 7. reorganize files
    #CleanDir_7TT2

    # 8. compute overlap
    Evaluate_7TT2

    # 9 apply to 7T T2 of the same test subjects 
    #Testing_7TT2_samesubject

    # 10 run # atlas vs. performance experiment
    #NumatlasVSPerformance

    # 11 compute dice of all evaluation experiments
    #ComputeNumatlasVSPerformanceDice

    # 12 run true # atlas vs. performance experiment
    #TrainNumatlasVSPerformance

    # 13 test true # atlas vs. performance experiment
    #TestNumatlasVSPerformance

    # 14 evaluate true # atlas vs. performance experiment
    #ComputeTrueNumatlasVSPerformanceDice

  #done
}

##################################################
function PrepareTrain()
{
    mkdir -p $EXPTRAINDIR $EXPTESTDIR

    # prepare manifest file
    rm -rf $EXPTRAINDIR/manifest.txt
    rm -rf $EXPTRAINDIR/xval.txt
    for id in $(cat $traintxt | awk -F, '{print $1}'); do

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

    for id in $(cat $traintxt | awk -F, '{print $1}' | sed 's/orig_//g' | sed 's/flip_//g' | uniq); do
      echo "orig_${id} flip_${id}" >> $EXPTRAINDIR/xval.txt
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
      -x $EXPTRAINDIR/xval.txt \
      -Q -z $CODEDIR/ashs-fast-z.sh \
      | tee -a $EXPTRAINDIR/ashs_train.out.txt
 
      #-x $EXPTRAINDIR/xval.txt  \

  fi
}

##################################################
function Testing()
{
  mkdir -p $EXPTESTDIR/dump
  PREFIX="xval${i}"
  for id in $(cat $testtxt | awk -F, '{print $1}'); do

    id="$id"
    T1=$(ls $RAWDATADIR/$id/*mprage.nii.gz)
    T2=$(ls $RAWDATADIR/$id/*tse.nii.gz)
    LSEG=$(ls $RAWDATADIR/$id/*_left_finalLW.nii.gz)
    RSEG=$(ls $RAWDATADIR/$id/*_right_finalLW.nii.gz)

    #for side in left right; do

        #RAWSEG=$ROOT/ABC_3T_Atlas_final/${id}_${side}_finalLW.nii.gz
        #if [[ ! -f $RAWSEG ]]; then
        #  echo "$id $side GT seg does not exist"
        #fi
        #OUTSEG=$RAWDATADIR/$id/${id}_${side}_finalLW.nii.gz
        #cp $RAWSEG $OUTSEG
      #done

    # submit job
    #qsub -cwd -o $EXPTESTDIR/dump -j y -V \
    #  -q all.q,basic.q \
    #  -l h_vmem=8.1G,s_vmem=8G \
    #  -N "${PREFIX}_${id}" \
    $ASHS_ROOT/bin/ashs_main.sh \
        -a $EXPTRAINDIR/final \
        -d -T \
        -I $id -g $T1 -f $T2 \
        -Q -s 1-7 \
        -z $CODEDIR/ashs-fast-z-test.sh \
        -w $EXPTESTDIR/$id \
        -r "$LSEG $RSEG" &
    sleep 1

  done
}

##################################################
function CleanDir()
{
  # move the stuff we want to keep to the result directory
  mkdir -p $EXPRESULTDIR
  cp $traintxt $EXPRESULTDIR/
  cp $testtxt $EXPRESULTDIR/
  for id in $(cat $testtxt | awk -F, '{print $1}'); do
 
    id="$id"
    mkdir -p $EXPRESULTDIR/$id/multiatlas
    mkdir -p $EXPRESULTDIR/$id/bootstrap
    ln -sf $EXPTESTDIR/$id/final $EXPRESULTDIR/$id/final
    ln -sf $EXPTESTDIR/$id/multiatlas/fusion \
       $EXPRESULTDIR/$id/multiatlas/fusion 
    ln -sf $EXPTESTDIR/$id/bootstrap/fusion \
       $EXPRESULTDIR/$id/bootstrap/fusion

  done

  # move core dump files to result directory
  mkdir -p $EXPRESULTDIR/core
  set +e
  COREDUMP=($(ls $ROOT/code | grep core.))
  set -e
  if [ ${#COREDUMP[*]} -gt 0 ]; then
    mv  $ROOT/code/core.* $EXPRESULTDIR/core
  fi

  # remove training and testing
  #rm -rf $EXPTRAINDIR
  #rm -rf $EXPTESTDIR
}

##################################################
function Evaluate()
{
  mkdir -p $EXPRESULTDIR/dump
  PREFIX="eval${i}"
  for id in $(cat $testtxt | awk -F, '{print $1}'); do

    id="$id"
    T1=$(ls $RAWDATADIR/$id/*mprage.nii.gz)
    T2=$(ls $RAWDATADIR/$id/*tse.nii.gz)
    LSEG=$(ls $RAWDATADIR/$id/*_left_finalLW.nii.gz)
    RSEG=$(ls $RAWDATADIR/$id/*_right_finalLW.nii.gz)
    ln -sf $T1 $EXPRESULTDIR/$id/mprage.nii.gz
    ln -sf $T2 $EXPRESULTDIR/$id/tse.nii.gz
    ln -sf $LSEG $EXPRESULTDIR/$id/left_gtseg.nii.gz
    ln -sf $RSEG $EXPRESULTDIR/$id/right_gtseg.nii.gz

    # submit jobs for evaluation
    qsub -cwd -o $EXPRESULTDIR/dump -j y -V \
      -q all.q \
      -l h_vmem=8.1G,s_vmem=8G \
      -N "${PREFIX}_${id}" \
      $0 Evaluate_sub $id $EXPRESULTDIR
    sleep 0.1

  done

  # Wait for completion
  qsub -cwd -o $EXPRESULTDIR/dump -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  rm -rf $EXPEVALDIR
  mkdir -p $EXPEVALDIR
  for run in multiatlas bootstrap; do
  for kind in heur corr_nogray corr_usegray; do
  for id in $(cat $testtxt | awk -F, '{print $1}'); do
  for side in left right; do
    if [[ ! -f $EXPEVALDIR/overlap_${run}_${kind}.csv ]]; then
      cat $EXPRESULTDIR/$id/$run/fusion/overlap_${run}_${kind}_${side}.txt \
        > $EXPEVALDIR/overlap_${run}_${kind}.csv
    else
      cat $EXPRESULTDIR/$id/$run/fusion/overlap_${run}_${kind}_${side}.txt \
        | tail -n 1 >> $EXPEVALDIR/overlap_${run}_${kind}.csv
    fi
  done
  done
  done
  done
}


function Evaluate_sub()
{
  id=$1
  EXPRESULTDIR=$2

  for run in multiatlas bootstrap; do
  for side in left right; do
    for kind in heur corr_nogray corr_usegray; do

      AUTOSEG=$EXPRESULTDIR/$id/$run/fusion/lfseg_${kind}_${side}.nii.gz
      GTSEG=$EXPRESULTDIR/$id/${side}_gtseg.nii.gz
      do_pair  $AUTOSEG $GTSEG

      # get GDSC of all gray matter
      GMALL=$(c3d $AUTOSEG -dup $GTSEG \
        -int 0 -reslice-identity \
        -foreach -replace 6 0 7 0 13 0 14 0 -endfor \
        -label-overlap \
        | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}')

     # write output
     echo "ID,SIDE,GMALL,$(echo ${EVALLABELS[*]} | sed 's/ /,/g')" > \
       $EXPRESULTDIR/$id/$run/fusion/overlap_${run}_${kind}_${side}.txt
     echo "$id,$side,$GMALL$FULLOVL" >> \
       $EXPRESULTDIR/$id/$run/fusion/overlap_${run}_${kind}_${side}.txt

    done
  done
  done
}

function do_pair()
{
  # Get a pair of segmentations
  seg_a=$1
  seg_b=$2

  # out dice file
  #out_dice_file=$3

  # Iterate over all relevant labels
  FULLOVL=""
  for ((i=0; i<${#EVALLABELS[*]}; i++)); do

    # Do the analysis on full-size meshes
    REPRULE=$(for lab in ${RANGES[i]}; do echo $lab 99; done)

    # Extract the binary images and compute overlap
    c3d \
      $seg_a -dup $seg_b -int 0 -reslice-identity \
      -foreach -replace $REPRULE -thresh 99 99 1 0 -endfor \
      -overlap 1 | tee $TMPDIR/ovl.txt

    # Get the full-extent overlap
    OVL=$(cat $TMPDIR/ovl.txt | grep OVL | awk -F '[ ,]+' '{print $6}')

    #echo $id ${LABELS[i]} full $OVL $DIST >> $out_file
    FULLOVL="${FULLOVL},${OVL}"

  done

<<COMMENT
  # Iterate over all relevant labels
  FULLOVL=""
  for ((i=0; i<${#EVALLABELS[*]}; i++)); do

    # compute DSC or GDSC
    if [[ ${EVALTYPE[i]} == "D" ]]; then

      # Do the analysis on full-size meshes
      REPRULE=$(for lab in ${RANGES[i]}; do echo $lab 99; done)

      # Extract the binary images and compute overlap
      c3d \
        $seg_a -dup $seg_b -int 0 -reslice-identity \
        -foreach -replace $REPRULE -thresh 99 99 1 0 -endfor \
        -overlap 1 | tee $TMPDIR/ovl.txt

      # Get the full-extent overlap
      OVL=$(cat $TMPDIR/ovl.txt | grep OVL | awk -F '[ ,]+' '{print $6}')

    elif [[ ${EVALTYPE[i]} == "G" ]]; then

      ii=9981
      REPRULE=$(for lab in ${RANGES[i]}; \
        do echo $lab $ii; \
        ii=$((ii+1)); \
        done)

      OVL=$(c3d \
        $seg_a -dup $seg_b -int 0 -reslice-identity \
        -foreach -replace $REPRULE \
        -dup -thresh 9981 inf 1 0 -multiply -endfor \
        -label-overlap | awk '{print $3}' | \
        awk '{printf("%s ",$1)}' | \
        awk '{printf("%.04f ",$3)}' | awk '{print $1}')

    fi

    #echo $id ${LABELS[i]} full $OVL $DIST >> $out_file
    FULLOVL="$FULLOVL $OVL"

  done
COMMENT
}

##################################################
function Entropy()
{
  PREFIX="EP${i}"
  for id in $(cat $testtxt); do

    #id="pp$id"

    # submit jobs for evaluation
    qsub -cwd -o $EXPRESULTDIR/dump -j y -V \
      -q all.q \
      -l h_vmem=4.1G,s_vmem=4G \
      -N "${PREFIX}_${id}" \
      $0 Entropy_sub $id $EXPRESULTDIR
    sleep 0.1

  done

  # Wait for completion
  qsub -cwd -o $EXPRESULTDIR/dump -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function Entropy_sub()
{
  id=$1
  EXPRESULTDIR=$2

  for side in left right; do
    for kind in corr_usegray; do

      AUTOSEG=$EXPRESULTDIR/$id/posteriors/bootstrap/fusion/lfseg_${kind}_${side}.nii.gz
      POSTPREFIX=$EXPRESULTDIR/$id/posteriors/bootstrap/fusion/posterior_${kind}_${side}

      # Get the sum of the posterior images
      #PSUM=$(c3d ${POSTPREFIX}_* \
      #  -accum -add -endaccum -probe 50% | awk '{print $NF}')
      PSUM=$(c3d ${POSTPREFIX}_* \
         -accum -add -endaccum \
         -dup -thresh 0.0001 inf 1 0 -lstat \
         | tail -n 1 | awk '{print $4}')

      # Compute the entropy
      c3d ${POSTPREFIX}_* \
        -foreach -stretch 0 $PSUM 0 1 -clip 0 1 -dup \
        -scale 1 -log -clip -100000 100000 -times \
        -scale -1.44269504089 -endfor \
        -accum -add -endaccum \
        -o $TMPDIR/${id}_entropy_${side}.nii.gz

      c3d $TMPDIR/${id}_entropy_${side}.nii.gz \
        $AUTOSEG \
        -thresh 1 inf 1 0 -lstat | tail -n 1 | \
        awk '{print $2","$4}' > \
        $EXPRESULTDIR/$id/final/entropy_${side}_${kind}_total.txt


      # write output
      #echo "ID ${EVALLABELS[*]}" > \
      #  $EXPRESULTDIR/$id/final/entropy_${side}_${kind}.txt
      #echo "$id $FULLOVL" >> \
      #  $EXPRESULTDIR/$id/final/entropy_${side}_${kind}.txt

    done
  done
}

##################################################
function Testing_7TT2()
{
  mkdir -p $EXPTEST7TT2DIR/dump
  PREFIX="test_7TT2_${i}"
  N=$(cat $SUBJTXT7TT2 | wc -l)
  #N=2
  for ((j=2;j<=${N};j++)); do

    aid=$(cat -A $SUBJTXT7TT2 | head -n $j | tail -n 1 | cut -d, -f 1)
    id=$(echo $aid | cut -d _ -f 1)
    scandate=$(echo $aid | cut -d _ -f 2)
    SCANDATADIR=/data/jux/wolk_group/Prisma3T/relong
    T1=$(ls $SCANDATADIR/$id/MRI7T/$scandate/processed/*mp2rinv2.nii.gz)
    T2=$(ls $SCANDATADIR/$id/MRI7T/$scandate/processed/T2_*.nii.gz)
    LSEG=$(ls $RAW7TDATADIR/$aid/${id}_left_LW5.nii.gz)
    RSEG=$(ls $RAW7TDATADIR/$aid/${id}_right_LW5.nii.gz)

    if [[ ! -d $EXPTEST7TT2DIR/$id ]]; then

    qsub -cwd -o $EXPTEST7TT2DIR/dump -j y -V \
      -q all.q \
      -l h_vmem=10.1G,s_vmem=10G \
      -N "${PREFIX}_${id}" \
       $ASHS_ROOT/bin/ashs_main.sh \
         -a $EXPTRAINDIR/final \
         -d -T \
         -I $id -g $T1 -f $T2 \
         -s 1-7 \
         -C $ANALYSISDIR/ashs_config_vacind_multimodal_7TT2.sh \
         -w $EXPTEST7TT2DIR/$id \
         -r "$LSEG $RSEG" 
    sleep 1

    fi

    #-z $CODEDIR/ashs-fast-z-test7T.sh

  done
}

##################################################
function CleanDir_7TT2()
{
  # move the stuff we want to keep to the result directory
  mkdir -p $EXPRESULT7TT2DIR
  N=$(cat $SUBJTXT7TT2 | wc -l)
  #N=2
  for ((j=2;j<=${N};j++)); do

    aid=$(cat -A $SUBJTXT7TT2 | head -n $j | tail -n 1 | cut -d, -f 1)
    id=$(echo $aid | cut -d _ -f 1)
    scandate=$(echo $aid | cut -d _ -f 2)
    SCANDATADIR=/data/jux/wolk_group/Prisma3T/relong

    mkdir -p $EXPRESULT7TT2DIR/$id/multiatlas
    mkdir -p $EXPRESULT7TT2DIR/$id/bootstrap
    ln -sf $EXPTEST7TT2DIR/$id/final $EXPRESULT7TT2DIR/$id/final
    ln -sf $EXPTEST7TT2DIR/$id/multiatlas/fusion \
       $EXPRESULT7TT2DIR/$id/multiatlas/fusion
    ln -sf $EXPTEST7TT2DIR/$id/bootstrap/fusion \
       $EXPRESULT7TT2DIR/$id/bootstrap/fusion

  done

  #rm -rf $EXPTEST7TT2DIR
}

##################################################
function Evaluate_7TT2()
{
  mkdir -p $EXPRESULT7TT2DIR/dump
  PREFIX="eval7TT2${i}"
  N=$(cat $SUBJTXT7TT2 | wc -l)
  #N=2
  for ((j=2;j<=${N};j++)); do

    aid=$(cat -A $SUBJTXT7TT2 | head -n $j | tail -n 1 | cut -d, -f 1)
    id=$(echo $aid | cut -d _ -f 1)
    if [[ $id == "122812" ]]; then
      continue
    fi
    scandate=$(echo $aid | cut -d _ -f 2)
    SCANDATADIR=/data/jux/wolk_group/Prisma3T/relong
    T1=$(ls $SCANDATADIR/$id/MRI7T/$scandate/processed/*mp2rinv2.nii.gz)
    T2=$(ls $SCANDATADIR/$id/MRI7T/$scandate/processed/T2_*.nii.gz)
    LSEG=$(ls $RAW7TDATADIR/$aid/${id}_left_LW5.nii.gz)
    RSEG=$(ls $RAW7TDATADIR/$aid/${id}_right_LW5.nii.gz)
    ln -sf $T1 $EXPRESULT7TT2DIR/$id/mprage.nii.gz
    ln -sf $T2 $EXPRESULT7TT2DIR/$id/tse.nii.gz
    ln -sf $LSEG $EXPRESULT7TT2DIR/$id/left_gtseg.nii.gz
    ln -sf $RSEG $EXPRESULT7TT2DIR/$id/right_gtseg.nii.gz

    # submit jobs for evaluation
    qsub -cwd -o $EXPRESULT7TT2DIR/dump -j y -V \
      -q all.q \
      -l h_vmem=8.1G,s_vmem=8G \
      -N "${PREFIX}_${id}" \
      $0 Evaluate_7TT2_sub $id $EXPRESULT7TT2DIR
    sleep 0.1

  done

  # Wait for completion
  qsub -cwd -o $EXPRESULT7TT2DIR/dump -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  rm -rf $EXPEVAL7TT2DIR
  mkdir -p $EXPEVAL7TT2DIR
  for run in multiatlas bootstrap; do
  for kind in heur corr_nogray corr_usegray; do
  for ((j=2;j<=${N};j++)); do
  for side in left right; do

    aid=$(cat -A $SUBJTXT7TT2 | head -n $j | tail -n 1 | cut -d, -f 1)
    id=$(echo $aid | cut -d _ -f 1)
    if [[ $id == "122812" ]]; then
      continue
    fi

    if [[ ! -f $EXPEVAL7TT2DIR/overlap_${run}_${kind}.csv ]]; then
      cat $EXPRESULT7TT2DIR/$id/$run/fusion/overlap_${run}_${kind}_${side}.txt \
        > $EXPEVAL7TT2DIR/overlap_${run}_${kind}.csv
    else
      cat $EXPRESULT7TT2DIR/$id/$run/fusion/overlap_${run}_${kind}_${side}.txt \
        | tail -n 1 >> $EXPEVAL7TT2DIR/overlap_${run}_${kind}.csv
    fi
  done
  done
  done
  done
}

function Evaluate_7TT2_sub()
{
  id=$1
  EXPRESULT7TT2DIR=$2

  for run in multiatlas bootstrap; do
  for side in left right; do
    for kind in heur corr_nogray corr_usegray; do

      AUTOSEG=$EXPRESULT7TT2DIR/$id/$run/fusion/lfseg_${kind}_${side}.nii.gz
      GTSEG=$EXPRESULT7TT2DIR/$id/${side}_gtseg.nii.gz
      c3d $GTSEG -replace 14 7 15 6 16 14 -o $TMPDIR/${id}_${side}_gtseg.nii.gz
      GTSEG=$TMPDIR/${id}_${side}_gtseg.nii.gz
      do_pair  $AUTOSEG $GTSEG

      # get GDSC of all gray matter
      GMALL=$(c3d $AUTOSEG -dup $GTSEG \
        -int 0 -reslice-identity \
        -foreach -replace 6 0 7 0 13 0 14 0 -endfor \
        -label-overlap \
        | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}')

     # write output
     echo "ID,SIDE,GMALL,$(echo ${EVALLABELS[*]} | sed 's/ /,/g')" > \
       $EXPRESULT7TT2DIR/$id/$run/fusion/overlap_${run}_${kind}_${side}.txt
     echo "$id,$side,$GMALL$FULLOVL" >> \
       $EXPRESULT7TT2DIR/$id/$run/fusion/overlap_${run}_${kind}_${side}.txt

    done
  done
  done
}

##################################################
function Testing_7TT2_samesubject()
{
  mkdir -p $EXPTEST7TT2SAMEDIR/dump
  PREFIX="test_7TT2_same_${i}"
  N=$(cat $testtxt | wc -l)
  #N=2
  for ((j=1;j<=${N};j++)); do

    aid=$(cat -A $testtxt | head -n $j | tail -n 1 | cut -d, -f 1)
    fliptype=$(echo $aid | cut -d _ -f 1)
    id=$(echo $aid | cut -d _ -f 2)
    #scandate=$(echo $aid | cut -d _ -f 3)
    SCANDATADIR=/data/jux/wolk_group/Prisma3T/relong
    set +e
    T2=($(ls $SCANDATADIR/$id/MRI7T/*/processed/T2_*.nii.gz))
    T2=${T2[0]}
    scandate=$(echo $T2 | cut -d "/" -f 9)
    T1=$(ls $SCANDATADIR/$id/MRI7T/$scandate/processed/*mp2rinv2.nii.gz)
    set -e 
    #LSEG=$(ls $RAW7TDATADIR/$aid/${id}_left_LW5.nii.gz)
    #RSEG=$(ls $RAW7TDATADIR/$aid/${id}_right_LW5.nii.gz)

    if [[ -f $T2 && -f $T1 && ! -f $EXPTEST7TT2SAMEDIR/${fliptype}_${id}/final/${fliptype}_${id}_left_lfseg_heur.nii.gz ]]; then

    mkdir -p $EXPTEST7TT2SAMEDIR/${fliptype}_${id}
    if [[ $fliptype == "flip" ]]; then
      /data/picsl/pauly/bin/c3d_affine_tool -sform $T2 -tran 447 0 0 -scale -1 1 1 \
        -sform $T2 -inv -mult -mult -mult \
        -o $EXPTEST7TT2SAMEDIR/${fliptype}_${id}/${fliptype}_${id}.mat
      c3d $T1 -pad 30x30x30vox 30x30x30vox 0 \
        -dup -reslice-matrix $EXPTEST7TT2SAMEDIR/${fliptype}_$id/${fliptype}_${id}.mat \
        -as FT1 \
        -as FT1 \
        -thresh 0.0001 inf 1 0 -trim 0x0x0vox \
        -push FT1 -int 0 -reslice-identity \
        -o $EXPTEST7TT2SAMEDIR/${fliptype}_${id}/${fliptype}_${id}_T1.nii.gz
      T1=$EXPTEST7TT2SAMEDIR/${fliptype}_${id}/${fliptype}_${id}_T1.nii.gz
      c3d $T2 -dup -reslice-matrix $EXPTEST7TT2SAMEDIR/${fliptype}_$id/${fliptype}_${id}.mat \
        -o $EXPTEST7TT2SAMEDIR/${fliptype}_${id}/${fliptype}_${id}_T2.nii.gz
      T2=$EXPTEST7TT2SAMEDIR/${fliptype}_${id}/${fliptype}_${id}_T2.nii.gz
    else
      ln -sf $T1 $EXPTEST7TT2SAMEDIR/${fliptype}_${id}/${fliptype}_${id}_T1.nii.gz
      ln -sf $T2 $EXPTEST7TT2SAMEDIR/${fliptype}_${id}/${fliptype}_${id}_T2.nii.gz

    fi

    qsub -cwd -o $EXPTEST7TT2SAMEDIR/dump -j y -V \
      -q all.q \
      -l h_vmem=10.1G,s_vmem=10G \
      -N "${PREFIX}_${id}" \
       $ASHS_ROOT/bin/ashs_main.sh \
         -a $EXPTRAINDIR/final \
         -d -T \
         -I ${fliptype}_${id} -g $T1 -f $T2 \
         -s 1-7 \
         -C $ANALYSISDIR/ashs_config_vacind_multimodal_7TT2.sh \
         -w $EXPTEST7TT2SAMEDIR/${fliptype}_${id}
    sleep 1

    fi

    #-z $CODEDIR/ashs-fast-z-test7T.sh

  done
}

##################################################
function NumatlasVSPerformance()
{
  mkdir -p $EXPVARYATLASDIR

  # prepare manifest file
  rm -rf $EXPVARYATLASDIR/manifest.txt
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
            >> $EXPVARYATLASDIR/manifest.txt

    done

  # run training
  if [[ 1 == 1 ]]; then
    time $ASHS_ROOT/bin/ashs_train.sh \
      -D $EXPVARYATLASDIR/manifest.txt \
      -L $ANALYSISDIR/Labels_3T_Atlas.txt \
      -w $EXPVARYATLASDIR/ashsrun \
      -C $ANALYSISDIR/ashs_config_vacind_multimodal.sh \
      -N \
      -s 6-7 \
      -x $EXPSETUPDIR/vary_Natlases_xval_fold${expid}.txt \
      -Q -z $CODEDIR/ashs-fast-z.sh \
      -B \
      | tee -a $EXPVARYATLASDIR/ashs_train.out.txt
  fi
}

##################################################
function ComputeNumatlasVSPerformanceDice()
{
  rm -rf $EXPVARYATLASDIR/dump
  mkdir -p $EXPVARYATLASDIR/dump $EXPVARYATLASDIR/evaluate
  PREFIX="evalNvsP${i}"
  N=$(cat $EXPSETUPDIR/vary_Natlases_xval_fold${expid}.txt | wc -l)
  for ((i=1;i<=${N};i++)); do

    ALLTESTS=($(cat $EXPSETUPDIR/vary_Natlases_xval_fold${expid}.txt | head -n $i | tail -n 1))
    NTEST=${#ALLTESTS[*]}
    NATLAS=$((46-$NTEST))
    
    # submit jobs for evaluation
    qsub -cwd -o $EXPVARYATLASDIR/dump -j y -V \
      -q all.q \
      -l h_vmem=4.1G,s_vmem=4G \
      -N "${PREFIX}_${i}_${NATLAS}" \
      $0 ComputeNumatlasVSPerformanceDice_sub \
      $i $NATLAS $EXPVARYATLASDIR/evaluate $testtxt $EXPVARYATLASDIR
    sleep 0.1

  done

  # Wait for completion
  qsub -cwd -o $EXPVARYATLASDIR/dump -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function ComputeNumatlasVSPerformanceDice_sub()
{
  i=$1
  idx=$(printf %04d $i)
  NATLAS=$2
  OUTDIR=$3
  mkdir -p $OUTDIR
  testtxt=$4
  EXPVARYATLASDIR=$5
  OUTCSV=$OUTDIR/eval_test_numselatlas$(printf %03d $NATLAS)_xval${idx}.csv
  echo "ID,SIDE,NATLAS,REP,GMALL,Hippo,CA1,CA2,CA3,CA,DG,SUB,Tail,ERC,BA35,BA36,PHC,HippoSulcus,CS,MISC1,MISC2" > $OUTCSV


  for id in $(cat $testtxt | awk -F, '{print $1}'); do
    for side in left right; do

      AUTOSEG=$EXPVARYATLASDIR/ashsrun/xval/xval${idx}/test/xval${idx}_test_${id}/final/xval${idx}_test_${id}_${side}_lfseg_corr_usegray.nii.gz
      GTSEG=$(ls $RAWDATADIR/$id/*_${side}_finalLW.nii.gz)

      GMALLTMP=$(c3d $GTSEG $AUTOSEG \
        -foreach -replace 6 0 7 0 13 0 14 0 -endfor \
        -label-overlap \
        |  awk '{print $3}' | awk '{printf("%s ", $1)}'  | awk '{print $3}' )
      
      HIPPO=$(c3d $GTSEG $AUTOSEG \
        -foreach -replace 1 999 2 999 3 999 4 999 5 999 8 999 -thresh 999 999 1 0 -endfor  \
        -overlap 1 | cut -d , -f 5)

      CA=$(c3d $GTSEG $AUTOSEG \
        -foreach -replace 1 999 2 999 4 999 -thresh 999 999 1 0 -endfor  \
        -overlap 1 | cut -d , -f 5)

      HippoSulcus=$(c3d $GTSEG $AUTOSEG \
        -foreach -replace 6 999 -thresh 999 999 1 0 -endfor  \
        -overlap 1 | cut -d , -f 5)
      MISC1=$(c3d $GTSEG $AUTOSEG \
        -foreach -replace 7 999 -thresh 999 999 1 0 -endfor  \
        -overlap 1 | cut -d , -f 5)
      CS=$(c3d $GTSEG $AUTOSEG \
        -foreach -replace 13 999 -thresh 999 999 1 0 -endfor  \
        -overlap 1 | cut -d , -f 5)
      MISC2=$(c3d $GTSEG $AUTOSEG \
        -foreach -replace 14 999 -thresh 999 999 1 0 -endfor  \
        -overlap 1 | cut -d , -f 5)

      SUB=($(c3d $GTSEG $AUTOSEG \
        -foreach -replace 6 0 7 0 13 0 14 0 -endfor \
        -label-overlap \
        |  awk '{print $4}' | awk '{printf("%s ", $1)}'  | awk '{print $6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}'))
      CA1="${SUB[0]}"
      CA2="${SUB[1]}"
      DG="${SUB[2]}"
      CA3="${SUB[3]}"
      TAIL="${SUB[4]}"
      SUB="${SUB[5]}"
      ERC="${SUB[6]}"
      BA35="${SUB[7]}"
      BA36="${SUB[8]}"
      PHC="${SUB[9]}"

      str="$id,$side,$NATLAS,0,$GMALLTMP,$HIPPO,$CA1,$CA2,$CA3,$CA,$DG,$SUB,$TAIL,$ERC,$BA35,$BA36,$PHC,$HippoSulcus,$CS,$MISC1,$MISC2"
      echo $str >> $OUTCSV

    done
  done
}

##################################################
function TrainNumatlasVSPerformance()
{
  mkdir -p $EXPTRUEVARYATLASDIR
  Nbegin=1
  #N=$(cat $EXPSETUPDIR/vary_Natlases_xval_fold${expid}.txt | wc -l)
  N=1
  for ((i=$Nbegin;i<=${N};i++)); do

    # prepare manifest file
    XVALDIR=$EXPTRUEVARYATLASDIR/xval/xval$(printf %04d $i)
    XVALTRAINDIR=$XVALDIR/train
    mkdir -p $XVALTRAINDIR
    rm -rf $XVALTRAINDIR/manifest.txt
    for id in $(cat $ANALYSISDIR/subjlist_withflip.txt); do

      id="$id"
      set +e
      exist=$(cat $EXPSETUPDIR/vary_Natlases_xval_fold${expid}.txt | head -n $i | tail -n 1 | grep $id)
      if [[ $exist != "" ]]; then
        continue
      fi
      MPRAGE=$(ls $RAWDATADIR/$id/*mprage.nii.gz)
      TSE=$(ls $RAWDATADIR/$id/*tse.nii.gz)
      LSEG=$(ls $RAWDATADIR/$id/*_left_finalLW.nii.gz)
      RSEG=$(ls $RAWDATADIR/$id/*_right_finalLW.nii.gz)
      echo "$id \
            $MPRAGE \
            $TSE \
            $LSEG \
            $RSEG" \
            >> $XVALTRAINDIR/manifest.txt

    done

    # run training
    time $ASHS_ROOT/bin/ashs_train.sh \
        -D $XVALTRAINDIR/manifest.txt \
        -L $ANALYSISDIR/Labels_3T_Atlas.txt \
        -w $XVALTRAINDIR \
        -C $ANALYSISDIR/ashs_config_vacind_multimodal.sh \
        -N \
        -s 1-7 \
        -Q -z $CODEDIR/ashs-fast-z.sh \
        -B \
        | tee -a $XVALTRAINDIR/ashs_train.out.txt &

    

  done
}


##################################################
# Main entrypoint
cmd=$0
#if [[ $# -lt 2 ]]; then
if [[ $cmd == $RUNCMD1 || $cmd == $RUNCMD2 ]]; then

  main $@

else
  cmd=$1
  shift
  $cmd $@
fi
