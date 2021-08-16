#!/bin/bash
#$ -S /bin/bash
set -x -e

##############################################################################
# INIT
#ANTSPATH=/home/avants/bin/ants
ANTSPATH=~pauly/bin/ants/
NEWANTSDIR=/share/apps/ANTs/2014-06-23/build/bin/
PATH=$ANTSPATH:$PATH

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

LABEL_IDS=(BKG             ERC  BA35 BA36 PHC  CS   OTS)
LABEL_MRG=("0 1 2 3 4 7 8" "10" "11" "12" "13" "14" "16")
LABEL_NEW=(0               1    2    3    4    5    6   )
KINDS="tse mprage ${LABEL_IDS[*]}"

# Directories
ROOT=/home/longxie/ASHS_PHC
WORKDIR=$ROOT/thickness_newlabel
ASHSRUNDIR=$ROOT/SuperResolution/ASHSexp/exp002/fullset/ashs

SUBJ_TXT=$WORKDIR/analysis_input/subj.txt

##############################################################################
# Parameters needs to be specify

# 1. Experiment number
expid=002
expdir=$WORKDIR/SCexp/exp${expid}
DUMPDIR=${expdir}/dump
mkdir -p ${expdir}/dump

# 2. ANTS parameters
WGT=1
REG_LABELS=(${LABEL_IDS[*]})
ANTs_t="SyN[0.25]"
ANTs_r="Gauss[0.5,0]"
ANTs_i="30x20x10"
ANTs_all_metrics="--use-all-metrics-for-convergence"

# 3. similarity type
SIM_TYPE="PRCCS_seg_dice"

#HOW ABOUT ONLY USE CS AT PRC!!!!

#######################################################################
function main()
{
  reset_dir
  #copy_data
  #pairwise
  #completeness
  similarity

}

#######################################################################
# Copy data
function copy_data()
{
  # Get the data
  mkdir -p ${expdir}/data

  # Load id number
  IDs=($(cat $SUBJ_TXT))

  # Submit job to copy data
  PREFIX=CP${expid}
  for side in left right; do
    for ((i=0;i<${#IDs[*]};i++)); do

      id=${IDs[i]}
      fn=$(ls $ASHSRUNDIR | grep $id)

      qsub -cwd -o $DUMPDIR -j y -N \
           "${PREFIX}_${id}_${side}" \
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
  SEG=$ASHSRUNDIR/${fn}/final/${fn}_${side}_lfseg_corr_usegray.nii.gz

  # Link the subfield images
  if [[ -f $SEG ]]; then
    ln -sf $ASHSRUNDIR/$fn/tse_to_chunktemp_${side}.nii.gz \
           ${expdir}/data/${fn}_${side}_tse.nii.gz
    ln -sf $ASHSRUNDIR/$fn/mprage_to_chunktemp_${side}.nii.gz \
           ${expdir}/data/${fn}_${side}_mprage.nii.gz
  fi

  # Generate a binary image for the label
  for ((i=0;i<${#LABEL_IDS[*]};i++)); do

    c3d $SEG -replace $(for k in ${LABEL_MRG[i]}; do echo $k 999; done) -thresh 999 999 1 0 \
      -o $TMPDIR/binary_${LABEL_IDS[i]}_${fn}_${side}.nii.gz

    WarpImageMultiTransform 3  \
      $TMPDIR/binary_${LABEL_IDS[i]}_${fn}_${side}.nii.gz \
      ${expdir}/data/${fn}_${side}_${LABEL_IDS[i]}.nii.gz \
      -R ${expdir}/data/${fn}_${side}_tse.nii.gz \
      $ASHSRUNDIR/$fn/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz \
      $ASHSRUNDIR/$fn/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
      $ASHSRUNDIR/$fn/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

  done

  # Vote in the ASHS template subject space to create segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo ${expdir}/data/${fn}_${side}_${sub}.nii.gz; done) \
    -vote -type ushort -o ${expdir}/data/${fn}_${side}_seg.nii.gz
}

####################################################################
# Pairwise registration
function pairwise()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=PW${expid}
  for side in left right; do

    for id_fix in $IDS; do

      qsub -cwd -o $DUMPDIR -j y \
           -N "${PREFIX}_${id_fix}_${side}" $0 \
           pairwise_qsub $id_fix $side
      sleep 0.1

    done

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function pairwise_qsub()
{
  id_fix=$1
  side=$2

  IDS=$(cat $SUBJ_TXT)
  fn_fix=$(ls $ASHSRUNDIR | grep $id_fix)
  for id_mov in $IDS; do

    if [[ $id_mov != $id_fix ]]; then

      fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)
      OUTDIR=${expdir}/pairwise/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}
      mkdir -p $OUTDIR

      if [[ -f $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz ]]; then

        echo "Seg file exists."

      else

      # Use ml_affine for nice affine alignment
      ~pauly/wolk/ashs/ext/Linux/bin/ml_affine \
        ${expdir}/data/${fn_fix}_${side}_seg.nii.gz \
        ${expdir}/data/${fn_mov}_${side}_seg.nii.gz \
        $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine.txt

      for sub in ${KINDS[*]}; do

        c3d \
        ${expdir}/data/${fn_fix}_${side}_tse.nii.gz \
        ${expdir}/data/${fn_mov}_${side}_${sub}.nii.gz \
        -reslice-matrix \
        $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine.txt \
        -o $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_${sub}.nii.gz

      done

      # Create seg
      c3d $(for sub in ${LABEL_IDS[*]}; do echo $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_${sub}.nii.gz; done) \
        -vote -type ushort \
        -o $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz

      rm -f $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine.txt

      fi
    fi

  done
}

##############################################################################
function completeness
{
  for side in left right; do

    IDS=$(ls ${expdir}/data | grep ${side}_tse | sed -e "s/_${side}_tse.nii.gz//")

    echo "" > ${expdir}/check_pairwise.txt

    for id_fix in $IDS; do

      for id_mov in $IDS; do

      if [[ $id_fix != $id_mov ]]; then

      # Check whether all the files exist
      OUTDIR=${expdir}/pairwise/${id_fix}_${side}/${id_mov}_to_${id_fix}
      if [[ -f $OUTDIR/${id_mov}_to_${id_fix}_${side}_reslice_seg.nii.gz ]]; then
          echo "ok"
      else
        echo "${id_mov} to ${id_fix} of ${side} is missing" \
          >> ${expdir}/check_pairwise.txt
      fi

      fi
      done
    done
  done
}

##############################################################################
function similarity()
{
  SIM_DIR=${expdir}/sim_${SIM_TYPE}
  rm -rf ${SIM_DIR}
  IDS=$(cat $SUBJ_TXT)

  PREFIX=SIM${expid}
  for side in left right; do

    mkdir -p ${SIM_DIR}/${side}

    for id_fix in $IDS; do

      qsub -cwd -o $DUMPDIR -j y \
           -N "${PREFIX}_${id_fix}" $0 \
           similarity_sub $id_fix $side
      sleep 0.1

    done

    qsub -cwd -o $DUMPDIR -j y \
         -hold_jid "${PREFIX}_*" -sync y -b y \
         sleep 1

    # Concadenate the results
    #echo '' > ${SIM_DIR}/sim_PRC_seg_msq_${side}.txt
    rm -rf ${SIM_DIR}/sim_${SIM_TYPE}_${side}.txt

    for id_fix in $IDS; do

        cat ${SIM_DIR}/${side}/${id_fix}.txt \
            >> ${SIM_DIR}/sim_${SIM_TYPE}_${side}.txt

    done

    rm -rf ${SIM_DIR}/${side}

  done
}

function similarity_sub()
{
  id_fix=$1
  side=$2
  fn_fix=$(ls $ASHSRUNDIR | grep $id_fix)
  fn_out=${expdir}/sim_${SIM_TYPE}/${side}/${id_fix}.txt
  fn_seg_fix=$TMPDIR/${fn_fix}_${side}_seg_tmp.nii.gz
  OVL=""

  if [ ${SIM_TYPE} == "PRC_seg_dice" ]; then

    c3d ${expdir}/data/${fn_fix}_${side}_seg.nii.gz \
      -thresh 2 3 1 0 \
      -o $fn_seg_fix

  elif [ ${SIM_TYPE} == "CS_seg_dice" ]; then

    c3d ${expdir}/data/${fn_fix}_${side}_seg.nii.gz \
      -thresh 5 5 1 0 \
      -o $fn_seg_fix

  elif [ ${SIM_TYPE} == "PRCCS_seg_dice" ]; then

    c3d ${expdir}/data/${fn_fix}_${side}_seg.nii.gz \
      -replace 1 0 4 0 6 0 \
      -o $fn_seg_fix

  fi


  # Go through all the other subjects
  IDS=$(cat $SUBJ_TXT)
  for id_mov in $IDS; do

    fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)
    OVL_tmp="NA"

    if [[ $id_mov == $id_fix ]]; then

      OVL="$OVL 1"

    else

      if [ ${SIM_TYPE} == "PRC_seg_dice" ]; then

        OVL_tmp=$(c3d ${expdir}/pairwise/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz \
          -thresh 2 3 1 0 \
          $fn_seg_fix -overlap 1 \
          | grep OVL | awk -F '[ ,]+' '{print $6}' )

      elif [ ${SIM_TYPE} == "CS_seg_dice" ]; then

        OVL_tmp=$(c3d ${expdir}/pairwise/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz \
          -thresh 5 5 1 0 \
          $fn_seg_fix -overlap 1 \
          | grep OVL | awk -F '[ ,]+' '{print $6}' )

      elif [ ${SIM_TYPE} == "PRCCS_seg_dice" ]; then

        OVL_tmp=$(c3d ${expdir}/pairwise/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz \
          -replace 1 0 4 0 6 0 \
          $fn_seg_fix -label-overlap \
          | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}' )

      fi

      OVL="$OVL $OVL_tmp"

    fi
  done

  echo $OVL > ${fn_out}

}

##############################################################################
function reset_dir()
{
  rm -rf ${expdir}/dump/*
}

##############################################################################
if [[ $1 == "copy_subject" ]]; then

  copy_subject $2 $3

elif [[ $1 == "pairwise_qsub" ]]; then

  pairwise_qsub $2 $3

elif [[ $1 == "similarity_sub" ]]; then

  similarity_sub $2 $3

else

  main

  exit

fi


