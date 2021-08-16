#!/bin/bash 
#$ -S /bin/bash
set -x -e

######################################################################################
# find similarity for each side seperately
######################################################################################

expid=01
expdir=exp${expid}

#ANTSPATH=/home/avants/bin/ants
ANTSPATH=~pauly/bin/ants/
PATH=$ANTSPATH:$PATH
C3D_BIN=/home/hwang3/ahead/turnkey/ext/Linux/bin

LABEL_IDS=(BKG CA DG SUB ERC BA35 BA36 CS)
LABEL_MRG=("0" "1 2 4" "3" "8" "9" "11" "12" "13")
KINDS="tse mprage ${LABEL_IDS[*]}"

# Relevant labels
LABEL_FG=(CA DG SUB ERC BA35 BA36)

ROOT=/home/longxie/ASHS/
ASHSRUNDIR=$ROOT/ashs
FIXUPDIR=$ROOT/ashs_output_truexval

#SUBJ_TXT=subj.txt
SUBJ_TXT=../analysis_input/subj.txt

mkdir -p ${expdir}/dump

####################################################################################
# Transform data to template (ASHS)
####################################################################################
function copy_subject()
{

  fn=$1
  side=$2

  ### SEG=$ASHSRUNDIR/$fn/bootstrap/fusion/lfseg_heur_${side}.nii.gz
  SEG=$FIXUPDIR/${fn}_seg_${side}.nii.gz

  # Link the subfield images
  if [[ -f $SEG ]]; then
    ln -sf $ASHSRUNDIR/$fn/tse_to_chunktemp_${side}.nii.gz ${expdir}/data/${fn}_${side}_tse.nii.gz
    ln -sf $ASHSRUNDIR/$fn/mprage_to_chunktemp_${side}.nii.gz ${expdir}/data/${fn}_${side}_mprage.nii.gz
  fi

  # Generate a binary image for the label
  for ((i=0;i<${#LABEL_IDS[*]};i++)); do

    c3d $SEG -replace $(for k in ${LABEL_MRG[i]}; do echo $k 999; done) -thresh 999 999 1 0 \
      -o $TMPDIR/binary_${LABEL_IDS[i]}_${fn}_${side}.nii.gz

    WarpImageMultiTransform 3  $TMPDIR/binary_${LABEL_IDS[i]}_${fn}_${side}.nii.gz \
      ${expdir}/data/${fn}_${side}_${LABEL_IDS[i]}.nii.gz \
      -R ${expdir}/data/${fn}_${side}_tse.nii.gz \
      $ASHSRUNDIR/$fn/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz \
      $ASHSRUNDIR/$fn/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
      $ASHSRUNDIR/$fn/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

  done

  # Vote in the subject space
  c3d $(for sub in ${LABEL_IDS[*]}; do echo ${expdir}/data/${fn}_${side}_${sub}.nii.gz; done) \
    -vote -type ushort -o ${expdir}/data/${fn}_${side}_seg.nii.gz

}

function copy_data() 
{

  # Get the data
  mkdir -p ${expdir}/data
  for id in $(cat $SUBJ_TXT); do
    fn=$(ls ${ASHSRUNDIR} | grep $id)
    for side in left right; do
    
      qsub -V -cwd -o ${expdir}/dump -j y \
        -N "CP${expid}_${fn}_${side}" $0 copy_subject $fn $side

    done
  done

  qsub -V -cwd -o ${expdir}/dump -j y -hold_jid "CP${expid}_*" -sync y -b y sleep 1

}


########################################################################
# Pairwise registration
########################################################################
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

      # Convert that to ITK format
      c3d_affine_tool \
        $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine.txt \
        -oitk $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine_itk.txt


      #WGT=$(echo ${#LABEL_IDS[*]} | awk '{print 1.0 / $1}')
      WGT=1

      CMD=""
      for sub in ${LABEL_IDS[*]}; do
        CMD="$CMD -m MSQ[${expdir}/data/${fn_fix}_${side}_${sub}.nii.gz,${expdir}/data/${fn_mov}_${side}_${sub}.nii.gz,$WGT]"
      done

      # Perform ANTs registration
      ANTS 3 $CMD \
        -t SyN[0.25] -r Gauss[0.5,0] -i 60x30x10 \
        -a $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine_itk.txt \
        --continue-affine 0 \
        -o $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_pairwise.nii.gz \
        | tee $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_antsoutput.txt

      for sub in ${LABEL_IDS[*]}; do

        WarpImageMultiTransform 3 \
          ${expdir}/data/${fn_mov}_${side}_${sub}.nii.gz \
          $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_${sub}.nii.gz \
          -R ${expdir}/data/${fn_fix}_${side}_tse.nii.gz \
          $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_pairwiseWarp.nii.gz \
          $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_pairwiseAffine.txt

      done

      # Create seg 
      c3d $(for sub in ${LABEL_IDS[*]}; do echo $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_${sub}.nii.gz; done) \
        -vote -type ushort \
        -o $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz

      rm -f $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine.txt
      rm -f $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine_itk.txt
      rm -f $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_pairwiseWarp.nii.gz
      rm -f $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_pairwiseInverseWarp.nii.gz
      rm -f $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_pairwiseAffine.txt

      fi

    fi

  done

}

function pairwise()
{

  PWREG_DIR=${expdir}/pairwise
  #rm -rf ${PWREG_DIR}
  IDS=$(cat $SUBJ_TXT)

  for side in left right; do

    for id_fix in $IDS; do

      qsub -V -v ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2 -pe serial 2 \
           -cwd -o ${expdir}/dump -j y \
           -N "PW${expid}_${id_fix}_${side}" $0 \
           pairwise_qsub $id_fix $side
#-pe serial 2 -l h_stack=128M \

    done

  done

#  qsub -V -cwd -o ${expdir}/dump -j y \
#   -hold_jid "PW${expid}_*" -sync y -b y sleep 1

}

####################################################################################
# prepration
####################################################################################
function make_images_qsub()
{

  id_fix=$1
  side=$2

  IDS=$(ls data | grep ${side}_tse | sed -e "s/_${side}_tse.nii.gz//")
  for id_mov in $IDS; do

    if [[ $id_mov != $id_fix ]]; then
      
      OUTDIR=${expdir}/pairwise/${id_fix}_${side}/${id_mov}_to_${id_fix}

      # Create PRC for each subject
      #c3d \
      #  $OUTDIR/${id_mov}_to_${id_fix}_${side}_reslice_seg.nii.gz \
      #  -replace 1 0 2 0 3 0 4 0 7 0 \
      #  -o $OUTDIR/${id_mov}_to_${id_fix}_${side}_reslice_seg_PRC.nii.gz

      # Create BA35 for each subject
      c3d \
        $OUTDIR/${id_mov}_to_${id_fix}_${side}_reslice_seg.nii.gz \
        -replace 1 0 2 0 3 0 4 0 6 0 7 0 \
        -o $OUTDIR/${id_mov}_to_${id_fix}_${side}_reslice_seg_BA35.nii.gz

      # Create seg for each subject
      #c3d $(for sub in ${LABEL_IDS[*]}; do echo ${expdir}/work/${id}_${side}_totemp_reslice_${sub}.nii.gz; done) \
      #   -vote -type ushort -o ${expdir}/work/${id}_${side}_totemp_reslice_seg.nii.gz

      # Create seg expectation for each subject
      #c3d \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_BKG.nii.gz -scale 0 \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_CA.nii.gz -scale 1 \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_DG.nii.gz -scale 2 \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_SUB.nii.gz -scale 3 \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_ERC.nii.gz -scale 4 \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_BA35.nii.gz -scale 5 \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_BA36.nii.gz -scale 6 \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_CS.nii.gz -scale 7 \
      #  -add -add -add -add -add -add -add \
      #  -o ${expdir}/work/${id}_${side}_totemp_reslice_seg_exp.nii.gz

      # Create PRC exp for each subject
      #c3d \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_BA35.nii.gz -scale 5 \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_BA36.nii.gz -scale 6 \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_CS.nii.gz -scale 7 \
      #  -add -add \
      #  -o ${expdir}/work/${id}_${side}_totemp_reslice_PRC_exp.nii.gz

      # Create BA35 label for each subject
      #c3d \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_seg.nii.gz \
      #  -replace 1 0 2 0 3 0 4 0 6 0 7 0 \
      #  -o ${expdir}/work/${id}_${side}_totemp_reslice_seg_BA35.nii.gz

      # Create PRC label for each subject
      #c3d \
      #  ${expdir}/work/${id}_${side}_totemp_reslice_seg.nii.gz \
      #  -replace 1 0 2 0 3 0 4 0 5 1 6 1 7 0 \
      #  -o ${expdir}/work/${id}_${side}_totemp_reslice_label_PRC.nii.gz

    fi

  done

}

function make_images()
{

  for side in left right; do

    IDS=$(ls data | grep ${side}_tse | sed -e "s/_${side}_tse.nii.gz//")

    for id_fix in $IDS; do

      #qsub -V -cwd -o ${expdir}/dump -j y \
      #     -N "MI${expid}_${id_fix}_${side}" $0 \
      #     make_images_qsub $id_fix $side
     
      make_images_qsub $id_fix $side

    done
  done

  #qsub -V -cwd -o ${expdir}/dump -j y \
  #  -hold_jid "MI${expid}_*" -sync y -b y sleep 1

}

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

########################################################################
# compute similarity
########################################################################
function sim()
{

  fn_seg_fix=$1
  fn_seg_mov=$2
  fn_out=$3
  sim_type=$4

  if [[ $sim_type == 'dice' ]]; then
    
    c3d ${fn_seg_fix} ${fn_seg_mov} -label-overlap \
      | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}' \
      > ${fn_out}

  elif [[ $sim_type == 'nmi' ]]; then

    c3d ${fn_seg_fix} ${fn_seg_mov} -nmi \
      | awk '{print $3}' > ${fn_out}

  elif [[ $sim_type == 'msq' ]]; then

    c3d ${fn_seg_fix} ${fn_seg_mov} -msq \
      | awk '{print $3}' > ${fn_out}

  fi

}

########################################################################
# Method 1
########################################################################
function sim_seg_dice_qsub()
{

  id_fix=$1
  side=$2

  IDS=$(cat $SUBJ_TXT)
  fn_fix=$(ls ashs | grep $id_fix)

  for id_mov in $IDS; do

    fn_out=${expdir}/sim_seg_dice/${side}/${id_mov}_to_${id_fix}.txt

    if [[ $id_mov == $id_fix ]]; then

      echo "1" > ${fn_out}

    else

      fn_mov=$(ls ashs | grep $id_mov)
      fn_seg_fix=data/${fn_fix}_${side}_seg.nii.gz
      fn_seg_mov=${expdir}/pairwise/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz
      sim_type=dice

      sim $fn_seg_fix $fn_seg_mov $fn_out $sim_type

    fi

  done

}

function sim_seg_dice()
{
  
  SIM_DIR=${expdir}/sim_seg_dice
  rm -rf ${SIM_DIR}
  IDS=$(cat $SUBJ_TXT)

  for side in left right; do

    mkdir -p ${SIM_DIR}/${side}
    echo '' > ${SIM_DIR}/sim_seg_dice_${side}.txt  

    for id_fix in $IDS; do

      qsub -V -cwd -o ${expdir}/dump -j y \
           -N "sim${expid}_${id_fix}" $0 \
           sim_seg_dice_qsub $id_fix $side    

    done

    qsub -V -cwd -o ${expdir}/dump -j y \
         -hold_jid "sim${expid}_*" -sync y -b y sleep 1


    # Concadenate the results
    #echo '' > ${SIM_DIR}/sim_PRC_seg_dice_${side}.txt
    rm -rf ${SIM_DIR}/sim_seg_dice_${side}.txt

    for id_fix in $IDS; do

      for id_mov in $IDS; do

        cat ${SIM_DIR}/${side}/${id_mov}_to_${id_fix}.txt \
            >> ${SIM_DIR}/sim_seg_dice_${side}.txt

      done

    done
    
    rm -rf ${SIM_DIR}/${side}

  done

}

########################################################################
# Method 2
########################################################################
function sim_PRC_seg_dice_qsub()
{

  id_fix=$1
  side=$2

  IDS=$(cat $SUBJ_TXT)
  fn_fix=$(ls ashs | grep $id_fix)

  for id_mov in $IDS; do

    fn_out=${expdir}/sim_PRC_seg_dice/${side}/${id_mov}_to_${id_fix}.txt

    if [[ $id_mov == $id_fix ]]; then

      echo "1" > ${fn_out}

    else

      fn_mov=$(ls ashs | grep $id_mov)
      c3d data/${fn_fix}_${side}_seg.nii.gz \
        -replace 1 0 2 0 3 0 4 0 7 0 \
        -o data/${fn_fix}_${side}_seg_PRC.nii.gz

      fn_seg_fix=data/${fn_fix}_${side}_seg_PRC.nii.gz
      fn_seg_mov=${expdir}/pairwise/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg_PRC.nii.gz
      sim_type=dice

      sim $fn_seg_fix $fn_seg_mov $fn_out $sim_type

    fi
  done

}

function sim_PRC_seg_dice()
{

  SIM_DIR=${expdir}/sim_PRC_seg_dice
  rm -rf ${SIM_DIR}
  IDS=$(cat $SUBJ_TXT)

  for side in left right; do

    mkdir -p ${SIM_DIR}/${side}
    #echo '' > ${SIM_DIR}/sim_PRC_seg_dice_${side}.txt

    for id_fix in $IDS; do

      qsub -V -cwd -o ${expdir}/dump -j y \
           -N "sim${expid}_${id_fix}" $0 \
           sim_PRC_seg_dice_qsub $id_fix $side

    done

    qsub -V -cwd -o ${expdir}/dump -j y \
         -hold_jid "sim${expid}_*" -sync y -b y sleep 1

    # Concadenate the results
    #echo '' > ${SIM_DIR}/sim_PRC_seg_msq_${side}.txt
    rm -rf ${SIM_DIR}/sim_PRC_seg_dice_${side}.txt

    for id_fix in $IDS; do

      for id_mov in $IDS; do

        cat ${SIM_DIR}/${side}/${id_mov}_to_${id_fix}.txt \
            >> ${SIM_DIR}/sim_PRC_seg_dice_${side}.txt

      done

    done

    rm -rf ${SIM_DIR}/${side}

  done

}

########################################################################
# Method 3
########################################################################
function sim_BA35_seg_dice_qsub()
{

  id_fix=$1
  side=$2

  IDS=$(cat $SUBJ_TXT)
  fn_fix=$(ls ashs | grep $id_fix)

  for id_mov in $IDS; do

    fn_out=${expdir}/sim_BA35_seg_dice/${side}/${id_mov}_to_${id_fix}.txt

    if [[ $id_mov == $id_fix ]]; then

      echo "1" > ${fn_out}

    else

      fn_mov=$(ls ashs | grep $id_mov)
      c3d data/${fn_fix}_${side}_seg.nii.gz \
        -replace 1 0 2 0 3 0 4 0 6 0 7 0 \
        -o data/${fn_fix}_${side}_seg_BA35.nii.gz

      fn_seg_fix=data/${fn_fix}_${side}_seg_BA35.nii.gz
      fn_seg_mov=${expdir}/pairwise/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg_BA35.nii.gz
      sim_type=dice

      sim $fn_seg_fix $fn_seg_mov $fn_out $sim_type

    fi
  done

}


function sim_BA35_seg_dice()
{

  SIM_DIR=${expdir}/sim_BA35_seg_dice
  rm -rf ${SIM_DIR}
  IDS=$(cat $SUBJ_TXT)

  for side in left right; do

    mkdir -p ${SIM_DIR}/${side}

    for id_fix in $IDS; do

      qsub -V -cwd -o ${expdir}/dump -j y \
           -N "sim${expid}_${id_fix}" $0 \
           sim_BA35_seg_dice_qsub $id_fix $side

    done

    qsub -V -cwd -o ${expdir}/dump -j y \
         -hold_jid "sim${expid}_*" -sync y -b y sleep 1

    # Concadenate the results
    #echo '' > ${SIM_DIR}/sim_PRC_seg_msq_${side}.txt
    rm -rf ${SIM_DIR}/sim_BA35_seg_dice_${side}.txt

    for id_fix in $IDS; do

      for id_mov in $IDS; do

        cat ${SIM_DIR}/${side}/${id_mov}_to_${id_fix}.txt \
            >> ${SIM_DIR}/sim_BA35_seg_dice_${side}.txt

      done

    done

    rm -rf ${SIM_DIR}/${side}

  done

}


########################################################################
# Method 4
########################################################################
function sim_PRC_label_dice_qsub()
{

  id_fix=$1
  side=$2

  IDS=$(cat $SUBJ_TXT)
  fn_fix=$(ls ashs | grep $id_fix)

  for id_mov in $IDS; do
    fn_mov=$(ls ashs | grep $id_mov)

    fn_seg_fix=${expdir}/work/${fn_fix}_${side}_totemp_reslice_label_PRC.nii.gz
    fn_seg_mov=${expdir}/work/${fn_mov}_${side}_totemp_reslice_label_PRC.nii.gz
    fn_out=${expdir}/sim_PRC_label_dice/${side}/${id_mov}_to_${id_fix}.txt
    sim_type=dice

    sim $fn_seg_fix $fn_seg_mov $fn_out $sim_type

  done

}

function sim_PRC_label_dice()
{

  SIM_DIR=${expdir}/sim_PRC_label_dice
  rm -rf ${SIM_DIR}
  IDS=$(cat $SUBJ_TXT)

  for side in left right; do

    mkdir -p ${SIM_DIR}/${side}
    echo '' > ${SIM_DIR}/sim_PRC_label_dice_${side}.txt

    for id_fix in $IDS; do

      qsub -V -cwd -o ${expdir}/dump -j y \
           -N "sim${expid}_${id_fix}" $0 \
           sim_PRC_label_dice_qsub $id_fix $side

    done

    qsub -V -cwd -o ${expdir}/dump -j y \
         -hold_jid "sim${expid}_*" -sync y -b y sleep 1


    # Concadenate the results
    #echo '' > ${SIM_DIR}/sim_PRC_label_dice_${side}.txt
    rm -rf ${SIM_DIR}/sim_PRC_label_dice_${side}.txt

    for id_fix in $IDS; do

      for id_mov in $IDS; do

        cat ${SIM_DIR}/${side}/${id_mov}_to_${id_fix}.txt \
            >> ${SIM_DIR}/sim_PRC_label_dice_${side}.txt

      done

    done

  done

}





function reset_dir()
{
  #rm -rf data
  #rm -rf ${expdir}/work ${expdir}/meshwarp ${expdir}/jacobian ${expdir}/segerrormap
  rm -rf ${expdir}/dump/*
}



if [[ $1 == "copy_subject" ]]; then

  copy_subject $2 $3


elif [[ $1 == "pairwise_qsub" ]]; then

  pairwise_qsub $2 $3

elif [[ $1 == "make_images_qsub" ]]; then

  make_images_qsub $2 $3

elif [[ $1 == "sim_seg_dice_qsub" ]]; then

  sim_seg_dice_qsub $2 $3

elif [[ $1 == "sim_PRC_seg_dice_qsub" ]]; then

  sim_PRC_seg_dice_qsub $2 $3

elif [[ $1 == "sim_BA35_seg_dice_qsub" ]]; then

  sim_BA35_seg_dice_qsub $2 $3

elif [[ $1 == "sim_PRC_label_dice_qsub" ]]; then

  sim_PRC_label_dice_qsub $2 $3

else

  #reset_dir
  #copy_data
  #pairwise
  #make_images
  completeness
  #sim_seg_dice
  #sim_PRC_seg_dice
  #sim_BA35_seg_dice
  
  #sim_PRC_label_dice

  #IDS=$(cat $SUBJ_TXT)
  #for id in $IDS; do
  #  fn_id=$(ls ashs | grep $id)
  #  echo $fn_id>>subj_full.txt
  #done

  exit

fi

