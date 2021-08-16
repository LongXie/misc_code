#!/bin/bash
#$ -S /bin/bash
set -e  -x

ANTSPATH=/data/picsl/pauly/bin/ants/
NEWANTSDIR=/share/apps/ANTs/2014-06-23/build/bin/
export ASHS_ROOT=/data/picsl/longxie/pkg/ashs/ashs-fast-beta
export PATH=$PATH:$ASHS_ROOT/bin
MATLAB_BIN=/share/apps/matlab/R2013a/bin/matlab
export PATH=$ANTSPATH:$PATH

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

# directories
expid="7T"
ROOT=/home/longxie/DeepLabelFusion/PMC_Prisma_3T_withflip
ANALYSISDIR=$ROOT/analysis_input
EXPSETUPDIR=$ROOT/expsetup
EXPROOT=$ROOT/exp_ASHST1/
CODEDIR=$ROOT/scripts
TEMPLATEDIR=$EXPROOT/template_all_7T
VOLUMEDIR=$EXPROOT/refseg_volume
RAWDATADIR=$ROOT/data
#RAWDATADIR=/data/picsl/lwisse/NACC_Prisma_Atlas_NotFinal
AUTOSEGDIR=$TEMPLATEDIR/autoseg_7T
#AUTOSEGDIR=$TEMPLATEDIR/autoseg_fold0

ASHSRUNDIR=$TEMPLATEDIR/train/atlas
#MATLABCODEDIR=$WORKDIR/matlabcode/
#SUBJ_TXT=$ANALYSISDIR/subjlist_withflip.txt
#SUBJ_TXT=/home/longxie/DeepLabelFusion/PMC_Prisma_3T_withflip/expsetup/test001_to_003.txt
SUBJ_TXT_7T=/home/longxie/ASHS_7TT2/analysis_input/7TAtlasIDs.csv
SUBJ_TXT=$TEMPLATEDIR/IDs.txt
ASHSTEMPDIR=$TEMPLATEDIR/train/template_build

expdir=$TEMPLATEDIR/MTL_template
DUMPDIR=${expdir}/dump
mkdir -p ${expdir}/dump
FLIPXFN=$ANALYSISDIR/flipx_itk.txt
FLIPYFN=$ANALYSISDIR/flipy_itk.txt
FLIPXFN=$ANALYSISDIR/flipz_itk.txt
IDENTITYFN=$ANALYSISDIR/identity_itk.txt
#INITTEMP_DIR=$expdir/singletemplate
EVALDIR=$expdir/evaluation

########################################
# 0. copy data
LABEL_IDS_ALL=(BKG              CA1 CA2 DG CA3 Tail SUB ERC BA35 BA36 PHC CS)
LABEL_MRG_ALL=("0 6 7 13 14 15" 1   2   3  4   5    8   9   10   11   12  16)
# Labels to get segmentation
LABEL_IDS=(BKG              CA1 CA2 DG CA3 Tail SUB ERC BA35 BA36 PHC CS)
LABEL_MRG=("0 6 7 13 14 15" 1   2   3  4   5    8   9   10   11   12  16)
LABEL_NEW=(0                1   2   3  4   5    6   7   8    9    10  11)
KINDS="tse mprage ${LABEL_IDS[*]}"
DATADIR=${expdir}/data

# 2. Group settings (can be single group)
# array indicating whether certain group exits (12 max)
groups=(1)

# 3. ANTS parameters
ITER=8
ANTs_start=3
WGT=1
REG_LABELS=(${LABEL_IDS[*]})
ANTs_t="SyN[0.25]"
ANTs_r="Gauss[0.5,0]"
ANTs_i="80x80x20"
ANTs_x="Y"
ANTs_all_metrics="--use-all-metrics-for-convergence"

# 4. Mesh extraction
subfield_th=0.5
template_th=0.0

MESH_EVAL=(       CA1 CA2 DG CA3 Tail SUB ERC BA35 BA36 PHC CS)
MESH_EVAL_RANGES=(1   2   3  4   5    6   7   8    9    10  11)

# Relevant labels
LABEL_FG=(CA1 CA2 DG CA3 Tail SUB ERC BA35 BA36 PHC CS)
EXTRAMESHES=(Hippo HippoNoDG Cortex MRG)
#                BKG CA1 CA2 DG CA3 Tail SUB ERC BA35 BA36 PHC CS
EXTRAMESHESDEF=("-1  1   1   1  1   1    1   -1  -1   -1   -1  -1" \
                "-1  1   1   -1 1   1    1   -1  -1   -1   -1  -1" \
                "-1  -1  -1  -1 -1  -1   -1  1   1    1    1   -1" \
                "-1  1   1  -1  1   1    1   1   1    1    1   -1")
ALLSF="${LABEL_FG[*]} ${EXTRAMESHES[*]}"


# 5. Thickness measurement
thick_p=1.2
thick_e=6



EVALLABELS=(Hippo         CA1 CA2 CA3 CA      DG SUB Tail ERC BA35 BA36 PHC CS)
RANGES=(    "1 2 3 4 5 6" 1   2   4   "1 2 4" 3  6   5    7   8    9    10  11)

TRAINFNS=($(ls $EXPSETUPDIR | grep train | grep .txt))

VOLLABELS=(Hippo          CA1 CA2 CA3 CA      DG SUB Tail ERC BA35 BA36 PHC HippoSulcus CS MISC1 MISC2)
VOLRANGES=("1 2 3 4 5 8"  1   2   4   "1 2 4" 3  8   5    9   10   11   12  14          16 14     15)

##################################################
function main()
{
  reset_dir

  # prepare to run ASHS to build template
  #Prepare
  #TrainASHS

  ###############################################################
  # apply the old one template pipeline to build a template for MTL labels
  #######################################
  #copy_data
  
  #initial_average

  #main_loop

  #make_images

  #warp_labels

  #warp_meshes

  #evaluation

  #gen_temp_label

  ################################################################
  # warp the error maps to the template space and sample to the mesh
  if [[ 1 == 0 ]]; then
  #for method in MV SVWV ashs_heur ashs_usegray UNET DLF; do
  for method in ashs_heur ashs_usegray UNET DLF; do
  #for method in ASHSUsegrayBootPaul ASHSHeurBootPaul; do
  #for method in ASHSUsegrayBootPaulExp02; do
  #for method in DLF ashs_usegray UNET ashs_heur; do
  #for method in UNET ashs_heur; do

    # warp error map of MV
    warp_errormap $method
 
    # sample error map on mesh
    sample_errormap $method

  done
  fi

  # compute difference map between two methods
  if [[ 1 == 0 ]]; then
  for method in ashs_heur ashs_usegray UNET; do

    diff_errormap $method DLF

  done
  fi

  #diff_errormap ASHSUsegrayBootPaulExp02 ASHSUsegrayBootPaul


  #################################################################
  # get volume of all labels from the ground truth
  get_volume
}

##################################################
function Prepare()
{
  mkdir -p $TEMPLATEDIR
  rm -rf $TEMPLATEDIR/manuifest.txt $TEMPLATEDIR/IDs.txt
  N=$(cat $SUBJ_TXT_7T | wc -l)
  SCANDATADIR=/data/jux/wolk_group/Prisma3T/relong
  RAW7TDATADIR=/home/longxie/ASHS_7TT2/rawdata/ABC_7T_Atlas_final

  # ierate
  for ((j=2;j<=${N};j++)); do

    aid=$(cat -A $SUBJ_TXT_7T | head -n $j | tail -n 1 | cut -d, -f 1)
    id=$(echo $aid | cut -d _ -f 1)
    scandate=$(echo $aid | cut -d _ -f 2)
    if [[ $id == "122812" ]]; then
      continue
    fi

    echo "$id \
          $(ls $SCANDATADIR/$id/MRI7T/$scandate/processed/*mp2rinv2.nii.gz) \
          $(ls $SCANDATADIR/$id/MRI7T/$scandate/processed/T2_*.nii.gz) \
          $(ls $RAW7TDATADIR/$aid/${id}_left_LW5.nii.gz) \
          $(ls $RAW7TDATADIR/$aid/${id}_right_LW5.nii.gz)" \
          >> $TEMPLATEDIR/manuifest.txt
    echo "$id" >> $TEMPLATEDIR/IDs.txt
  done
}

##################################################
function TrainASHS()
{
  # run training
  mkdir -p $TEMPLATEDIR/train
  if [[ 1 == 1 ]]; then
    time $ASHS_ROOT/bin/ashs_train.sh \
      -D $TEMPLATEDIR/manuifest.txt \
      -L /home/longxie/ASHS_7TT2/analysis_input/7TAtlas_snaplabel.txt \
      -w $TEMPLATEDIR/train \
      -C /home/longxie/ASHS_7TT2/analysis_input/ashs_user_config_7T_3TT1.sh \
      -N \
      -s 1-3 \
      -Q -z  /home/longxie/ASHS_7TT2/analysis_input/ashs-fast-train-z-7T.sh
  fi
}

#######################################################################
# Copy data
function copy_data()
{
  # perform initial affine registration between left and flip right temps
  mkdir -p $DATADIR/tempreg
  ln -sf $ASHSTEMPDIR/refspace_meanseg_left.nii.gz \
     $DATADIR/tempreg/refspace_meanseg_left.nii.gz
  ln -sf $ASHSTEMPDIR/refspace_meanseg_right.nii.gz \
     $DATADIR/tempreg/refspace_meanseg_right.nii.gz
  #c3d $ASHSTEMPDIR/refspace_meanseg_right.nii.gz -flip z \
  #  -o $DATADIR/tempreg/refspace_meanseg_right_flip.nii.gz

  # register fliped right to left
  #ml_affine $DATADIR/tempreg/refspace_meanseg_left.nii.gz \
  #  $DATADIR/tempreg/refspace_meanseg_right_flip.nii.gz \
  #  $DATADIR/tempreg/right_to_left_mlaffine.txt
  #c3d_affine_tool $DATADIR/tempreg/right_to_left_mlaffine.txt \
  #  -oitk $DATADIR/tempreg/right_to_left_mlaffine_itk.txt

  # combine affine transformations
  #$ANTSPATH/ComposeMultiTransform 3 \
  #  $DATADIR/tempreg/right_to_left_init_affine.txt \
  #  -R $DATADIR/tempreg/right_to_left_mlaffine_itk.txt \
  #  $DATADIR/tempreg/right_to_left_mlaffine_itk.txt \
  #  $FLIPYFN

  greedy -d 3 \
    -i $DATADIR/tempreg/refspace_meanseg_left.nii.gz \
       $DATADIR/tempreg/refspace_meanseg_right.nii.gz \
    -moments 2 \
    -o $DATADIR/tempreg/right_to_left_init_affine.txt
  c3d_affine_tool $DATADIR/tempreg/right_to_left_init_affine.txt \
    -oitk $DATADIR/tempreg/right_to_left_init_affine.txt

  # Load id number
  IDs=($(cat $SUBJ_TXT))

  # Submit job to copy data
  PREFIX=CP${expid}
  for side in left right; do
    for ((i=0;i<${#IDs[*]};i++)); do

      id=${IDs[i]}

      qsub -cwd -o $DUMPDIR -j y \
           -l h_vmem=2.1G,s_vmem=2G \
           -N "${PREFIX}_${id}_${side}" \
           $0 copy_subject $id $side
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
  SEG=$ASHSRUNDIR/${fn}/seg_${side}.nii.gz

  if [[ ! -f $SEG ]]; then
    echo "$SEG does not exist."
    exit
  fi

  if [[ $side == "left" ]]; then
    ADD_ON_TRANS=$IDENTITYFN
    ln -sf $ASHSRUNDIR/$fn/tse_to_chunktemp_${side}.nii.gz \
           $DATADIR/${fn}_${side}_tse.nii.gz
    ln -sf $ASHSRUNDIR/$fn/mprage_to_chunktemp_${side}.nii.gz \
           $DATADIR/${fn}_${side}_mprage.nii.gz
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
    WarpImageMultiTransform 3  \
      $ASHSRUNDIR/$fn/tse_to_chunktemp_${side}.nii.gz \
      $DATADIR/${fn}_${side}_tse.nii.gz \
      -R $DATADIR/tempreg/refspace_meanseg_left.nii.gz \
      $ADD_ON_TRANS
    WarpImageMultiTransform 3  \
      $ASHSRUNDIR/$fn/mprage_to_chunktemp_${side}.nii.gz \
      $DATADIR/${fn}_${side}_mprage.nii.gz \
      -R $DATADIR/tempreg/refspace_meanseg_left.nii.gz \
      $ADD_ON_TRANS
  fi

  # Generate a binary image for the label
  for ((i=0;i<${#LABEL_IDS_ALL[*]};i++)); do

    c3d $SEG -replace $(for k in ${LABEL_MRG_ALL[i]}; do echo $k 999; done) -thresh 999 999 1 0 \
      -o $TMPDIR/binary_${LABEL_IDS_ALL[i]}_${fn}_${side}.nii.gz

    if [[ -f $ASHSRUNDIR/$fn/flirt_t2_to_t1/greedy_t2_to_t1_chunk_${side}_inv.mat ]]; then

    greedy -d 3 \
      -rf $DATADIR/${fn}_${side}_tse.nii.gz \
      -rm $TMPDIR/binary_${LABEL_IDS_ALL[i]}_${fn}_${side}.nii.gz \
          $DATADIR/${fn}_${side}_${LABEL_IDS_ALL[i]}.nii.gz \
      -r $ADD_ON_TRANS \
         $ASHSRUNDIR/$fn/ants_t1_to_temp/greedy_t1_to_template_warp.nii.gz \
         $ASHSRUNDIR/$fn/affine_t1_to_template/t1_to_template_affine.mat \
         $ASHSRUNDIR/$fn/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat \
         $ASHSRUNDIR/$fn/flirt_t2_to_t1/greedy_t2_to_t1_chunk_${side}_inv.mat \

    else

    greedy -d 3 \
      -rf $DATADIR/${fn}_${side}_tse.nii.gz \
      -rm $TMPDIR/binary_${LABEL_IDS_ALL[i]}_${fn}_${side}.nii.gz \
          $DATADIR/${fn}_${side}_${LABEL_IDS_ALL[i]}.nii.gz \
      -r $ADD_ON_TRANS \
         $ASHSRUNDIR/$fn/ants_t1_to_temp/greedy_t1_to_template_warp.nii.gz \
         $ASHSRUNDIR/$fn/affine_t1_to_template/t1_to_template_affine.mat \
         $ASHSRUNDIR/$fn/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat

    fi
         
  done

  # Vote in the ASHS template subject space to create segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $DATADIR/${fn}_${side}_${sub}.nii.gz; done) \
    -vote -type ushort -o $DATADIR/${fn}_${side}_seg.nii.gz

}

#########################################################################
function initial_average()
{
  # Compute initial average for all the groups
  mkdir -p $expdir/work

  # Compute initial average for each subfield mask
  PREFIX=AVG${expid}
  # Average all the input images
  for kind in $KINDS; do
    qsub -cwd -o $DUMPDIR -j y -N \
         "${PREFIX}_${kind}" $0 \
         average_subfield $kind
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # Compute the initial mask for the segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo ${expdir}/work/template_fullchunk_${sub}.nii.gz; done | grep -v BKG) \
    -mean -thresh 1e-5 inf 1 0 -trim 5vox \
    -o ${expdir}/work/template_mask.nii.gz

  # Trim every template component using the mask
  for kind in $KINDS; do
    c3d ${expdir}/work/template_mask.nii.gz \
      ${expdir}/work/template_fullchunk_${kind}.nii.gz -reslice-identity \
      -o ${expdir}/work/template_${kind}.nii.gz

  done
}

# Average subfields
function average_subfield()
{
  kind=$1

  AverageImages 3 \
       ${expdir}/work/template_fullchunk_${kind}.nii.gz \
       0 ${expdir}/data/*_*_${kind}.nii.gz
}

##############################################################################
# Build template using ANTs
function shape_update_to_template()
{
  nouse=$1

  # Borrowed from ANTS buildtemplateparallel

  # Average the warp fields
  local TEMPWARP=${expdir}/work/template_warp.nii.gz
  $NEWANTSDIR/AverageImages 3 $TEMPWARP 0 \
    ${expdir}/work/*_totempWarp.nii.gz

  # Scale by -0.25 (gradient step)
  MultiplyImages 3 $TEMPWARP -0.25 $TEMPWARP

  # Create a new template affine
  local TEMPAFF=${expdir}/work/template_Affine.txt
  if [[ -f $TEMPAFF ]]; then rm -rf $TEMPAFF; fi

  cat ${expdir}/work/*_totempAffine.txt \
    | grep '^Parameters:' | awk '\
    BEGIN { for(i=0;i<12;i++) x[i]=0} \
    { for(i=0;i<12;i++) x[i]=x[i]+$(i+2) } \
    END { \
      printf "Transform: MatrixOffsetTransformBase_double_3_3\nParameters: "; \
      for(i=0;i<12;i++) printf "%f ",x[i]/NR; \
      printf "\nFixedParameters: 0 0 0\n";}' > $TEMPAFF

  # Compose the warps
  WarpImageMultiTransform 3 $TEMPWARP $TEMPWARP -i $TEMPAFF \
    -R ${expdir}/work/template_tse.nii.gz

  TEMPWARPFULL=${expdir}/template_fullwarp.nii.gz
  $ANTSPATH/ComposeMultiTransform 3 \
    $TEMPWARPFULL -R ${expdir}/work/template_tse.nii.gz \
    -i $TEMPAFF $TEMPWARP $TEMPWARP $TEMPWARP $TEMPWARP

  # Apply this warp to all the template derivatives
  for kind in $KINDS; do

    WarpImageMultiTransform 3 \
      ${expdir}/work/template_${kind}.nii.gz \
      ${expdir}/work/template_${kind}.nii.gz \
      $TEMPWARPFULL

  done
}

function ants_iter()
{
  id=$1
  side=$2
  doants=$3

  # Before we vote, use ml_affine for nice affine alignment
  ml_affine \
    ${expdir}/work/template_seg.nii.gz \
    ${expdir}/data/${id}_${side}_seg.nii.gz \
    ${expdir}/work/${id}_${side}_mlaffine.txt

  # Do we do ANTS or not?
  if [[ $doants -eq 0 ]]; then

    for sub in $KINDS; do

      c3d \
        ${expdir}/work/template_tse.nii.gz \
        ${expdir}/data/${id}_${side}_${sub}.nii.gz \
        -reslice-matrix \
        ${expdir}/work/${id}_${side}_mlaffine.txt \
        -o ${expdir}/work/${id}_${side}_totemp_reslice_${sub}.nii.gz

    done

  else

    # Convert that to ITK format
    c3d_affine_tool \
      ${expdir}/work/${id}_${side}_mlaffine.txt \
      -oitk ${expdir}/work/${id}_${side}_totempAffine.txt

    #WGT=$(echo ${#REG_LABELS[*]} | awk '{print 1.0 / $1}')

    CMD=""
    for sub in ${LABEL_IDS[*]}; do
      CMD="$CMD -w 1 -i $expdir/work/template_${sub}.nii.gz ${expdir}/data/${id}_${side}_${sub}.nii.gz "
    done

    c3d $expdir/work/template_${sub}.nii.gz -dup \
      ${expdir}/data/${id}_${side}_${sub}.nii.gz \
      -int 0 -reslice-identity \
      -add -binarize -dilate 1 10x10x10vox \
      -o $TMPDIR/${id}_${side}_totemp_mask.nii.gz

    greedy -d 3 \
      $CMD \
      -it ${expdir}/work/${id}_${side}_mlaffine.txt \
      -n 120x120x40 \
      -gm $TMPDIR/${id}_${side}_totemp_mask.nii.gz \
      -s 0.6mm 0.1mm \
      -e 0.5 \
      -o ${expdir}/work/${id}_${side}_totempWarp.nii.gz

    
<< COMMENT
    CMD=""
    for sub in ${REG_LABELS[*]}; do
      CMD="$CMD -m MSQ[${expdir}/work/template_${sub}.nii.gz,${expdir}/data/${id}_${side}_${sub}.nii.gz,$WGT]"
    done

    if [[ $ANTs_x == "Y" ]]; then
       ANTs_mask="-x ${expdir}/work/template_mask.nii.gz"
    else
       ANTs_mask=""
    fi

    ANTS 3 $CMD \
      -t $ANTs_t \
      -r $ANTs_r \
      -i $ANTs_i \
      $ANTs_mask \
      -a ${expdir}/work/${id}_${side}_mlaffine_itk.txt \
      --continue-affine 0 \
      $ANTs_all_metrics \
      -o ${expdir}/work/${id}_${side}_totemp.nii.gz \
      | tee ${expdir}/work/${id}_${side}_antsoutput.txt
COMMENT

    for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        ${expdir}/data/${id}_${side}_${sub}.nii.gz \
        ${expdir}/work/${id}_${side}_totemp_reslice_${sub}.nii.gz \
        -R ${expdir}/work/template_tse.nii.gz \
        ${expdir}/work/${id}_${side}_totempWarp.nii.gz \
        ${expdir}/work/${id}_${side}_totempAffine.txt

    done

  fi
}

function main_loop
{
  # Main iteration loop
  PREFIX=ANTs${expid}
  for ((iter=0;iter<$ITER;iter++)); do

    # Create the segmentation for the template
    c3d $(for sub in ${LABEL_IDS[*]}; do echo ${expdir}/work/template_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o ${expdir}/work/template_seg.nii.gz

    # Back up template
    ITDIR=${expdir}/work/$(printf iter_%02d $iter)
    mkdir -p $ITDIR
    cp -a ${expdir}/work/template_*.nii.gz $ITDIR/

    # Do ants?
    if [[ iter -lt $ANTs_start ]]; then doants=0; else doants=1; fi

    IDS=$(ls ${expdir}/data | grep _left_tse | sed -e "s/_left_tse.nii.gz//")

    # Run ANTS for each image
    for id in $IDS; do
      #Submit ANTS job
      for side in left right; do
        qsub -cwd -o $DUMPDIR -j y -N \
             "${PREFIX}_${id}_${side}" $0 \
             ants_iter $id $side $doants
      done
    done

    # Wait for completion
    qsub -cwd -o $DUMPDIR -j y \
      -hold_jid "${PREFIX}_*" -sync y -b y \
      sleep 1

    # If this is the last iteration, we don't want to recompute the template
    if [[ $iter -lt $((ITER-1)) ]]; then

      # Compute average images
      for kind in $KINDS; do

        if [[ $kind = "tse" ]]; then NORM=1; else NORM=0; fi
        AverageImages 3 ${expdir}/work/template_${kind}.nii.gz $NORM ${expdir}/work/*_totemp_reslice_${kind}.nii.gz

      done

      # Perform shape averaging
      if [[ $doants -eq 1 ]]; then
        shape_update_to_template "NOINPUT"
      fi

    fi
  done
}

##############################################################################
# Make images for further analysis
function make_images()
{
  # Generate one compound label for all other subfields
  ITDIR=${expdir}/work/$(printf iter_%02d $((ITER-1)) )
  for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

    LABELDEF=(${EXTRAMESHESDEF[i]})
    echo $LABELDEF

    c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "${expdir}/work/template_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
        -accum -add -endaccum \
        -o $ITDIR/template_${EXTRAMESHES[i]}.nii.gz

  done
  
  # for each subject
  IDS=$(ls ${expdir}/data | grep left_tse \
        | sed -e "s/_left_tse.nii.gz//")

  for side in left right; do
    for id in $IDS; do

    # Create seg for each subject
    c3d $(for sub in ${LABEL_IDS[*]}; do echo ${expdir}/work/${id}_${side}_totemp_reslice_${sub}.nii.gz; done) -vote -type ushort \
      -o ${expdir}/work/${id}_${side}_totemp_reslice_seg.nii.gz

    done
  done
}

##############################################################################
# Measurement dice overlap when warp template back to original subject space
# Warp labels back to subject space and vote
function warp_labels()
{
  PREFIX=WL${expid}
  mkdir -p ${expdir}/labelwarp

  # Iterate over side and subject
  IDS=$(ls ${expdir}/data | grep left_tse | sed -e "s/_left_tse.nii.gz//")
  for side in left right; do

    for id in $IDS; do

      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
        -N "${PREFIX}_${id}_${side}" \
        $0 warp_labels_subj $id $side

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
    -hold_jid "${PREFIX}_*" -sync y -b y \
    sleep 1
}

function warp_labels_subj()
{
  id=$1
  side=$2

  # add on transform
  if [[ $side == "left" ]]; then
    ADD_ON_TRANS=$IDENTITYFN
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
  fi

  # Apply the warp to each of the labels
  for sub in ${LABEL_IDS[*]}; do


    if [[ -f $ASHSRUNDIR/$fn/flirt_t2_to_t1/greedy_t2_to_t1_chunk_${side}_inv.mat ]]; then
    greedy -d 3 \
      -rf $ASHSRUNDIR/${id}/seg_${side}.nii.gz \
      -rm $expdir/work/template_${sub}.nii.gz \
          ${expdir}/labelwarp/${id}_${side}__${sub}_tempfit.nii.gz \
      -r $ASHSRUNDIR/$id/flirt_t2_to_t1/greedy_t2_to_t1_chunk_${side}_inv.mat,-1 \
         $ASHSRUNDIR/$id/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat,-1 \
         $ASHSRUNDIR/$id/affine_t1_to_template/t1_to_template_affine.mat,-1 \
         $ASHSRUNDIR/$id/ants_t1_to_temp/greedy_t1_to_template_invwarp.nii.gz \
         $ADD_ON_TRANS,-1 \
         ${expdir}/work/${id}_${side}_totempAffine.txt,-1 \
         ${expdir}/work/${id}_${side}_totempInverseWarp.nii.gz 

    else

    greedy -d 3 \
      -rf $ASHSRUNDIR/${id}/seg_${side}.nii.gz \
      -rm $expdir/work//template_${sub}.nii.gz \
          ${expdir}/labelwarp/${id}_${side}__${sub}_tempfit.nii.gz \
      -r $ASHSRUNDIR/$id/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat,-1 \
         $ASHSRUNDIR/$id/affine_t1_to_template/t1_to_template_affine.mat,-1 \
         $ASHSRUNDIR/$id/ants_t1_to_temp/greedy_t1_to_template_invwarp.nii.gz \
         $ADD_ON_TRANS,-1 \
         ${expdir}/work/${id}_${side}_totempAffine.txt,-1 \
         ${expdir}/work/${id}_${side}_totempInverseWarp.nii.gz

    fi


  done

  # Vote in original subject space to get the tempfit segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo ${expdir}/labelwarp/${id}_${side}_${sub}_tempfit.nii.gz; done) -vote -type ushort \
    -o ${expdir}/labelwarp/${id}_${side}_seg_tempfit.nii.gz


  # change ASHS labels to the current formate
  for ((i=0;i<${#LABEL_IDS[*]};i++)); do
    c3d $ASHSRUNDIR/$id/seg_${side}.nii.gz \
      -replace $(for k in ${LABEL_MRG[i]}; do echo $k 999; done) \
      -thresh 999 999 1 0 \
      -o ${expdir}/labelwarp/${id}_${side}_${LABEL_IDS[i]}_ashs.nii.gz
  done

  c3d $(for sub in ${LABEL_IDS[*]}; do echo ${expdir}/labelwarp/${id}_${side}_${sub}_ashs.nii.gz; done) -vote -type ushort \
    -o ${expdir}/labelwarp/${id}_${side}_seg_ashs.nii.gz
}

#########################################################
# Extract meshes for the subfields and apply warps to these meshes. This allows us to# perform statistical analysis on the mesh boundaries
function warp_meshes()
{
  PREFIX=WM${expid}
  mkdir -p ${expdir}/meshwarp
  mkdir -p ${expdir}/jacobian
  mkdir -p ${expdir}/meshwarp/template

  # Generate meshes for the individual subfields
  for sub in ${LABEL_FG[*]}; do

    vtklevelset \
      ${expdir}/work/template_${sub}.nii.gz \
      ${expdir}/meshwarp/template_${sub}.vtk \
      $subfield_th

  done

  # Generate one mesh for all non-DG non-CS subfields
  ITDIR=${expdir}/work/$(printf iter_%02d $((ITER-1)) )
  for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

    vtklevelset \
      $ITDIR/template_${EXTRAMESHES[i]}.nii.gz \
      ${expdir}/meshwarp/template_${EXTRAMESHES[i]}.vtk \
      $template_th

  done

  IDS=$(ls ${expdir}/data | grep left_tse | sed -e "s/_left_tse.nii.gz//")

  for id in ${IDS[*]}; do
    for side in left right; do

    # Submit job for this subject
    qsub  -cwd -o $DUMPDIR -j y \
        -q all.q \
        -N "${PREFIX}_${id}_${side}" \
        $0 warp_meshes_subj $id $side

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  cp ${expdir}/meshwarp/template*.vtk \
     ${expdir}/meshwarp/template/
}

function warp_meshes_subj()
{
  id=$1
  side=$2

  # add on transform
  if [[ $side == "left" ]]; then
    ADD_ON_TRANS=$IDENTITYFN
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
  fi

  # Compose the transformation between the template and the subject
  greedy -d 3 \
    -rf ${expdir}/work/template_tse.nii.gz \
    -rc $TMPDIR/compose.nii \
    -r ${expdir}/work/${id}_${side}_totempWarp.nii.gz \
       ${expdir}/work/${id}_${side}_totempAffine.txt \
       $ADD_ON_TRANS \
       $ASHSRUNDIR/$id/ants_t1_to_temp/greedy_t1_to_template_warp.nii.gz \
       $ASHSRUNDIR/$id/affine_t1_to_template/t1_to_template_affine.mat \
       $ASHSRUNDIR/$id/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat \
       $ASHSRUNDIR/$id/flirt_t2_to_t1/greedy_t2_to_t1_chunk_${side}_inv.mat

  # Split the transformations into xyz
  c3d -mcs $TMPDIR/compose.nii -oo $TMPDIR/comp%02d.nii

  # Apply the warp to each of the meshes
  for sub in $ALLSF; do

    # Warp the subfield into subject space
    warpmesh -w ants -m ras \
      $expdir/meshwarp/template_${sub}.vtk \
      $expdir/meshwarp/${id}_${side}_${sub}_tempfit.vtk \
      $TMPDIR/comp00.nii $TMPDIR/comp01.nii $TMPDIR/comp02.nii

    # Extract the thickness of the subfield
    cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
      -T $expdir/meshwarp/${id}_${side}_${sub}_thickmap.vtk \
      -p $thick_p -e $thick_e \
      $expdir/meshwarp/${id}_${side}_${sub}_tempfit.vtk \
      $expdir/meshwarp/skel_${id}_${side}_${sub}.vtk

  done
}

##############################################################################
# Evaluate the performance
# compute dice of mask
# compute dice of none DG none CS overlap
# compute mesh difference

function evaluation()
{
  EVALTMPDIR=$EVALDIR/tmp
  rm -rf $EVALTMPDIR
  mkdir -p $EVALTMPDIR

  PREFIX=EV${expid}
  IDS=$(ls ${expdir}/data | grep left_tse | sed -e "s/_left_tse.nii.gz//")
  for id in ${IDS[*]}; do
  for side in left right; do

      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
        -q all.q \
        -N "${PREFIX}_${id}_${side}" \
        $0 eval_subj $id $side $EVALTMPDIR
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # Combine scores
  IDs=($(cat $SUBJ_TXT))
  echo "ID Side ${EVALLABELS[*]}" | \
    sed 's/ /,/g'  \
    > $EVALDIR/overlap.csv

  for side in left right; do
    for id in ${IDs[*]}; do

      cat $EVALTMPDIR/${id}_${side}_overlap.txt \
          >> $EVALDIR/overlap.csv

    done
  done
}

function eval_subj()
{
  id=$1
  side=$2
  EVALTMPDIR=$3
  ALLOVL=""

  ###################################
  # Compute dice overlap in template space
  do_pair $expdir/work/template_seg.nii.gz \
    $expdir/work/${id}_${side}_totemp_reslice_seg.nii.gz
  ALLOVL="$ALLOVL $FULLOVL"

  ###################################
  # save to temporary file
  echo $id,$side$ALLOVL > $EVALTMPDIR/${id}_${side}_overlap.txt
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
}

##############################################################################
function gen_temp_label()
{
  MESHTEMPDIR=$expdir/meshwarp/template
  MESHLABELDIR=$MESHTEMPDIR/label
  rm -rf $MESHLABELDIR
  mkdir -p $MESHLABELDIR
  ITDIR=${expdir}/work/$(printf iter_%02d $((ITER-1)) )

  for ((ii=0;ii<${#EXTRAMESHES[*]};ii++)); do
    LABELDEF=(${EXTRAMESHESDEF[ii]})
    MESHES=""
    for ((i=0;i<${#LABEL_IDS[*]};i++)); do
      if [[ ${LABELDEF[$i]} == "1" ]]; then
        sub=${LABEL_IDS[$i]}
        mesh_image_sample \
          $MESHTEMPDIR/template_${EXTRAMESHES[ii]}.vtk \
          $ITDIR/template_${sub}.nii.gz \
          $MESHLABELDIR/template_${EXTRAMESHES[ii]}_${sub}.vtk \
          PROB
        MESHES="$MESHES $MESHLABELDIR/template_${EXTRAMESHES[ii]}_${sub}.vtk"
      fi
    done

    # merge prob arrays
    mesh_merge_arrays \
      -r $MESHTEMPDIR/template_${EXTRAMESHES[ii]}.vtk \
      $MESHLABELDIR/template_${EXTRAMESHES[ii]}_withlabel.vtk \
      PROB $MESHES

    # run matlab script to generate label
    $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
      addpath('/home/longxie/ASHS_PHC/thickness_newlabel/matlabcode');
      compute_mesh_label('$MESHLABELDIR/template_${EXTRAMESHES[ii]}_withlabel.vtk');
MATCODE

  done
}


##############################################################################
function warp_errormap()
{
  method=$1

  PREFIX=WEM${expid}
  mkdir -p ${expdir}/errormap_warp/$method

  # Iterate over side and subject
  IDS=$(ls ${expdir}/data | grep left_tse | sed -e "s/_left_tse.nii.gz//")
  for side in left right; do
    for id in $IDS; do

      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
        -N "${PREFIX}_${id}_${side}_${method}" \
        $0 warp_errormap_subj $id $side ${method}

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
    -hold_jid "${PREFIX}_*" -sync y -b y \
    sleep 1

  # average all the erromaps
  c3d ${expdir}/errormap_warp/$method/*/*errormap_${method}_totemp.nii.gz \
    -mean \
    -o ${expdir}/errormap_warp/$method/mean_errormap_${method}.nii.gz \
    -smooth 0.3x0.3x0.3vox \
    -o ${expdir}/errormap_warp/$method/mean_errormap_${method}_smoothed.nii.gz
}

function warp_errormap_subj()
{
  id=$1
  side=$2
  method=$3
  OUTDIR=${expdir}/errormap_warp/$method/${id}_${side}
  mkdir -p $OUTDIR
  LOGTXT=${expdir}/errormap_warp/$method/log.txt

  # compute errormap
  if [[ $method == "MV" ]]; then
    AUTOSEG=$AUTOSEGDIR/other_LFmethods/${id}_${side}_MV.nii.gz
  elif [[ $method == "SVWV" ]]; then
    AUTOSEG=$AUTOSEGDIR/other_LFmethods/${id}_${side}_GaussSVWV[0.05]_rp3_rs4.nii.gz
  elif [[ $method == "ashs_usegray" ]]; then
    AUTOSEG=$(ls /home/longxie/DeepLabelFusion/PMC_Prisma_3T_withflip/exp_ASHST1/exp000/testing_7TT2/$id/multiatlas/fusion/lfseg_corr_usegray_${side}.nii.gz)
  elif [[ $method == "ashs_heur" ]]; then
    AUTOSEG=$(ls /home/longxie/DeepLabelFusion/PMC_Prisma_3T_withflip/exp_ASHST1/exp000/testing_7TT2/$id/multiatlas/fusion/lfseg_heur_${side}.nii.gz)
  elif [[ $method == "UNET" ]]; then
    AUTOSEG=$AUTOSEGDIR/unet_run002/test_seg_${id}_${side}_native.nii.gz
  elif [[ $method == "DLF" ]]; then
    AUTOSEG=$AUTOSEGDIR/deeplabelfusion_run395/test_seg_${id}_${side}_native.nii.gz 
  #elif [[ $method == "ASHSUsegrayBootPaul" ]]; then
  #  AUTOSEG=$(ls /data/picsl/pauly/wolk/abc_prisma/exp01/xval/xval*/test/xval*orig_${id}/final/xval*_test_orig_${id}_${side}_lfseg_corr_usegray.nii.gz)
  #  if [[ $AUTOSEG == "" ]]; then
  #    echo "Can not find $id $side ASHSUsegrayBootPaul" >> $LOGTXT
  #    exit
  #  fi
  #elif [[ $method == "ASHSHeurBootPaul" ]]; then
  #  AUTOSEG=$(ls /data/picsl/pauly/wolk/abc_prisma/exp01/xval/xval*/test/xval*orig_${id}/final/xval*_test_orig_${id}_${side}_lfseg_heur.nii.gz)
  #  if [[ $AUTOSEG == "" ]]; then
  #    echo "Can not find $id $side ASHSUsegrayBootPaul" >> $LOGTXT
  #    exit
  #  fi
  #elif [[ $method == "ASHSUsegrayBootPaulExp02" ]]; then
  #  AUTOSEG=$(ls /data/picsl/pauly/wolk/abc_prisma/exp02/xval/xval*/test/xval*orig_${id}/final/xval*_test_orig_${id}_${side}_lfseg_corr_usegray.nii.gz)
  #  if [[ $AUTOSEG == "" ]]; then
  #    echo "Can not find $id $side ASHSUsegrayBootPaul" >> $LOGTXT
  #    exit
  #  fi
  else
    echo "does not support $method"
    exit
  fi
  REFSEG=$ASHSRUNDIR/$id/seg_${side}.nii.gz
  c3d $AUTOSEG -dup $REFSEG \
    -int 0 -reslice-identity \
    -o $OUTDIR/${id}_${side}_refseg.nii.gz \
    -scale -1 -add \
    -binarize \
    -o $OUTDIR/${id}_${side}_errormap_${method}.nii.gz 
  ln -sf $AUTOSEG $OUTDIR/${id}_${side}_${method}_autoseg.nii.gz

  # warp error map to template
  if [[ $side == "left" ]]; then
    ADD_ON_TRANS=$IDENTITYFN
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
  fi

  # Compose the transformation between the template and the subject
  if [[ -f $ASHSRUNDIR/$id/flirt_t2_to_t1/greedy_t2_to_t1_chunk_${side}_inv.mat ]]; then

  greedy -d 3 \
    -rf ${expdir}/work/template_tse.nii.gz \
    -rm $OUTDIR/${id}_${side}_errormap_${method}.nii.gz \
        $OUTDIR/${id}_${side}_errormap_${method}_totemp.nii.gz \
    -r ${expdir}/work/${id}_${side}_totempWarp.nii.gz \
       ${expdir}/work/${id}_${side}_totempAffine.txt \
       $ADD_ON_TRANS \
       $ASHSRUNDIR/$id/ants_t1_to_temp/greedy_t1_to_template_warp.nii.gz \
       $ASHSRUNDIR/$id/affine_t1_to_template/t1_to_template_affine.mat \
       $ASHSRUNDIR/$id/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat \
       $ASHSRUNDIR/$id/flirt_t2_to_t1/greedy_t2_to_t1_chunk_${side}_inv.mat

  else

  greedy -d 3 \
    -rf ${expdir}/work/template_tse.nii.gz \
    -rm $OUTDIR/${id}_${side}_errormap_${method}.nii.gz \
        $OUTDIR/${id}_${side}_errormap_${method}_totemp.nii.gz \
    -r ${expdir}/work/${id}_${side}_totempWarp.nii.gz \
       ${expdir}/work/${id}_${side}_totempAffine.txt \
       $ADD_ON_TRANS \
       $ASHSRUNDIR/$id/ants_t1_to_temp/greedy_t1_to_template_warp.nii.gz \
       $ASHSRUNDIR/$id/affine_t1_to_template/t1_to_template_affine.mat \
       $ASHSRUNDIR/$id/flirt_t2_to_t1/flirt_t2_to_t1_inv.mat

  fi
}

##################################################
function sample_errormap()
{
  method=$1
  for sub in $ALLSF; do
      ERRORMESHDIR=${expdir}/errormap_warp/meshes/$method
      mkdir -p $ERRORMESHDIR

      cp $expdir/meshwarp/template_${sub}.vtk \
        $ERRORMESHDIR/${sub}_errormap_${method}.vtk


      METHODDIR=${expdir}/errormap_warp/$method
      mesh_image_sample \
        $ERRORMESHDIR/${sub}_errormap_${method}.vtk \
        ${expdir}/errormap_warp/$method/mean_errormap_${method}.nii.gz \
        $ERRORMESHDIR/${sub}_errormap_${method}.vtk \
        ${method}_orig

      mesh_image_sample \
        $ERRORMESHDIR/${sub}_errormap_${method}.vtk \
        ${expdir}/errormap_warp/$method/mean_errormap_${method}_smoothed.nii.gz \
        $ERRORMESHDIR/${sub}_errormap_${method}.vtk \
        ${method}_smoothed
  done
}

##################################################
function diff_errormap()
{
  method1=$1
  method2=$2
  for sub in $ALLSF; do
      ERRORMESHDIR=${expdir}/errormap_warp/meshes/diff_${method1}-${method2}
      mkdir -p $ERRORMESHDIR

      c3d ${expdir}/errormap_warp/$method1/mean_errormap_${method1}.nii.gz \
        ${expdir}/errormap_warp/$method2/mean_errormap_${method2}.nii.gz -scale -1 \
        -add \
        -o $ERRORMESHDIR/diff_mean_errormap_${method1}-${method2}.nii.gz

      cp ${expdir}/errormap_warp/meshes/${method1}/${sub}_errormap_${method1}.vtk \
        $ERRORMESHDIR/${sub}_diff_mean_errormap_${method1}-${method2}.vtk

      mesh_image_sample \
        $ERRORMESHDIR/${sub}_diff_mean_errormap_${method1}-${method2}.vtk \
        $ERRORMESHDIR/diff_mean_errormap_${method1}-${method2}.nii.gz \
        $ERRORMESHDIR/${sub}_diff_mean_errormap_${method1}-${method2}.vtk \
        ${method1}-${method2}_diff
  done
}

##################################################
function get_volume()
{
  mkdir -p $VOLUMEDIR
  IDs=($(cat $SUBJ_TXT))
  echo "ID Side ${VOLLABELS[*]}" | \
    sed 's/ /,/g'  \
    > $VOLUMEDIR/refseg_volume.csv

  for id in ${IDs[*]}; do
    for side in left right; do

      REFSEG=$TEMPLATEDIR/train/atlas/$id/seg_${side}.nii.gz
      do_single_volume $REFSEG
      echo "$id,$side$FULLVOL" >> $VOLUMEDIR/refseg_volume.csv
    done
  done
}

function do_single_volume()
{
  # Get a pair of segmentations
  seg=$1

  # out dice file
  #out_dice_file=$3

  # Iterate over all relevant labels
  FULLVOL=""
  for ((i=0; i<${#VOLLABELS[*]}; i++)); do

    # Do the analysis on full-size meshes
    REPRULE=$(for lab in ${VOLRANGES[i]}; do echo $lab 99; done)

    # Extract the binary images and compute overlap
    c3d \
      $seg -replace $REPRULE -thresh 99 99 1 0 \
      -dup -lstat | tee $TMPDIR/vol.txt

    # Get the full-extent overlap
    VOL=$(cat $TMPDIR/vol.txt | grep "    1    " | awk -F '[ ,]+' '{print $8}')

    #echo $id ${LABELS[i]} full $OVL $DIST >> $out_file
    FULLVOL="${FULLVOL},${VOL}"

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
#if [[ $cmd == $RUNCMD1 || $cmd == $RUNCMD2 ]]; then

  main $@

else
  cmd=$1
  shift
  $cmd $@
fi
