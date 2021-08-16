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
ASHSRUNDIR=$ROOT/ASHS_OTS/fullset/ashs
ASHSTEMPDIR=$ROOT/ASHS_OTS/headtailatlas/run02/final/template
GTSEGDIR=$ROOT/rawdata/atlas2014/segs_OTS/
MATLABCODEDIR=$WORKDIR/matlabcode/
GTGROUPDIR=$WORKDIR/group/
SUBJ_TXT=$WORKDIR/analysis_input/subj_xval_vertical.txt
SUBJ_WITHTEMP_TXT=$WORKDIR/analysis_input/subj_xval_vertical_withtemp.txt
FLIPXFN=$WORKDIR/analysis_input/flipx_itk.txt
IDENTITYFN=$WORKDIR/analysis_input/identity_itk.txt

##############################################################################
# Parameters needs to be specify
# Experiment number
expid=700
expdir=$WORKDIR/exp/exp${expid}
DUMPDIR=${expdir}/dump
mkdir -p ${expdir}/dump

########################################
# 0. copy data
CLEANUPDIR=$expdir/cleanup

LABEL_IDS_ALL=(BKG                ERC  BA35 BA36 PHC  aCS  pCS)
LABEL_MRG_ALL=("0 1 2 3 4 7 8 16" "10" "11" "12" "13" "15" "14")
# Labels to get segmentation
LABEL_IDS=(BKG                ERC  BA35 BA36 PHC  aCS  pCS)
LABEL_MRG=("0 1 2 3 4 7 8 16" "10" "11" "12" "13" "15" "14")
LABEL_NEW=(0                  1    2    3    4    5    6)
MESH_LABEL=(ERC  BA35 BA36 PHC)
KINDS="tse mprage ${LABEL_IDS[*]}"
DATADIR=${expdir}/data


########################################
# 1 pairwise registration and clustering
# 1.1 ANTs parameters
REUSE="N"
REUSEID="504"
REUSEEXPDIR=$WORKDIR/exp/exp${REUSEID}
WGT_1=1
REG_LABELS_1=(${LABEL_IDS[*]})
ANTs_t_1="SyN[0.25]"
ANTs_r_1="Gauss[1.5,0.5]"
ANTs_i_1="15x8x0"
ANTs_x_1="Y"
ANTs_G_1=""
ANTs_all_metrics_1="--use-all-metrics-for-convergence"
PWDIR=${expdir}/pairwise/

# 1.2 similarity type
SIM_TYPE="PRCCS_seg_dice"
SIM_DIR=$expdir/PWsim/sim_${SIM_TYPE}

# 1.3 clustering
INITGRP_DIR=$expdir/clustering/InitGroup/

# 1.4 prepare statitics documents
GPSTATDIR=$expdir/stat/groups/
ALLSTATDIR=$expdir/stat/all/

#########################################
# 2 generate three templates
# Relevant labels
LABEL_FG=(           ERC BA35 BA36 PHC aCS pCS)
EXTRAMESHES=(PRC Anterior MRG)
#                BKG ERC BA35 BA36 PHC aCS pCS
EXTRAMESHESDEF=("-1  -1   1    1   -1  -1  -1" \
                "-1   1   1    1   -1  -1  -1" \
                "-1   1   1    1    1  -1  -1")

# 2.1 initilization
groups=(1 2 3)
INITTEMP_DIR=$expdir/MultiTemps/

# 2.2 ANTs
ITER=8
ANTs_start=3
WGT_2=1
REG_LABELS_2=(${LABEL_IDS[*]})
ANTs_t_2="SyN[0.25]"
ANTs_r_2="Gauss[1.5,0.5]"
ANTs_i_2="80x80x20"
ANTs_x_2="Y"
ANTs_all_metrics_2="--use-all-metrics-for-convergence"

# 2.3 warp labels

# 2.4 warp meshes
subfield_th=0.5
template_th=0.0
thick_p=1.2
thick_e=6

# 2.6 evaluation
EVALLABELS=(ERC BA35 BA36 PRC   PHC aCS pCS ALL)
RANGES=(    1   2    3    "2 3" 4   5   6   "1 2 3 4")
MESH_EVAL=(       ERC BA35 BA36 PRC   PHC aCS pCS)
MESH_EVAL_RANGES=(1   2    3    "2 3" 4   5   6)
INITEVALDIR=$INITTEMP_DIR/evaluation

# 2.6 

# 2.7 thickness analysis
#ANGRP[0]="MRG Anterior"
ANGRP[0]="MRG"
ANGRP[1]="ERC BA35 BA36 PHC"

GRPNM[0]="merged"
GRPNM[1]="all"

# 2.8 Mesh labels


#########################################
# 3 insert three templates into the graph

# 3.3 
FINALGRAPHDIR=$expdir/clustering/finalgraph/

#########################################
# 4 build graph and generate super template
FINALTEMPDIR=$expdir/FinalTemp

# 4.1 
ALLPATHSDIR=$FINALTEMPDIR/paths/allpaths

# 4.3 
BESTPATHSDIR=$FINALTEMPDIR/paths/bestpaths

# 4.5 template propagation
WGT_4=1
#REG_LABELS_4=(${LABEL_IDS[*]})
REG_LABELS_4=(BKG                ERC  BA35 BA36 PHC)
ANTs_t_4="SyN[0.25]"
ANTs_r_4="Gauss[2,1]"
ANTs_i_4="40x20x0"
ANTs_x_4="Y"
ANTs_G_4=""
ANTs_SMOOTH_4=""
#ANTs_SMOOTH_4="--gaussian-smoothing-sigmas 3x2x1"
ANTs_all_metrics_4="--use-all-metrics-for-convergence"

groups_temp=(3)

# 4.6
FINALITER=6
ANTs_start_final=2
WGT_3=1
REG_LABELS_3=(${LABEL_IDS[*]})
ANTs_t_3="SyN[0.25]"
ANTs_r_3="Gauss[1.5,0.5]"
ANTs_i_3="80x80x20"
ANTs_x_3="Y"
ANTs_all_metrics_3="--use-all-metrics-for-convergence"

# exp1: ANTs_i_3="60x30x10"
# exp2: ANTs_i_3="200x150x100"
# exp3: ANTs_i_3="40x20x5", 8, 3

##############################################################################
function main()
{
  reset_dir

  #######################################
  # 0. preparation
  # perform manual grouping
  #ManualGrouping

  #PreprocessLabels
  #copy_data

  #######################################
  # 1 cluster subjects into three groups
  # 1.1 pairwise registration
  #pairwise
  #completeness

  # 1.2 prepare the group info
  #PrepGroup

  ########################################
  # 2 build three templates
  # 2.1 initialization
  #initialization

  # 2.2 construct three templates
  #main_loop

  # 2.3 make images
  #make_images

  # 2.4 wapr labels back to subject space
  #warp_labels

  # 2.5 warp meshes back to subject space and measure thickness
  #warp_meshes

  # 2.6 evaluation
  #evaluation

  # 2.7 displacement
  #disp_stats

  # 2.8 within template statistical analysis
  #thick_stats

  # 2.9 generate labels for each template mesh
  #mean_thickness

  ########################################
  # 3 insert three templates into the graph
  # 3.1 perform pairwise registration
  #Insertion

  # 3.2 measure similarity between templates and subjects
  #TempSim

  ########################################
  # 4 build graph and generate super template


  # 4.1 find paths between varient templates
  #FindPath

  # 4.4 propagate three templates to build super template
  #PropagateTemplatesNew

  # 4.5 warp all the subjects to the super template space
  WarpSuperTemplate
  AverageFinalTemplate

  # 4.6 register to the final super template
  RegisterSuperTemplate

  # 4.6 make final images
  MakeFinalImages

  # 4.7 warp super template labels
  WarpFinalLabel
  #WarpFinalLabelToData

  # 4.8 warp super template meshes
  WarpFinalMesh
  
  # 4.9 evaluate super template
  EvalFinalTemp
  #EvalFinalTempData

  # 4.10 displacement
  DispFinalStat

  # 4.11 thickness analysis
  #ThickFinalStat

  # 4.12 generate labels for each template mesh
  #FinalMeanThickness

}

##########################################################
function ManualGrouping()
{
  # Load id number
  IDs=($(cat $SUBJ_TXT))

  # Submit job to copy data
  PREFIX=PL${expid}
  for side in left right; do
    for ((i=0;i<${#IDs[*]};i++)); do

      id=${IDs[i]}

      itksnap -g  $ROOT/rawdata/atlas2014/nii/${id}_tse_TP1.nii.gz \
        -s $GTSEGDIR/${id}_${side}_sf_full.nii.gz


    done
  done
}

##########################################################
function PreprocessLabels()
{
  mkdir -p $CLEANUPDIR

  # Load id number
  IDs=($(cat $SUBJ_TXT))

  # Submit job to copy data
  PREFIX=PL${expid}
  for side in left right; do
    for ((i=0;i<${#IDs[*]};i++)); do

      id=${IDs[i]}

      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${id}_${side}" \
           $0 PreprocessLabels_sub $id $side
      sleep 0.1

     done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function PreprocessLabels_sub()
{
  id=$1
  side=$2

  in=$GTSEGDIR/${id}_${side}_sf_full.nii.gz
  out5=$TMPDIR/${id}_${side}_lfseg_corr_nogray_cleanup.nii.gz
  out=$CLEANUPDIR/${id}_${side}_dividedCS.nii.gz

  ## fifth step: clean CS
  $C3DPATH/c3d $in -as S \
    -thresh 10 13 1 0 -dilate 1 100x100x0 \
    -push S -thresh 14 14 1 0 -multiply -popas CSSEL \
    -push S -thresh 14 14 1 0 -push CSSEL -scale -1 -add \
    -scale -1 -shift 1 -push S -multiply \
    -o $out5

  $C3DPATH/c3d $out5 -as S \
      -thresh 10 12 1 0 -dilate 1 100x100x0 \
      -push S -times \
      -thresh 14 14 1 0 \
      -push S -add \
      -o $out
}

###########################################################
# Copy data
function copy_data()
{
  # perform initial affine registration between left and flip right temps
  mkdir -p $DATADIR/tempreg
  ln -sf $ASHSTEMPDIR/refspace_meanseg_left.nii.gz \
     $DATADIR/tempreg/refspace_meanseg_left.nii.gz
  c3d $ASHSTEMPDIR/refspace_meanseg_right.nii.gz -flip x \
    -o $DATADIR/tempreg/refspace_meanseg_right_flip.nii.gz

  # register fliped right to left
  ml_affine $DATADIR/tempreg/refspace_meanseg_left.nii.gz \
    $DATADIR/tempreg/refspace_meanseg_right_flip.nii.gz \
    $DATADIR/tempreg/right_to_left_mlaffine.txt
  c3d_affine_tool $DATADIR/tempreg/right_to_left_mlaffine.txt \
    -oitk $DATADIR/tempreg/right_to_left_mlaffine_itk.txt

  # combine affine transformations
  $ANTSPATH/ComposeMultiTransform 3 \
    $DATADIR/tempreg/right_to_left_init_affine.txt \
    -R $DATADIR/tempreg/right_to_left_mlaffine_itk.txt \
    $DATADIR/tempreg/right_to_left_mlaffine_itk.txt \
    $FLIPXFN

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
           -q all.q,basic.q \
           -l h_vmem=2.1G,s_vmem=2G \
           -N "${PREFIX}_${id}_${side}" \
           $0 copy_subject $id $side
      sleep 0.1

     done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function copy_subject()
{
  id=$1
  side=$2
  fn=$(ls $ASHSRUNDIR | grep $id)

  # ASHS segmentation
  SEG=$CLEANUPDIR/${id}_${side}_dividedCS.nii.gz
  if [[ ! -f $SEG ]]; then
    echo "$SEG does not exist."
    exit
  fi

  # Link the subfield images
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
      -R $ASHSRUNDIR/$fn/tse_to_chunktemp_left.nii.gz \
      $ADD_ON_TRANS
    WarpImageMultiTransform 3  \
      $ASHSRUNDIR/$fn/mprage_to_chunktemp_${side}.nii.gz \
      $DATADIR/${fn}_${side}_mprage.nii.gz \
      -R $ASHSRUNDIR/$fn/mprage_to_chunktemp_left.nii.gz \
      $ADD_ON_TRANS
  fi

  # Generate a binary image for the label
  for ((i=0;i<${#LABEL_IDS_ALL[*]};i++)); do

    c3d $SEG -replace $(for k in ${LABEL_MRG_ALL[i]}; do echo $k 999; done) -thresh 999 999 1 0 \
      -o $TMPDIR/binary_${LABEL_IDS_ALL[i]}_${fn}_${side}.nii.gz

    WarpImageMultiTransform 3  \
      $TMPDIR/binary_${LABEL_IDS_ALL[i]}_${fn}_${side}.nii.gz \
      $DATADIR/${fn}_${side}_${LABEL_IDS_ALL[i]}.nii.gz \
      -R $DATADIR/${fn}_${side}_tse.nii.gz \
      $ADD_ON_TRANS \
      $ASHSRUNDIR/$fn/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz \
      $ASHSRUNDIR/$fn/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
      $ASHSRUNDIR/$fn/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

  done

  # Vote in the ASHS template subject space to create segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $DATADIR/${fn}_${side}_${sub}.nii.gz; done) \
    -vote -type ushort -o $DATADIR/${fn}_${side}_seg.nii.gz
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
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${id_fix}_${side}" $0 \
           pairwise_sub $id_fix $side
      sleep 0.1

    done

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function pairwise_sub()
{
  id_fix=$1
  side_fix=$2
  idside_fix=${id_fix}_${side_fix}
  IDS=$(cat $SUBJ_TXT)
  fn_fix=$(ls $ASHSRUNDIR | grep $id_fix)
  fnside_fix=${fn_fix}_${side_fix}

  # preprocess for similarity measure
  fn_seg_fix=$TMPDIR/${fnside_fix}_seg_tmp.nii.gz
  c3d $DATADIR/${fnside_fix}_seg.nii.gz \
    -replace 1 0 4 0 6 0 \
    -o $fn_seg_fix

  for id_mov in $IDS; do
  for side_mov in left right; do

    idside_mov=${id_mov}_${side_mov}
    if [[ $idside_mov != $idside_fix ]]; then

      fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)
      fnside_mov=${fn_mov}_${side_mov}
      OUTDIR=$PWDIR/${fnside_fix}/${fnside_mov}_to_${fnside_fix}
      #TMPDIR=$OUTDIR
      mkdir -p $PWDIR/${fnside_fix}

      # perform registration
      if [[ -f $OUTDIR/${fnside_mov}_to_${fnside_fix}_sim.txt ]]; then

        echo "Seg file exists."

      else

        mkdir -p $OUTDIR

        # Use ml_affine for nice affine alignment
        /data/picsl/pauly/wolk/ashs/ext/Linux/bin/ml_affine \
          $DATADIR/${fnside_fix}_seg.nii.gz \
          $DATADIR/${fnside_mov}_seg.nii.gz \
          $TMPDIR/${fnside_mov}_to_${fnside_fix}_mlaffine.txt

        # Convert that to ITK format
        c3d_affine_tool \
          $TMPDIR/${fnside_mov}_to_${fnside_fix}_mlaffine.txt \
          -oitk \
          $TMPDIR/${fnside_mov}_to_${fnside_fix}_mlaffine_itk.txt

        CMD=""
        for sub in ${REG_LABELS_1[*]}; do
          CMD="$CMD -m MSQ[$DATADIR/${fnside_fix}_${sub}.nii.gz,$DATADIR/${fnside_mov}_${sub}.nii.gz,$WGT_1]"
        done

        if [[ $ANTs_x_1 == "Y" ]]; then
          c3d $DATADIR/${fnside_fix}_seg.nii.gz -dup \
            $DATADIR/${fnside_mov}_seg.nii.gz \
            -int 0 -reslice-identity \
            -add -binarize -dilate 1 10x10x10vox \
            -o $TMPDIR/${fnside_mov}_to_${fnside_fix}_mask.nii.gz
          ANTs_mask_1="-x $TMPDIR/${fnside_mov}_to_${fnside_fix}_mask.nii.gz"
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
           -a $TMPDIR/${fnside_mov}_to_${fnside_fix}_mlaffine_itk.txt \
           --continue-affine 0 \
           $ANTs_all_metrics_1 \
           -o $TMPDIR/${fnside_mov}_to_${fnside_fix}_pairwise.nii.gz

        for sub in ${KINDS[*]}; do

          WarpImageMultiTransform 3 \
            $DATADIR/${fnside_mov}_${sub}.nii.gz \
            $TMPDIR/${fnside_mov}_to_${fnside_fix}_reslice_${sub}.nii.gz \
            -R $DATADIR/${fnside_fix}_tse.nii.gz \
            $TMPDIR/${fnside_mov}_to_${fnside_fix}_pairwiseWarp.nii.gz \
            $TMPDIR/${fnside_mov}_to_${fnside_fix}_pairwiseAffine.txt

        done

        # Create seg
        c3d $(for sub in ${LABEL_IDS[*]}; do echo $TMPDIR/${fnside_mov}_to_${fnside_fix}_reslice_${sub}.nii.gz; done) \
          -vote -type ushort \
          -o $TMPDIR/${fnside_mov}_to_${fnside_fix}_reslice_seg.nii.gz

        # measure similarity
        OVL=$(c3d $TMPDIR/${fnside_mov}_to_${fnside_fix}_reslice_seg.nii.gz \
          -replace 1 0 4 0 6 0 \
          $fn_seg_fix -label-overlap \
          | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}' )
        echo "$OVL" > \
          $OUTDIR/${fnside_mov}_to_${fnside_fix}_sim.txt


      fi
    fi

  done
  done
}

##############################################################################
function completeness
{
  rm -rf $PWDIR/check_pairwise.txt

  for side_fix in left right; do

    IDS=$(ls $DATADIR | grep ${side_fix}_tse | sed -e "s/_${side_fix}_tse.nii.gz//")

    for id_fix in $IDS; do

      for side_mov in left right; do
      for id_mov in $IDS; do

      idside_fix=${id_fix}_${side_fix}
      idside_mov=${id_mov}_${side_mov}

      if [[ $idside_fix != $idside_mov ]]; then

      # Check whether all the files exist
      OUTDIR=$PWDIR/${idside_fix}/${idside_mov}_to_${idside_fix}
      if [[ -f $OUTDIR/${idside_mov}_to_${idside_fix}_sim.txt ]]; then
          echo "ok"
      else
        echo "${idside_mov} to ${idside_fix} of ${side} is missing" \
          >> $PWDIR/check_pairwise.txt
      fi

      fi
      done
      done
    done
  done

  if [ -f $PWDIR/check_pairwise.txt ]; then
    echo "pairwise registration has missing data"
    exit
  fi
}




###########################################################
function PrepGroup()
{
  rm -rf $INITGRP_DIR
  mkdir -p $INITGRP_DIR

  # get IDs and other info
  IDS=($(cat $SUBJ_TXT))
  LGRP=($(cat $GTGROUPDIR/group_xval_left.txt))
  RGRP=($(cat $GTGROUPDIR/group_xval_right.txt))

  # prepare the files
  for side in left right; do
    for ((i=0;i<${#IDS[*]};i++)); do

      id=${IDS[i]}
      if [[ $side == "left" ]]; then
        grp=${LGRP[i]}
      else
        grp=${RGRP[i]}
      fi

      # output files
      echo ${id}_${side} >> $INITGRP_DIR/subjID_${grp}.txt
      echo ${grp} >> $INITGRP_DIR/group.txt

    done
  done
}

###########################################################
# Take the one that is most similar to others as the initial template
function average_subfield()
{
  grp=$1
  kind=$2

  GRP_IDS=($(cat $INITGRP_DIR/subjID_${grp}.txt))  
  for ((i=0;i<${#GRP_IDS[*]};i++)); do
    IDS[$i]=$(echo ${GRP_IDS[i]} | cut -d '_' -f 1)
    SIDES[$i]=$(echo ${GRP_IDS[i]} | cut -d '_' -f 2)
  done

  AverageImages 3 $INITTEMP_DIR/group${grp}/work/template_fullchunk_${grp}_${kind}.nii.gz 0 $(for ((i=0;i<${#GRP_IDS[*]};i++)); do echo $DATADIR/*${IDS[i]}*${SIDES[i]}_${kind}.nii.gz; done)
}

function initialization()
{
  # Initialization for all the groups
  for grp in ${groups[*]}; do

    GROUPDIR=$INITTEMP_DIR/group${grp}
    mkdir -p ${GROUPDIR}/work

    # Compute initial average for each subfield mask to generate mask
    PREFIX=AVG${expid}
    # Average all the input images
    for kind in $KINDS; do
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -N "${PREFIX}_${kind}_${grp}" $0 \
           average_subfield $grp $kind
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # Compute the mask for all the groups and cp init template
  for grp in ${groups[*]}; do

    # Compute the initial mask for the segmentation
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITTEMP_DIR/group${grp}/work/template_fullchunk_${grp}_${sub}.nii.gz; done | grep -v BKG) \
      -mean -thresh 1e-5 inf 1 0 -trim 5vox \
      -o $INITTEMP_DIR/group${grp}/work/template_mask_${grp}.nii.gz

    # Trim every template component using the mask
    for kind in $KINDS; do

      c3d $INITTEMP_DIR/group${grp}/work/template_mask_${grp}.nii.gz \
        $INITTEMP_DIR/group${grp}/work/template_fullchunk_${grp}_${kind}.nii.gz \
        -reslice-identity \
        -o $INITTEMP_DIR/group${grp}/work/template_${grp}_${kind}.nii.gz

    done
  done
}

##############################################################################
# Build template using ANTs
function shape_update_to_template()
{

  grp=$1

  # Borrowed from ANTS buildtemplateparallel

  # Average the warp fields
  local TEMPWARP=$INITTEMP_DIR/group${grp}/work/template_${grp}_warp.nii.gz
  $NEWANTSDIR/AverageImages 3 $TEMPWARP 0 \
    $INITTEMP_DIR/group${grp}/work/*_${grp}_totempWarp.nii.gz

  # Scale by -0.25 (gradient step)
  MultiplyImages 3 $TEMPWARP -0.25 $TEMPWARP

  # Create a new template affine
  local TEMPAFF=$INITTEMP_DIR/group${grp}/work/template_${grp}_Affine.txt
  if [[ -f $TEMPAFF ]]; then rm -rf $TEMPAFF; fi

  cat $INITTEMP_DIR/group${grp}/work/*_${grp}_totempAffine.txt \
    | grep '^Parameters:' | awk '\
    BEGIN { for(i=0;i<12;i++) x[i]=0} \
    { for(i=0;i<12;i++) x[i]=x[i]+$(i+2) } \
    END { \
      printf "Transform: MatrixOffsetTransformBase_double_3_3\nParameters: "; \
      for(i=0;i<12;i++) printf "%f ",x[i]/NR; \
      printf "\nFixedParameters: 0 0 0\n";}' > $TEMPAFF

  # Compose the warps
  WarpImageMultiTransform 3 $TEMPWARP $TEMPWARP -i $TEMPAFF \
    -R $INITTEMP_DIR/group${grp}/work/template_${grp}_tse.nii.gz

  TEMPWARPFULL=$INITTEMP_DIR/group${grp}/work/template_${grp}_fullwarp.nii.gz
  $ANTSPATH/ComposeMultiTransform 3 \
    $TEMPWARPFULL -R $INITTEMP_DIR/group${grp}/work/template_${grp}_tse.nii.gz \
    -i $TEMPAFF $TEMPWARP $TEMPWARP $TEMPWARP $TEMPWARP

  # Apply this warp to all the template derivatives
  for kind in $KINDS; do

    WarpImageMultiTransform 3 \
      $INITTEMP_DIR/group${grp}/work/template_${grp}_${kind}.nii.gz \
      $INITTEMP_DIR/group${grp}/work/template_${grp}_${kind}.nii.gz \
      $TEMPWARPFULL

  done
}

function ants_iter()
{
  id=$1
  side=$2
  grp=$3
  doants=$4

  # Before we vote, use ml_affine for nice affine alignment
  ~pauly/wolk/ashs/ext/Linux/bin/ml_affine \
    $INITTEMP_DIR/group${grp}/work/template_${grp}_seg.nii.gz \
    $DATADIR/${id}_${side}_seg.nii.gz \
    $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_mlaffine.txt

  # Do we do ANTS or not?
  if [[ $doants -eq 0 ]]; then

    for sub in $KINDS; do

      c3d \
        $INITTEMP_DIR/group${grp}/work/template_${grp}_tse.nii.gz \
        $DATADIR/${id}_${side}_${sub}.nii.gz \
        -reslice-matrix \
        $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_mlaffine.txt \
        -o $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp_reslice_${sub}.nii.gz

    done

  else

    # Convert that to ITK format
    c3d_affine_tool \
      $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_mlaffine.txt \
      -oitk $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_mlaffine_itk.txt

    #WGT=$(echo ${#REG_LABELS[*]} | awk '{print 1.0 / $1}')

    CMD=""
    for sub in ${REG_LABELS_2[*]}; do
      CMD="$CMD -m MSQ[$INITTEMP_DIR/group${grp}/work/template_${grp}_${sub}.nii.gz,$DATADIR/${id}_${side}_${sub}.nii.gz,$WGT_2]"
    done

    if [[ $ANTs_x_2 == "Y" ]]; then
       ANTs_mask_2="-x $INITTEMP_DIR/group${grp}/work/template_mask_${grp}.nii.gz"
    else
       ANTs_mask_2=""
    fi

    ANTS 3 $CMD \
      -t $ANTs_t_2 \
      -r $ANTs_r_2 \
      -i $ANTs_i_2 \
      $ANTs_mask_2 \
      -a $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_mlaffine_itk.txt \
      --continue-affine 0 \
      $ANTs_all_metrics_2 \
      -o $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp.nii.gz 

    for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        $DATADIR/${id}_${side}_${sub}.nii.gz \
        $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp_reslice_${sub}.nii.gz \
        -R $INITTEMP_DIR/group${grp}/work/template_${grp}_tse.nii.gz \
        $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempWarp.nii.gz \
        $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempAffine.txt

    done

  fi
}

function main_loop()
{
  # Main iteration loop
  PREFIX=ANTs${expid}
  for ((iter=0;iter<$ITER;iter++)); do
    for grp in ${groups[*]}; do

      # Create the segmentation for the template
      c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITTEMP_DIR/group${grp}/work/template_${grp}_${sub}.nii.gz; done) \
        -vote -type ushort \
        -o $INITTEMP_DIR/group${grp}/work/template_${grp}_seg.nii.gz

      # Back up template
      ITDIR=$INITTEMP_DIR/group${grp}/work/$(printf iter_%02d $iter)
      mkdir -p $ITDIR
      cp -a $INITTEMP_DIR/group${grp}/work/template_*.nii.gz $ITDIR/

      # Do ants?
      if [[ $iter -lt $ANTs_start ]]; then doants=0; else doants=1; fi

      GRP_IDS=($(cat $INITGRP_DIR/subjID_${grp}.txt))

      # Run ANTS for each image
      for grpid in ${GRP_IDS[*]}; do

        id=$(echo $grpid | cut -d '_' -f 1)
        side=$(echo $grpid | cut -d '_' -f 2)
        id=$(ls $ASHSRUNDIR | grep $id)

        #Submit ANTS job
        qsub -cwd -o $DUMPDIR -j y \
             -q all.q,basic.q \
             -N "${PREFIX}_${id}_${side}_${grp}" \
             $0 ants_iter $id $side $grp $doants

      done
    done
    
    # Wait for completion
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -hold_jid "${PREFIX}_*" -sync y -b y \
         sleep 1

    # If this is the last iteration, we don't want to recompute the template
    if [[ $iter -lt $((ITER-1)) ]]; then
      for grp in ${groups[*]}; do

        # Compute average images
        for kind in $KINDS; do

          if [[ $kind = "tse" ]]; then NORM=1; else NORM=0; fi
          AverageImages 3 $INITTEMP_DIR/group${grp}/work/template_${grp}_${kind}.nii.gz $NORM $INITTEMP_DIR/group${grp}/work/*_${grp}_totemp_reslice_${kind}.nii.gz

        done

        # Perform shape averaging
        if [[ $doants -eq 1 ]]; then
          shape_update_to_template $grp
        fi

      done

    fi
  done
}

##############################################################################
# Make images for further analysis
function make_images()
{
  for grp in ${groups[*]}; do

    # Generate one compound label for all other subfields
    ITDIR=$INITTEMP_DIR/group${grp}/work/$(printf iter_%02d $((ITER-1)) )
    for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

      LABELDEF=(${EXTRAMESHESDEF[i]})
      echo $LABELDEF

      c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$INITTEMP_DIR/group${grp}/work/template_${grp}_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
        -accum -add -endaccum \
        -o $ITDIR/template_${grp}_${EXTRAMESHES[i]}.nii.gz

    done

    # for each subject
    #IDS=$(ls ${expdir}/data | grep ${side}_${grp}_tse \
    #  | sed -e "s/_${side}_${grp}_tse.nii.gz//")
    GRP_IDS=($(cat $INITGRP_DIR/subjID_${grp}.txt))

    for grpid in ${GRP_IDS[*]}; do

      id=$(echo $grpid | cut -d '_' -f 1)
      side=$(echo $grpid | cut -d '_' -f 2)
      id=$(ls $ASHSRUNDIR | grep $id)

      # Create seg for each subject
      c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp_reslice_${sub}.nii.gz; done) -vote -type ushort \
        -o $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp_reslice_seg.nii.gz

    done
  done
}

##############################################################################
# Measurement dice overlap when warp template back to original subject space
# Warp labels back to subject space and vote
function warp_labels()
{
  PREFIX=WL${expid}
  for grp in ${groups[*]}; do

    mkdir -p $INITTEMP_DIR/group${grp}/labelwarp

    GRP_IDS=($(cat $INITGRP_DIR/subjID_${grp}.txt))

    for grpid in ${GRP_IDS[*]}; do

      id=$(echo $grpid | cut -d '_' -f 1)
      side=$(echo $grpid | cut -d '_' -f 2)
      #id=$(ls $ASHSRUNDIR | grep $id)

      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
        -q all.q,basic.q \
        -N "${PREFIX}_${id}_${side}_${grp}" \
        $0 warp_labels_subj $id $side $grp

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
    -q all.q,basic.q \
    -hold_jid "${PREFIX}_*" -sync y -b y \
    sleep 1
}

function warp_labels_subj()
{
  id=$1
  side=$2
  grp=$3
  SEG=$CLEANUPDIR/${id}_${side}_dividedCS.nii.gz
  id=$(ls $ASHSRUNDIR | grep $id)

  # add on transform
  if [[ $side == "left" ]]; then
    ADD_ON_TRANS=$IDENTITYFN
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
  fi

  # Apply the warp to each of the labels
  for sub in ${LABEL_IDS[*]}; do

    antsApplyTransforms -d 3 \
      -r $SEG \
      -i $INITTEMP_DIR/group${grp}/work/template_${grp}_${sub}.nii.gz \
      -o $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${sub}_tempfit.nii.gz \
      -t [$ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt,1] \
      -t [$ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt,1] \
      -t $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii \
      -t [$ADD_ON_TRANS,1] \
      -t [$INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempAffine.txt,1] \
      -t $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempInverseWarp.nii.gz \
      -n BSpline

  done

  # Vote in original subject space to get the tempfit segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${sub}_tempfit.nii.gz; done) -vote -type ushort \
    -o $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_tempfit.nii.gz

  # change ASHS labels to the current formate
  for ((i=0;i<${#LABEL_IDS[*]};i++)); do
    c3d $SEG \
      -replace $(for k in ${LABEL_MRG[i]}; do echo $k 999; done) \
      -thresh 999 999 1 0 \
      -o $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${LABEL_IDS[i]}_ashs.nii.gz
  done

  c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${sub}_ashs.nii.gz; done) -vote -type ushort \
    -o $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_ashs.nii.gz

  for ((i=0;i<${#LABEL_IDS[*]};i++)); do
    rm $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${LABEL_IDS[i]}_ashs.nii.gz
    rm $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${LABEL_IDS[i]}_tempfit.nii.gz
  done

  # change manual segmentation to the current formate
  if [[ -d $ASHSRUNDIR/$id/refseg ]]; then
    for ((i=0;i<${#LABEL_IDS[*]};i++)); do
      c3d $ASHSRUNDIR/${id}/refseg/refseg_${side}.nii.gz \
        -replace $(for k in ${LABEL_MRG[i]}; do echo $k 999; done) \
        -thresh 999 999 1 0 \
        -o $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${LABEL_IDS[i]}_refseg.nii.gz
    done

    c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${sub}_refseg.nii.gz; done) -vote -type ushort \
      -o $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_refseg.nii.gz

    for ((i=0;i<${#LABEL_IDS[*]};i++)); do
      rm $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${LABEL_IDS[i]}_refseg.nii.gz
    done
  fi
}

#########################################################
# Extract meshes for the subfields and apply warps to these meshes. This allows us to# perform statistical analysis on the mesh boundaries
function warp_meshes()
{
  PREFIX=WM${expid}
  for grp in ${groups[*]}; do

    mkdir -p $INITTEMP_DIR/group${grp}/meshwarp
    mkdir -p $INITTEMP_DIR/group${grp}/jacobian
    mkdir -p $INITTEMP_DIR/group${grp}/meshwarp/template

    # Generate meshes for the individual subfields
    for sub in ${LABEL_FG[*]}; do

      vtklevelset \
        $INITTEMP_DIR/group${grp}/work/template_${grp}_${sub}.nii.gz \
        $INITTEMP_DIR/group${grp}/meshwarp/template_${grp}_${sub}.vtk \
        $subfield_th

    done

    # Generate one mesh for all non-DG non-CS subfields
    ITDIR=$INITTEMP_DIR/group${grp}/work/$(printf iter_%02d $((ITER-1)) )
    for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

      vtklevelset \
        $ITDIR/template_${grp}_${EXTRAMESHES[i]}.nii.gz \
        $INITTEMP_DIR/group${grp}/meshwarp/template_${grp}_${EXTRAMESHES[i]}.vtk \
        $template_th

      done

    GRP_IDS=($(cat $INITGRP_DIR/subjID_${grp}.txt))

    for grpid in ${GRP_IDS[*]}; do

      id=$(echo $grpid | cut -d '_' -f 1)
      side=$(echo $grpid | cut -d '_' -f 2)
      id=$(ls $ASHSRUNDIR | grep $id)

      # Submit job for this subject
      qsub  -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -N "${PREFIX}_${id}_${side}_${grp}" \
          $0 warp_meshes_subj $id $side $grp

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  for grp in ${groups[*]}; do
    cp $INITTEMP_DIR/group${grp}/meshwarp/template*.vtk \
       $INITTEMP_DIR/group${grp}/meshwarp/template/
  done
}

function warp_meshes_subj()
{
  id=$1
  side=$2
  grp=$3
  ALLSF="${LABEL_FG[*]} ${EXTRAMESHES[*]}"

  # add on transform
  if [[ $side == "left" ]]; then
    ADD_ON_TRANS=$IDENTITYFN
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
  fi

  # Compose the transformation between the template and the subject
  $ANTSPATH/ComposeMultiTransform 3 \
     $TMPDIR/compose.nii \
     -R $INITTEMP_DIR/group${grp}/work/template_${grp}_tse.nii.gz \
     $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempWarp.nii.gz \
     $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempAffine.txt \
     $ADD_ON_TRANS \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempWarp.nii \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
     $ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

  # Split the transformations into xyz
  c3d -mcs $TMPDIR/compose.nii -oo $TMPDIR/comp%02d.nii

  # Apply the warp to each of the meshes
  for sub in $ALLSF; do

    # Warp the subfield into subject space
    warpmesh -w ants -m ras \
      $INITTEMP_DIR/group${grp}/meshwarp/template_${grp}_${sub}.vtk \
      $INITTEMP_DIR/group${grp}/meshwarp/${id}_${side}_${grp}_${sub}_tempfit.vtk \
      $TMPDIR/comp00.nii $TMPDIR/comp01.nii $TMPDIR/comp02.nii

    # Extract the thickness of the subfield
    cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
      -T $INITTEMP_DIR/group${grp}/meshwarp/${id}_${side}_${grp}_${sub}_thickmap.vtk \
      -p $thick_p -e $thick_e \
      $INITTEMP_DIR/group${grp}/meshwarp/${id}_${side}_${grp}_${sub}_tempfit.vtk \
      $TMPDIR/group${grp}/meshwarp/skel_${id}_${side}_${grp}_${sub}.vtk

  done

  # Compute the Jacobian map
  #ANTSJacobian 3 $TMPDIR/compose.nii \
  #  $INITTEMP_DIR/group${grp}/jacobian/${id}_${side}_${grp}_totemp 1
}

##############################################################################
# Evaluate the performance
# compute dice of mask
# compute dice of none DG none CS overlap
# compute mesh difference

function evaluation()
{
  INITEVALTMPDIR=$INITEVALDIR/tmp
  rm -rf $INITEVALTMPDIR
  mkdir -p $INITEVALTMPDIR

  PREFIX=EV${expid}
  for grp in ${groups[*]}; do

    #IDS=$(ls ${expdir}/data | grep ${side}_${grp}_tse | sed -e "s/_${side}_${grp}_tse.nii.gz//")
    GRP_IDS=($(cat $INITGRP_DIR/subjID_${grp}.txt))

    for grpid in ${GRP_IDS[*]}; do

      id=$(echo $grpid | cut -d '_' -f 1)
      side=$(echo $grpid | cut -d '_' -f 2)
      id=$(ls $ASHSRUNDIR | grep $id)  

      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
        -q all.q,basic.q \
        -N "${PREFIX}_${id}_${side}_${grp}" \
        $0 eval_subj $id $side $grp $INITEVALTMPDIR
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # Combine scores
  IDs=($(cat $SUBJ_TXT))
  echo "ID Side ${EVALLABELS[*]} ${MESH_EVAL[*]} ${EVALLABELS[*]}" \
    > $INITEVALDIR/overlap.txt

  echo "ID Side ${EVALLABELS[*]} ${MESH_EVAL[*]}" \
    > $INITEVALDIR/overlap_refseg.txt

  for side in left right; do

    echo "ID ${EVALLABELS[*]} ${MESH_EVAL[*]} ${EVALLABELS[*]}" \
      > $INITEVALDIR/overlap_${side}.txt

    echo "ID ${EVALLABELS[*]} ${MESH_EVAL[*]}" \
      > $INITEVALDIR/overlap_refseg_${side}.txt
  
  done

  for side in left right; do
    for ids in ${IDs[*]}; do

      id=$(ls $ASHSRUNDIR | grep $ids)
      cat $INITEVALTMPDIR/${id}_${side}_overlap.txt \
          >> $INITEVALDIR/overlap_${side}.txt

      cat $INITEVALTMPDIR/${id}_${side}_overlap_combined.txt \
          >> $INITEVALDIR/overlap.txt

      if [[ -f $INITEVALTMPDIR/${id}_${side}_refseg_overlap.txt ]]; then
        cat $INITEVALTMPDIR/${id}_${side}_refseg_overlap.txt \
          >> $INITEVALDIR/overlap_refseg_${side}.txt

        cat $INITEVALTMPDIR/${id}_${side}_refseg_overlap_combined.txt \
          >> $INITEVALDIR/overlap_refseg.txt

      fi

    done
  done
}

function eval_subj()
{
  id=$1
  side=$2
  grp=$3
  EVALTMPDIR=$4
  ALLOVL=""

  ###################################
  # Compute dice overlap in subject space
  FITTED=$INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_tempfit.nii.gz
  ASHSSEG=$INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_ashs.nii.gz
  do_pair $ASHSSEG $FITTED
  ALLOVL="$ALLOVL $FULLOVL"

  ###################################
  # Compute H distance
  MESH_DIST=""
  for ((i=0;i<${#MESH_EVAL[*]};i++)); do

    c3d $ASHSSEG \
      -replace $(for k in ${MESH_EVAL_RANGES[i]}; do echo $k 999; done) \
      -thresh 999 999 1 0 \
      -trim 2mm -resample 100x100x300% \
      -o $TMPDIR/mask_seg.nii.gz

    vtklevelset $TMPDIR/mask_seg.nii.gz $TMPDIR/mask_seg.vtk 0.5
    SIZE=$(cat $TMPDIR/mask_seg.vtk | wc -l)
    if [[ $SIZE -gt 6 ]]; then

      MESH_DIST_TMP="$(meshdiff $INITTEMP_DIR/group${grp}/meshwarp/${id}_${side}_${grp}_${MESH_EVAL[i]}_tempfit.vtk $TMPDIR/mask_seg.vtk | grep RESULT | awk '{print $7}')"

      if [[ $MESH_DIST_TMP == "" ]]; then
        MESH_DIST="$MESH_DIST NA"
      else
        MESH_DIST="$MESH_DIST $MESH_DIST_TMP"
      fi

    else

      MESH_DIST="$MESH_DIST NA"

    fi

  done

  ALLOVL="$ALLOVL $MESH_DIST"

  ###################################
  # Compute dice overlap in template space
  do_pair $INITTEMP_DIR/group${grp}/work/template_${grp}_seg.nii.gz \
    $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp_reslice_seg.nii.gz
  ALLOVL="$ALLOVL $FULLOVL"

  ###################################
  # save to temporary file
  echo $id $ALLOVL > $EVALTMPDIR/${id}_${side}_overlap.txt
  echo $id $side $ALLOVL > $EVALTMPDIR/${id}_${side}_overlap_combined.txt

  ###################################
  ALLOVL=""
  # Compute dice overlap between warp seg and manual seg
  FITTED=$INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_tempfit_ref.nii.gz
  REFSEG=$INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_refseg.nii.gz
  if [[ -f $REFSEG ]]; then
  # compute dice
  c3d $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_tempfit.nii.gz -replace 5 6 -o $FITTED
  do_pair $REFSEG $FITTED
  ALLOVL="$ALLOVL $FULLOVL"

  # Compute H distance
  MESH_DIST=""
  for ((i=0;i<${#MESH_EVAL[*]};i++)); do

    c3d $REFSEG \
      -replace $(for k in ${MESH_EVAL_RANGES[i]}; do echo $k 999; done) \
      -thresh 999 999 1 0 \
      -trim 2mm -resample 100x100x300% \
      -o $TMPDIR/mask_seg.nii.gz

    vtklevelset $TMPDIR/mask_seg.nii.gz $TMPDIR/mask_seg.vtk 0.5
    SIZE=$(cat $TMPDIR/mask_seg.vtk | wc -l)
    if [[ $SIZE -gt 6 ]]; then

      MESH_DIST_TMP="$(meshdiff $INITTEMP_DIR/group${grp}/meshwarp/${id}_${side}_${grp}_${MESH_EVAL[i]}_tempfit.vtk $TMPDIR/mask_seg.vtk | grep RESULT | awk '{print $7}')"

      if [[ $MESH_DIST_TMP == "" ]]; then
        MESH_DIST="$MESH_DIST NA"
      else
        MESH_DIST="$MESH_DIST $MESH_DIST_TMP"
      fi

    else

      MESH_DIST="$MESH_DIST NA"

    fi

  done

  ALLOVL="$ALLOVL $MESH_DIST"

  # save to temporary file
  echo $id $ALLOVL > $EVALTMPDIR/${id}_${side}_refseg_overlap.txt
  echo $id $side $ALLOVL > $EVALTMPDIR/${id}_${side}_refseg_overlap_combined.txt
  fi
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
    FULLOVL="$FULLOVL $OVL"

  done
}

##############################################################################
# Display statistics
ALLSF="${LABEL_FG[*]} ${EXTRAMESHES[*]}"
function disp_stats()
{
  PREFIX=STATS${expid}
  for grp in ${groups[*]}; do

    rm -rf $INITTEMP_DIR/group${grp}/meshwarp/analysis
    mkdir -p $INITTEMP_DIR/group${grp}/meshwarp/analysis

    for side in left right; do
      for sub in $ALLSF; do

        # Submit job for this subject
        qsub -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -N "${PREFIX}_${side}_${grp}_${sub}" \
          $0 disp_stats_sub $side $grp $sub

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function disp_stats_sub()
{
  side=$1
  grp=$2
  sub=$3

  #IDS=$(ls ${expdir}/data | grep ${side}_${grp}_tse | sed -e "s/_${side}_${grp}_tse.nii.gz//")
  IDS=""
  GRP_IDS=$(cat $INITGRP_DIR/subjID_${grp}.txt)
  for grpid in ${GRP_IDS[*]}; do
      id=$(echo $grpid | cut -d '_' -f 1)
      grpside=$(echo $grpid | cut -d '_' -f 2)
      if [[ $side == $grpside ]]; then
        IDS="$IDS $id"
      fi
  done

  # Displacement analysis
  MESHES=$(for id in $IDS; do \
    echo $(find $INITTEMP_DIR/group${grp}/meshwarp | grep $id | grep $side | grep ${sub}_tempfit ); done)

  meshdisp $MESHES \
    $INITTEMP_DIR/group${grp}/meshwarp/analysis/disp_${side}_${sub}.vtk

  # Thickness analysis
  MESHES=$(for id in $IDS; do \
    echo $(find $INITTEMP_DIR/group${grp}/meshwarp | grep $id | grep $side | grep ${sub}_thickmap ); done)

  mesh_merge_arrays \
    -r $INITTEMP_DIR/group${grp}/meshwarp/template_${grp}_${sub}.vtk \
    $INITTEMP_DIR/group${grp}/meshwarp/analysis/thick_${side}_${sub}.vtk \
    Thickness $MESHES
}

##############################################################################
# For thickness analysis, we can run analysis on multiple groups of meshes at once
function thick_stats()
{
  PREFIX=THICK${expid}
  for grp in ${groups[*]}; do

    rm -rf $INITTEMP_DIR/group${grp}/meshwarp/analysis/design*

    for side in left right; do

     for design in $(ls $GPSTATDIR | grep "design_.*txt" | grep $side | grep $grp ); do

        exp=$(echo $design | sed -e "s/^design_${side}_${grp}_//" | sed -e "s/\.txt//")

        for ((igrp=0;igrp<${#ANGRP[*]};igrp++)); do

          # Submit job for this subject
          qsub -cwd -o $DUMPDIR -j y \
            -q all.q,basic.q \
            -N "${PREFIX}_${exp}_${side}_${grp}_${GRPNM[igrp]}" \
            $0 thick_stats_sub $exp $side $grp $igrp
          #thick_stats_sub $exp $side $grp $igrp

        done
      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function thick_stats_sub()
{
  exp=$1
  side=$2
  grp=$3
  igrp=$4

  MYGRP=${ANGRP[igrp]}
  GNAME=${GRPNM[igrp]}

  # Create the work directory for this analysis
  WORK=$INITTEMP_DIR/group${grp}/meshwarp/analysis/design_${exp}_group_${GNAME}
  mkdir -p $WORK

  # Merge the meshes for this analysis
  for sub in $MYGRP; do

    # Get the list of subjects
    DESIGNTXT=$GPSTATDIR/design_${side}_${grp}_${exp}.txt
    SUBJ=$(cat $DESIGNTXT | awk '{print $1}')

    # Generate the design matrix for meshglm
    cat $DESIGNTXT | awk '{$1=""; print}' > $WORK/design_${side}.txt

    MESHES=$(for id in $SUBJ; do \
      echo $(find $INITTEMP_DIR/group${grp}/meshwarp | grep ${id} | grep ${side}_${grp}_${sub}_thickmap.vtk); done)

    mesh_merge_arrays -r \
      $INITTEMP_DIR/group${grp}/meshwarp/template_${grp}_${sub}.vtk \
      $WORK/thick_${side}_${grp}_${sub}.vtk Thickness $MESHES

  done

  # Go through the list of contrasts
  for con in $(ls $GPSTATDIR | grep "contrast_${exp}_.*\.txt"); do

    # Get the suffix for this contrast
    suffix=$(echo $con | sed -e "s/^contrast_${exp}_//" -e "s/\.txt//")

    # Create the directory for this contrast
    CWORK=$WORK/contrast_${suffix}
    mkdir -p $CWORK

    FULLNM="design_${exp}_group_${GNAME}_con_${suffix}"

    # Copy the contrast
    cp $GPSTATDIR/$con $CWORK/contrast.txt

    # Build the list of meshes to include
    MESHPARAM=""
    for sub in $MYGRP; do
      MESHPARAM="$MESHPARAM -m $WORK/thick_${side}_${grp}_${sub}.vtk $CWORK/thickstat_${FULLNM}_${side}_${grp}_${sub}.vtk"

    done

    meshglm $MESHPARAM \
      -g $WORK/design_${side}.txt $CWORK/contrast.txt \
      -a Thickness \
      -d 4 \
      -p 1000 \
      -s T \
      -t 2.4 \
      -e

  done
}

##############################################################################
function mean_thickness()
{
  PREFIX=MT${expid}
  for grp in ${groups[*]}; do
    MESHGRPDIR=$INITTEMP_DIR/group${grp}/meshwarp
    rm -rf $MESHGRPDIR/analysis/labels
    mkdir -p $MESHGRPDIR/analysis/labels
    for side in left right; do

      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=2.1G,s_vmem=2G \
           -N "${PREFIX}_${grp}_${side}" \
           $0 mean_thickness_sub $grp $side
      sleep 0.1
      #mean_thickness_sub $grp $side

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # generate csv file with summary measurements!
  IDS=($(cat $SUBJ_TXT))
  rm -rf $INITTEMP_DIR/analysis
  mkdir -p $INITTEMP_DIR/analysis

  # generate header
  header="ID"
  for side in left right; do
    for ((i=0;i<${#MESH_LABEL[*]};i++)); do
      sub=${MESH_LABEL[$i]}
      header="${header},${sub}_${side}"
    done
    header="${header},aBA35_${side},pBA35_${side}"
  done
  echo "$header" > $INITTEMP_DIR/analysis/summary_thick.csv

  # combine files
  for ((i=0;i<${#IDS[*]};i++)); do

    # id
    id=${IDS[$i]}

    # load left side
    side=left
    GRP=($(cat $INITGRP_DIR/group_${side}.txt))
    grp=${GRP[$i]} 
    THICK_LEFT=$(cat $INITTEMP_DIR/group${grp}/meshwarp/analysis/labels/thick_${side}_${grp}_MRG.txt | grep $id)

    # load right side
    side=right
    GRP=($(cat $INITGRP_DIR/group_${side}.txt))
    grp=${GRP[$i]}
    THICK_RIGHT=$(cat $INITTEMP_DIR/group${grp}/meshwarp/analysis/labels/thick_${side}_${grp}_MRG.txt | grep $id | cut -d ' ' -f 2)

    echo "${THICK_LEFT}${THICK_RIGHT}" >> \
      $INITTEMP_DIR/analysis/summary_thick.csv

  done

  # perform statistical analysis
  Rscript $WORKDIR/scripts/summary_thickness_analysis.R \
    $WORKDIR/analysis_input/demog.csv \
    $INITTEMP_DIR/analysis/summary_thick.csv \
    $INITTEMP_DIR/analysis/mean_thickness_stats.csv
}

function mean_thickness_sub()
{
  grp=$1
  side=$2

  # sample on probability maps (don`t run BKG)
  MESHGRPDIR=$INITTEMP_DIR/group${grp}/meshwarp
  WORKGRPDIR=$INITTEMP_DIR/group${grp}/work
  for ((ii=0;ii<${#EXTRAMESHES[*]};ii++)); do
    LABELDEF=(${EXTRAMESHESDEF[ii]})
    MESHES=""
    for ((i=0;i<${#LABEL_IDS[*]};i++)); do
      if [[ ${LABELDEF[$i]} == "1" ]]; then
        sub=${LABEL_IDS[$i]}
        mesh_image_sample \
          $MESHGRPDIR/template_${grp}_${EXTRAMESHES[ii]}.vtk \
          $WORKGRPDIR/template_${grp}_${sub}.nii.gz \
          $MESHGRPDIR/analysis/labels/template_${side}_${grp}_${EXTRAMESHES[ii]}_${sub}.vtk \
          PROB
        MESHES="$MESHES $MESHGRPDIR/analysis/labels/template_${side}_${grp}_${EXTRAMESHES[ii]}_${sub}.vtk"
      fi
    done

    # merge prob arrays
    mesh_merge_arrays \
      -r $MESHGRPDIR/analysis/thick_${side}_${EXTRAMESHES[ii]}.vtk \
      $MESHGRPDIR/analysis/labels/thick_${side}_${grp}_${EXTRAMESHES[ii]}.vtk \
      PROB $MESHES

    # run matlab script to generate label
      $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
        addpath('/data/picsl/longxie/ASHS/thickness_newlabel/matlabcode');
        compute_mesh_label('$MESHGRPDIR/analysis/labels/thick_${side}_${grp}_${EXTRAMESHES[ii]}.vtk');
MATCODE

  done

  # prepare the subjID_${side}_${grp}.txt
  GRP_IDS=$(cat $INITGRP_DIR/subjID_${grp}.txt)
  rm -f $MESHGRPDIR/analysis/labels/subjID_${side}_${grp}.txt
  for grpid in ${GRP_IDS[*]}; do
      id=$(echo $grpid | cut -d '_' -f 1)
      grpside=$(echo $grpid | cut -d '_' -f 2)
      if [[ $side == $grpside ]]; then
        echo $id >> $MESHGRPDIR/analysis/labels/subjID_${side}_${grp}.txt
      fi
  done

  # run matlab script to generate label
  mkdir -p $TMPDIR/${id}
  source $PKGDIR/matlab_batch.sh \
    $TMPDIR/${id}/ \
    $MATLABCODEDIR \
    compute_mean_thickness_apBA35 \
    $MESHGRPDIR/analysis/labels/thick_${side}_${grp}_MRG.vtk \
    2 \
    $MESHGRPDIR/analysis/labels/subjID_${side}_${grp}.txt
}

##############################################################################
function Insertion()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=IS${expid}
  for grp in ${groups[*]}; do      

    # iter directory
    ITDIR=$INITTEMP_DIR/group${grp}/work/$(printf iter_%02d $((ITER-1)) )

    # link templates to data folder
    for kind in $KINDS seg; do
        
      ln -sf $ITDIR/template_${grp}_${kind}.nii.gz \
        $DATADIR/template${grp}_${kind}.nii.gz

    done

    for id in $IDS; do

      fn_id=$(ls $ASHSRUNDIR | grep $id)

      for side in left right; do
        # perform registration from templates to subjects
        fnside=${fn_id}_${side}
        qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -N "${PREFIX}_${fnside}_${grp}" $0 \
           Insertion_sub template${grp} $fnside
        sleep 0.1

        # perform registration from subjects to templates
        fnside=${fn_id}_${side}
        qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -N "${PREFIX}_${grp}_${fnside}" $0 \
           Insertion_sub $fnside template${grp}
        sleep 0.1
      done

    done
  done

  # perform registration between templates
  for grp_src in ${groups[*]}; do
    for grp_tg in ${groups[*]}; do

      if [[ $grp_src -eq $grp_tg ]]; then
        continue
      fi
        
      qsub -cwd -o $DUMPDIR -j y \
        -q all.q,basic.q \
        -N "${PREFIX}_${grp_src}_${grp_tg}" $0 \
        Insertion_sub \
        template${grp_src} template${grp_tg}

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function Insertion_sub()
{
  id_fix=$1
  id_mov=$2
  OUTDIR=$PWDIR/${id_fix}/${id_mov}_to_${id_fix}

  if [[ $id_mov != $id_fix ]]; then

    if [[ -d $OUTDIR ]]; then

      echo "Seg file exists."
      rm -rf $OUTDIR

    fi

    mkdir -p $OUTDIR

    # Use ml_affine for nice affine alignment
    /data/picsl/pauly/wolk/ashs/ext/Linux/bin/ml_affine \
      $DATADIR/${id_fix}_seg.nii.gz \
      $DATADIR/${id_mov}_seg.nii.gz \
      $TMPDIR/${id_mov}_to_${id_fix}_mlaffine.txt

    # Convert that to ITK format
    c3d_affine_tool \
      $TMPDIR/${id_mov}_to_${id_fix}_mlaffine.txt \
      -oitk $TMPDIR/${id_mov}_to_${id_fix}_mlaffine_itk.txt

    CMD=""
    for sub in ${REG_LABELS_1[*]}; do
      CMD="$CMD -m MSQ[$DATADIR/${id_fix}_${sub}.nii.gz,$DATADIR/${id_mov}_${sub}.nii.gz,$WGT_1]"
    done

    # Perform ANTs registration
    ANTS 3 $CMD \
       -t $ANTs_t_1 \
       -r $ANTs_r_1 \
       -i $ANTs_i_1 \
       -a $TMPDIR/${id_mov}_to_${id_fix}_mlaffine_itk.txt \
       --continue-affine 0 \
       $ANTs_all_metrics_1 \
       -o $TMPDIR/${id_mov}_to_${id_fix}_pairwise.nii.gz 

    for sub in ${KINDS[*]}; do

      WarpImageMultiTransform 3 \
        $DATADIR/${id_mov}_${sub}.nii.gz \
        $TMPDIR/${id_mov}_to_${id_fix}_reslice_${sub}.nii.gz \
        -R $DATADIR/${id_fix}_tse.nii.gz \
        $TMPDIR/${id_mov}_to_${id_fix}_pairwiseWarp.nii.gz \
        $TMPDIR/${id_mov}_to_${id_fix}_pairwiseAffine.txt

    done

    # Create seg
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $TMPDIR/${id_mov}_to_${id_fix}_reslice_${sub}.nii.gz; done) \
        -vote -type ushort \
        -o $TMPDIR/${id_mov}_to_${id_fix}_reslice_seg.nii.gz

    # measure similarity
    OVL=$(c3d $TMPDIR/${id_mov}_to_${id_fix}_reslice_seg.nii.gz \
      -replace 1 0 4 0 6 0 \
      $DATADIR/${id_fix}_seg.nii.gz \
      -replace 1 0 4 0 6 0 \
      -label-overlap \
      | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}' )
    echo "$OVL" > \
      $OUTDIR/${id_mov}_to_${id_fix}_sim.txt


  fi
}

##############################################################################
function TempSim()
{
  mkdir -p $SIM_DIR
  IDS=$(cat $SUBJ_WITHTEMP_TXT)

  # Concadenate the results
  rm -rf ${SIM_DIR}/sim_${SIM_TYPE}_withtemp.txt

  for side_fix in left right; do
  for id_fix in $IDS; do

    if [[ $id_fix == "template1" ]] || \
       [[ $id_fix == "template2" ]] || \
       [[ $id_fix == "template3" ]]; then
      if [[ $side_fix == "left" ]]; then
        continue
      fi
      fn_fix=$id_fix
      idside_fix=${id_fix}
      fnside_fix=${fn_fix}
    else
      fn_fix=$(ls $ASHSRUNDIR | grep $id_fix)
      idside_fix=${id_fix}_${side_fix}
      fnside_fix=${fn_fix}_${side_fix}
    fi

    OVL=""
    for side_mov in left right; do
      for id_mov in $IDS; do

      if [[ $id_mov == "template1" ]] || \
         [[ $id_mov == "template2" ]] || \
         [[ $id_mov == "template3" ]]; then
        if [[ $side_mov == "left" ]]; then
          continue
        fi
        fn_mov=$id_mov
        idside_mov=${id_mov}
        fnside_mov=${fn_mov}
      else
        fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)
        idside_mov=${id_mov}_${side_mov}
        fnside_mov=${fn_mov}_${side_mov}
      fi

      if [[ $fnside_fix == $fnside_mov ]]; then

        OVL="$OVL 1"

      else

        SIMFILE=$PWDIR/${fnside_fix}/${fnside_mov}_to_${fnside_fix}/${fnside_mov}_to_${fnside_fix}_sim.txt
        OVL="$OVL $(cat $SIMFILE)"

      fi

      done
    done

    echo "$OVL" >> ${SIM_DIR}/sim_${SIM_TYPE}_withtemp.txt

  done
  done
}

############################################################
function FindPath()
{
  # Submit job to copy data
  PREFIX=FP${expid}
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -l h_vmem=4.1G,s_vmem=4G \
       -N "${PREFIX}_NOINPUT" \
       $0 FindPath_sub "NOINPUT"

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function FindPath_sub()
{
  INPUT=$1

  # run matlab script to generate label
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('/data/picsl/longxie/ASHS/thickness_newlabel/matlabcode');
    find_paths_in_PWsimMatrix_combined('$INITGRP_DIR/group.txt','${SIM_DIR}/sim_${SIM_TYPE}_withtemp.txt','$BESTPATHSDIR');
MATCODE
}

############################################################
function PropagateTemplatesNew()
{
  PREFIX=PTN${expid}
  for grp_fix in ${groups[*]}; do
  #for grp_fix in 1; do
    for grp_mov in ${groups[*]}; do
    #for grp_mov in 3; do

    if [[ $grp_fix == $grp_mov ]]; then
      continue
    fi

    # submit jobs
    qsub -cwd -o $DUMPDIR -j y \
      -q all.q,basic.q \
      -l h_vmem=4.1G,s_vmem=4G \
      -N "${PREFIX}_${grp_fix}_${grp_mov}" $0 \
      PropagateTemplatesNew_sub $grp_fix $grp_mov
    sleep 0.1

    done
  done

  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function PropagateTemplatesNew_sub()
{
  grp_fix=$1
  grp_mov=$2
  IDS_ALL=($(cat $SUBJ_TXT))

  # get the path
  paths=($(cat $BESTPATHSDIR/template${grp_mov}.txt))
  N=${#IDS_ALL[*]}
  N_grps=${#groups[*]}
  N_all=${#paths[*]}
  idx=$((grp_fix-1))
  cur_path=${paths[$idx]}

  # generate IDS array with templates
  IDS=""
  SIDES=""
  for ((i=0;i<$N;i++)); do
    id=${IDS_ALL[i]}
    IDS[$i]=$(ls $ASHSRUNDIR | grep $id)
    SIDES[$i]="left"
  done
  for ((i=0;i<$N;i++)); do
    id=${IDS_ALL[i]}
    idx=$((i+N))
    IDS[$idx]=$(ls $ASHSRUNDIR | grep $id)
    SIDES[$idx]="right"
  done
  for ((i=0;i<$N_grps;i++)); do
    k=$((i+1))
    idx=$((2*N))
    idx=$((idx+i))
    IDS[$idx]="template${k}"
    SIDES[$idx]="NULL"
  done

  # loop
  i=1
  # initial template to target template
  TRANSDIR=$FINALTEMPDIR/template${grp_fix}/tempwarps
  TRANSITERDIR=$TRANSDIR/${grp_mov}_to_${grp_fix}/${i}_template${grp_mov}/
  mkdir -p $TRANSITERDIR
  for sub in ${LABEL_IDS[*]}; do
    ln -sf $DATADIR/template${grp_mov}_${sub}.nii.gz \
      $TRANSITERDIR/template${grp_mov}_reslice_${sub}.nii.gz
  done
  ln -sf $DATADIR/template${grp_mov}_seg.nii.gz \
    $TRANSITERDIR/template${grp_mov}_reslice_seg.nii.gz

  # reslice seg
  c3d  $DATADIR/template${grp_fix}_tse.nii.gz \
    $TRANSITERDIR/template${grp_mov}_reslice_seg.nii.gz \
    -interp 0 -reslice-identity \
    -o $TRANSITERDIR/${i}_template${grp_mov}_reslice_to_template${grp_fix}_seg.nii.gz
  for ((ii=0;ii<${#EXTRAMESHES[*]};ii++)); do
    ln -sf $INITTEMP_DIR/group${grp_mov}/meshwarp/template_${grp_mov}_${EXTRAMESHES[ii]}.vtk $TRANSITERDIR/${i}_template${grp_mov}_${EXTRAMESHES[ii]}.vtk
  done

  # go through the path
  warps=""
  invwarps=""
  mkdir -p $TRANSDIR
  while true; do

    # from subject
    idx_from=$(echo $cur_path | cut -d , -f $i)
    id_from=${IDS[$((idx_from-1))]}
    if [[ $id_from == "template1" ]] || \
       [[ $id_from == "template2" ]] || \
       [[ $id_from == "template3" ]]; then
      idside_from=$id_from
    else
      side_from=${SIDES[$((idx_from-1))]}
      idside_from=${id_from}_${side_from}
    fi

    # to subject
    set +e
    i=$((i+1))
    idx_to=$(echo $cur_path | cut -d , -f $i)
    id_to=${IDS[$((idx_to-1))]}
    if [[ $id_to == "template1" ]] || \
       [[ $id_to == "template2" ]] || \
       [[ $id_to == "template3" ]]; then
      idside_to=$id_to
    else
      side_to=${SIDES[$((idx_to-1))]}
      idside_to=${id_to}_${side_to}
    fi
    set -e

    if [[ $idx_to == "" ]]; then
      break
    fi

    ####################################
    PRETRANSITERDIR=$TRANSDIR/${grp_mov}_to_${grp_fix}/$((i-1))_${idside_from}
    TRANSITERDIR=$TRANSDIR/${grp_mov}_to_${grp_fix}/${i}_${idside_to}
    #TRANSITERTMPDIR=$TRANSITERDIR/tmp
    mkdir -p $TRANSITERDIR

    # link the original subject seg
    for sub in ${LABEL_IDS[*]}; do
    ln -sf $DATADIR/${idside_to}_${sub}.nii.gz \
      $TRANSITERDIR/${idside_to}_${sub}.nii.gz
    done
    ln -sf $DATADIR/${idside_to}_seg.nii.gz \
      $TRANSITERDIR/${idside_to}_seg.nii.gz

    # perform registration between steps
    #if [ ! -f $TRANSITERDIR/${idside_from}_to_${idside_to}_pairwiseWarp.nii.gz ] || \
    #   [ ! -f $TRANSITERDIR/${idside_from}_to_${idside_to}_pairwiseInverseWarp.nii.gz ]; then

    # Use ml_affine for nice affine alignment
    /data/picsl/pauly/wolk/ashs/ext/Linux/bin/ml_affine \
      $TRANSITERDIR/${idside_to}_seg.nii.gz \
      $PRETRANSITERDIR/template${grp_mov}_reslice_seg.nii.gz \
      $TRANSITERDIR/${idside_from}_to_${idside_to}_mlaffine.txt

    # Convert that to ITK format
    c3d_affine_tool \
      $TRANSITERDIR/${idside_from}_to_${idside_to}_mlaffine.txt \
      -oitk \
      $TRANSITERDIR/${idside_from}_to_${idside_to}_mlaffine_itk.txt

    CMD=""
    for sub in ${REG_LABELS_4[*]}; do
      CMD="$CMD -m MSQ[$TRANSITERDIR/${idside_to}_${sub}.nii.gz,$PRETRANSITERDIR/template${grp_mov}_reslice_${sub}.nii.gz,$WGT_4]"
    done

    if [[ $ANTs_x_4 == "Y" ]]; then
      c3d $TRANSITERDIR/${idside_to}_seg.nii.gz -dup \
          $PRETRANSITERDIR/template${grp_mov}_reslice_seg.nii.gz \
          -int 0 -reslice-identity \
          -add -binarize -dilate 1 10x10x10vox \
          -o $TRANSITERDIR/${idside_from}_to_${idside_to}_mask.nii.gz
      ANTs_mask_1="-x $TRANSITERDIR/${idside_from}_to_${idside_to}_mask.nii.gz"
    else
      ANTs_mask_1=""
    fi

    # Perform ANTs registration
    ANTS 3 $CMD \
         -t $ANTs_t_4 \
         -r $ANTs_r_4 \
         -i $ANTs_i_4 \
         $ANTs_G_4 \
         $ANTs_SMOOTH_4 \
         $ANTs_mask_4 \
         -a $TRANSITERDIR/${idside_from}_to_${idside_to}_mlaffine_itk.txt \
         --continue-affine 0 \
         $ANTs_all_metrics_4 \
         -o $TRANSITERDIR/${idside_from}_to_${idside_to}_pairwise.nii.gz

    #fi

    # compose warps
    warps="$TRANSITERDIR/${idside_from}_to_${idside_to}_pairwiseWarp.nii.gz $TRANSITERDIR/${idside_from}_to_${idside_to}_mlaffine_itk.txt $warps"

    invwarps="$invwarps -i $TRANSITERDIR/${idside_from}_to_${idside_to}_mlaffine_itk.txt $TRANSITERDIR/${idside_from}_to_${idside_to}_pairwiseInverseWarp.nii.gz"

    $ANTSPATH/ComposeMultiTransform 3 \
      $TRANSITERDIR/template${grp_mov}_InitWarp.nii.gz \
      -R $TRANSITERDIR/${idside_to}_seg.nii.gz \
      $warps

    # apply init warps
    for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        $DATADIR/template${grp_mov}_${sub}.nii.gz \
        $TRANSITERDIR/template${grp_mov}_reslice_${sub}.nii.gz \
        -R $TRANSITERDIR/${idside_to}_seg.nii.gz \
        $TRANSITERDIR/template${grp_mov}_InitWarp.nii.gz

    done

    # Create the segmentation for the template
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $TRANSITERDIR/template${grp_mov}_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $TRANSITERDIR/template${grp_mov}_reslice_seg.nii.gz

    # resliced to output template
    c3d $DATADIR/template${grp_fix}_tse.nii.gz \
      $TRANSITERDIR/template${grp_mov}_reslice_seg.nii.gz \
      -interp 0 -reslice-identity \
      -o $TRANSITERDIR/${i}_template${grp_mov}_reslice_to_template${grp_fix}_seg.nii.gz

    # generate a surface mesh (TODO)
    $ANTSPATH/ComposeMultiTransform 3 \
      $TRANSITERDIR/template${grp_mov}_InitInverseWarp.nii \
      -R $DATADIR/template${grp_mov}_tse.nii.gz \
      $invwarps

    c3d -mcs $TRANSITERDIR/template${grp_mov}_InitInverseWarp.nii \
      -oo $TMPDIR/comp%02d.nii

    warpmesh -w ants -m ras \
      $TRANSDIR/${grp_mov}_to_${grp_fix}/1_template${grp_mov}/1_template${grp_mov}_MRG.vtk \
      $TRANSITERDIR/${i}_template${grp_mov}_MRG.vtk \
      $TMPDIR/comp00.nii $TMPDIR/comp01.nii $TMPDIR/comp02.nii

  done

  # compose transform
  $ANTSPATH/ComposeMultiTransform 3 \
    $TRANSDIR/template${grp_mov}_InitWarp.nii.gz \
    -R $DATADIR/template${grp_fix}_tse.nii.gz \
    $warps

  $ANTSPATH/ComposeMultiTransform 3 \
    $TRANSDIR/template${grp_mov}_InitInverseWarp.nii.gz \
    -R $DATADIR/template${grp_mov}_tse.nii.gz \
    $invwarps

  # apply init warps
  for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        $DATADIR/template${grp_mov}_${sub}.nii.gz \
        $TRANSDIR/template${grp_mov}_inittotemp_reslice_${sub}.nii.gz \
        -R $DATADIR/template${grp_fix}_tse.nii.gz \
        $TRANSDIR/template${grp_mov}_InitWarp.nii.gz

  done

  # Create the segmentation for the template
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $TRANSDIR/template${grp_mov}_inittotemp_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $TRANSDIR/template${grp_mov}_inittotemp_reslice_seg.nii.gz
}

##############################################################################
function WarpSuperTemplate()
{
  PREFIX=WST${expid}
  for grp_subj in ${groups[*]}; do

    GRP_IDS=($(cat $INITGRP_DIR/subjID_${grp_subj}.txt))

    for grpid in ${GRP_IDS[*]}; do

      id=$(echo $grpid | cut -d '_' -f 1)
      side=$(echo $grpid | cut -d '_' -f 2)
      fn_id=$(ls $ASHSRUNDIR | grep $id)

      # for all three templates
      for grp_temp in ${groups_temp[*]}; do

        # Submit job for this subject
        qsub -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -N "${PREFIX}_${side}_${grp_temp}_${id}_${grp_subj}" $0 \
          WarpSuperTemplate_sub $side $grp_temp $id $grp_subj

      done
    done
  done

  # wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function AverageFinalTemplate()
{
  # generate the average template
  for grp_temp in ${groups_temp[*]}; do

    FINALWORKDIR=$FINALTEMPDIR/template${grp_temp}/finalwork/
    mkdir -p $FINALWORKDIR

    # Compute average images
    for kind in $KINDS; do

      if [[ $kind = "tse" ]]; then NORM=1; else NORM=0; fi
      AverageImages 3 $FINALWORKDIR/template_${kind}.nii.gz $NORM $FINALTEMPDIR/template${grp_temp}/work/*_totemp_reslice_${kind}.nii.gz

    done

    # Create the segmentation for the template
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALWORKDIR/template_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $FINALWORKDIR/template_seg.nii.gz

    # generate mask
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALWORKDIR/template_${sub}.nii.gz; done | grep -v BKG) \
      -mean -thresh 1e-5 inf 1 0 -trim 5vox \
      -o $FINALWORKDIR/template_mask.nii.gz

  done
}

function WarpSuperTemplate_sub()
{
  side=$1
  grp_temp=$2
  id=$3
  grp_subj=$4
  fn_id=$(ls $ASHSRUNDIR | grep $id)

  FINALTEMPLATEWORKDIR=$FINALTEMPDIR/template${grp_temp}/work
  mkdir -p $FINALTEMPLATEWORKDIR

  for sub in $KINDS; do

    warps="$INITTEMP_DIR/group${grp_subj}/work/${fn_id}_${side}_${grp_subj}_totempWarp.nii.gz $INITTEMP_DIR/group${grp_subj}/work/${fn_id}_${side}_${grp_subj}_totempAffine.txt"

    invwarps="-i $INITTEMP_DIR/group${grp_subj}/work/${fn_id}_${side}_${grp_subj}_totempAffine.txt $INITTEMP_DIR/group${grp_subj}/work/${fn_id}_${side}_${grp_subj}_totempInverseWarp.nii.gz"

    if [[ $grp_subj -ne $grp_temp ]]; then
      warps="$FINALTEMPDIR/template${grp_temp}/tempwarps/template${grp_subj}_InitWarp.nii.gz $warps"

      invwarps="$invwarps $FINALTEMPDIR/template${grp_temp}/tempwarps/template${grp_subj}_InitInverseWarp.nii.gz"
    fi

    # compose warps
    $ANTSPATH/ComposeMultiTransform 3 \
      $FINALTEMPLATEWORKDIR/${fn_id}_${side}_totempWarp.nii.gz \
      -R $DATADIR/template${grp_temp}_tse.nii.gz \
      $warps  

    # compose inverse warps
    $ANTSPATH/ComposeMultiTransform 3 \
      $FINALTEMPLATEWORKDIR/${fn_id}_${side}_totempInverseWarp.nii.gz \
      -R $DATADIR/${fn_id}_${side}_tse.nii.gz \
      $invwarps

    WarpImageMultiTransform 3 \
      $DATADIR/${fn_id}_${side}_${sub}.nii.gz \
      $FINALTEMPLATEWORKDIR/${fn_id}_${side}_totemp_reslice_${sub}.nii.gz \
      -R $DATADIR/template${grp_temp}_tse.nii.gz \
      $FINALTEMPLATEWORKDIR/${fn_id}_${side}_totempWarp.nii.gz

    done
}

##############################################################################
function shape_update_to_finaltemplate()
{
  grp_temp=$1
  FINALWORKDIR=$FINALTEMPDIR/template${grp_temp}/finalwork/

  # Borrowed from ANTS buildtemplateparallel

  # Average the warp fields
  local TEMPWARP=$FINALWORKDIR/template_warp.nii.gz
  $NEWANTSDIR/AverageImages 3 $TEMPWARP 0 \
    $FINALWORKDIR/*_totempWarp.nii.gz

  # Scale by -0.25 (gradient step)
  MultiplyImages 3 $TEMPWARP -0.25 $TEMPWARP

  # Create a new template affine
  local TEMPAFF=$FINALWORKDIR/template_Affine.txt
  if [[ -f $TEMPAFF ]]; then rm -rf $TEMPAFF; fi

  cat $FINALWORKDIR/*_totempAffine.txt \
    | grep '^Parameters:' | awk '\
    BEGIN { for(i=0;i<12;i++) x[i]=0} \
    { for(i=0;i<12;i++) x[i]=x[i]+$(i+2) } \
    END { \
      printf "Transform: MatrixOffsetTransformBase_double_3_3\nParameters: "; \
      for(i=0;i<12;i++) printf "%f ",x[i]/NR; \
      printf "\nFixedParameters: 0 0 0\n";}' > $TEMPAFF

  # Compose the warps
  WarpImageMultiTransform 3 $TEMPWARP $TEMPWARP -i $TEMPAFF \
    -R $FINALWORKDIR/template_tse.nii.gz

  TEMPWARPFULL=$FINALWORKDIR/template_fullwarp.nii.gz
  $ANTSPATH/ComposeMultiTransform 3 \
    $TEMPWARPFULL -R $FINALWORKDIR/template_tse.nii.gz \
    -i $TEMPAFF $TEMPWARP $TEMPWARP $TEMPWARP $TEMPWARP

  # Apply this warp to all the template derivatives
  for kind in $KINDS; do

    WarpImageMultiTransform 3 \
      $FINALWORKDIR/template_${kind}.nii.gz \
      $FINALWORKDIR/template_${kind}.nii.gz \
      $TEMPWARPFULL

  done
}

function RegisterSuperTemplate()
{
  # choose the best template to be the super template
  PREFIX=RST${expid}
  for grp_temp in ${groups_temp[*]}; do
    FINALWORKDIR=$FINALTEMPDIR/template${grp_temp}/finalwork/
    for ((iter=0;iter<$FINALITER;iter++)); do

      # Back up template
      ITDIR=$FINALWORKDIR/$(printf iter_%02d $iter)
      mkdir -p $ITDIR
      cp -a $FINALWORKDIR/template_*.nii.gz $ITDIR/

      # Do ants?
      if [[ $iter -lt $ANTs_start_final ]]; then doants=0; else doants=1; fi

      IDS=$(cat $SUBJ_TXT)

      # Run ANTS for each image
      for id in $IDS; do

        id=$(ls $ASHSRUNDIR | grep $id)

        for side in left right; do

          #Submit ANTS job
          qsub -cwd -o $DUMPDIR -j y \
             -q all.q,basic.q \
             -N "${PREFIX}_${id}_${side}_${grp_temp}" $0 \
             RegisterSuperTemplate_sub $id $side $grp_temp $doants

        done
      done


      # Wait for completion
      qsub -cwd -o $DUMPDIR -j y \
        -q all.q,basic.q \
        -hold_jid "${PREFIX}_*" -sync y -b y \
        sleep 1

      # If this is the last iteration, we don't want to recompute the template
      if [[ $iter -lt $((FINALITER-1)) ]]; then
        # Compute average images
        for kind in $KINDS; do

          if [[ $kind = "tse" ]]; then NORM=1; else NORM=0; fi
          AverageImages 3 $FINALWORKDIR/template_${kind}.nii.gz $NORM $FINALWORKDIR/*_totemp_reslice_${kind}.nii.gz

        done

        # Create the segmentation for the template
        c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALWORKDIR/template_${sub}.nii.gz; done) \
          -vote -type ushort \
          -o $FINALWORKDIR/template_seg.nii.gz

        # Perform shape averaging
        if [[ $doants -eq 1 ]]; then
          shape_update_to_finaltemplate $grp_temp
        fi

      fi

    done
  done
}

function RegisterSuperTemplate_sub()
{
  id=$1
  side=$2
  grp_temp=$3
  doants=$4

  SRCWORKDIR=$FINALTEMPDIR/template${grp_temp}/work
  FINALWORKDIR=$FINALTEMPDIR/template${grp_temp}/finalwork

  if [[ ! -f $SRCWORKDIR/${id}_${side}_totemp_reslice_seg.nii.gz ]]; then
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $SRCWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $SRCWORKDIR/${id}_${side}_totemp_reslice_seg.nii.gz
  fi

  # Before we vote, use ml_affine for nice affine alignment
  ~pauly/wolk/ashs/ext/Linux/bin/ml_affine \
    $FINALWORKDIR/template_seg.nii.gz \
    $SRCWORKDIR/${id}_${side}_totemp_reslice_seg.nii.gz \
    $FINALWORKDIR/${id}_${side}_mlaffine.txt

  # Do we do ANTS or not?
  if [[ $doants -eq 0 ]]; then

    for sub in $KINDS; do

      c3d \
        $FINALWORKDIR/template_tse.nii.gz \
        $SRCWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz \
        -reslice-matrix \
        $FINALWORKDIR/${id}_${side}_mlaffine.txt \
        -o $FINALWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz

    done

  else

    # Convert that to ITK format
    c3d_affine_tool \
      $FINALWORKDIR/${id}_${side}_mlaffine.txt \
      -oitk $FINALWORKDIR/${id}_${side}_mlaffine_itk.txt

    CMD=""
    for sub in ${REG_LABELS_3[*]}; do
      CMD="$CMD -m MSQ[$FINALWORKDIR/template_${sub}.nii.gz,$SRCWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz,$WGT_3]"
    done

    if [[ $ANTs_x_3 == "Y" ]]; then
       ANTs_mask_3="-x $FINALWORKDIR/template_mask.nii.gz"
    else
       ANTs_mask_3=""
    fi

    ANTS 3 $CMD \
      -t $ANTs_t_3 \
      -r $ANTs_r_3 \
      -i $ANTs_i_3 \
      $ANTs_mask_3 \
      -a $FINALWORKDIR/${id}_${side}_mlaffine_itk.txt \
      --continue-affine 0 \
      $ANTs_all_metrics_3 \
      -o $FINALWORKDIR/${id}_${side}_totemp.nii.gz

    for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        $SRCWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz \
        $FINALWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz \
        -R $FINALWORKDIR/template_tse.nii.gz \
        $FINALWORKDIR/${id}_${side}_totempWarp.nii.gz \
        $FINALWORKDIR/${id}_${side}_totempAffine.txt

    done

    c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $FINALWORKDIR/${id}_${side}_totemp_reslice_seg.nii.gz

  fi
}

##############################################################################
function MakeFinalImages()
{
  for grp_temp in ${groups_temp[*]}; do
    FINALWORKDIR=$FINALTEMPDIR/template${grp_temp}/finalwork/

    # Generate one compound label for all other subfields
    ITDIR=$FINALWORKDIR/$(printf iter_%02d $((FINALITER-1)) )

    for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

      LABELDEF=(${EXTRAMESHESDEF[i]})
      echo $LABELDEF

      c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$ITDIR/template_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
        -accum -add -endaccum \
        -o $ITDIR/template_${EXTRAMESHES[i]}.nii.gz

    done
  done
}

##############################################################################
function WarpFinalLabel()
{
  PREFIX=WFL${expid}
  # Iterate over side and subject
  for side in left right; do

    IDS=$(cat $SUBJ_TXT)

    for id in $IDS; do

      #id=$(ls $ASHSRUNDIR | grep $id)

      for grp_temp in ${groups_temp[*]}; do

        mkdir -p $FINALTEMPDIR/template${grp_temp}/labelwarp
        # Submit job for this subject
        qsub -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -N "${PREFIX}_${id}_${side}_${grp_temp}" \
          $0 WarpFinalLabel_sub $id $side $grp_temp

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
    -q all.q,basic.q \
    -hold_jid "${PREFIX}_*" -sync y -b y \
    sleep 1
}

function WarpFinalLabel_sub()
{
  id=$1
  side=$2
  grp_temp=$3
  SEG=$CLEANUPDIR/${id}_${side}_dividedCS.nii.gz
  id=$(ls $ASHSRUNDIR | grep $id)

  FINALDIR=$FINALTEMPDIR/template${grp_temp}
  FINALWORKDIR=$FINALDIR/finalwork/

  # add on transform
  if [[ $side == "left" ]]; then
    ADD_ON_TRANS=$IDENTITYFN
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
  fi

  # Apply the warp to each of the labels
  for sub in ${LABEL_IDS[*]}; do

    antsApplyTransforms -d 3 \
      -r $SEG \
      -i $FINALWORKDIR/template_${sub}.nii.gz \
      -o $FINALDIR/labelwarp/${id}_${side}_${sub}_tempfit.nii.gz \
      -t [$ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt,1] \
      -t [$ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt,1] \
      -t $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii \
      -t [$ADD_ON_TRANS,1] \
      -t $FINALDIR/work/${id}_${side}_totempInverseWarp.nii.gz \
      -t [$FINALWORKDIR/${id}_${side}_totempAffine.txt,1] \
      -t $FINALWORKDIR/${id}_${side}_totempInverseWarp.nii.gz \
      -n BSpline

  done

  # Vote in original subject space to get the tempfit segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALDIR/labelwarp/${id}_${side}_${sub}_tempfit.nii.gz; done) -vote -type ushort \
    -o $FINALDIR/labelwarp/${id}_${side}_seg_tempfit.nii.gz

  # change ASHS labels to the current formate
  for ((i=0;i<${#LABEL_IDS[*]};i++)); do
    c3d $SEG \
      -replace $(for k in ${LABEL_MRG[i]}; do echo $k 999; done) \
      -thresh 999 999 1 0 \
      -o $FINALDIR/labelwarp/${id}_${side}_${LABEL_IDS[i]}_ashs.nii.gz
  done

  c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALDIR/labelwarp/${id}_${side}_${sub}_ashs.nii.gz; done) -vote -type ushort \
    -o $FINALDIR/labelwarp/${id}_${side}_seg_ashs.nii.gz

  for ((i=0;i<${#LABEL_IDS[*]};i++)); do
    rm $FINALDIR/labelwarp/${id}_${side}_${LABEL_IDS[i]}_ashs.nii.gz
    rm $FINALDIR/labelwarp/${id}_${side}_${LABEL_IDS[i]}_tempfit.nii.gz
  done

  # change manual segmentation to the current formate
  if [[ -d $ASHSRUNDIR/$id/refseg ]]; then
    for ((i=0;i<${#LABEL_IDS[*]};i++)); do
      c3d $ASHSRUNDIR/${id}/refseg/refseg_${side}.nii.gz \
        -replace $(for k in ${LABEL_MRG[i]}; do echo $k 999; done) \
        -thresh 999 999 1 0 \
        -o $FINALDIR/labelwarp/${id}_${side}_${LABEL_IDS[i]}_refseg.nii.gz
    done

    c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALDIR/labelwarp/${id}_${side}_${sub}_refseg.nii.gz; done) -vote -type ushort \
      -o $FINALDIR/labelwarp/${id}_${side}_seg_refseg.nii.gz

    for ((i=0;i<${#LABEL_IDS[*]};i++)); do
      rm $FINALDIR/labelwarp/${id}_${side}_${LABEL_IDS[i]}_refseg.nii.gz
    done
  fi
}

##############################################################################
function WarpFinalLabelToData()
{
  PREFIX=WFL${expid}
  # Iterate over side and subject
  for side in left right; do

    IDS=$(cat $SUBJ_TXT)

    for id in $IDS; do

      id=$(ls $ASHSRUNDIR | grep $id)

      for grp_temp in ${groups_temp[*]}; do

        mkdir -p $FINALTEMPDIR/template${grp_temp}/labelwarp_todata
        # Submit job for this subject
        qsub -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -N "${PREFIX}_${id}_${side}_${grp_temp}" \
          $0 WarpFinalLabelToData_sub $id $side $grp_temp

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
    -q all.q,basic.q \
    -hold_jid "${PREFIX}_*" -sync y -b y \
    sleep 1
}

function WarpFinalLabelToData_sub()
{
  id=$1
  side=$2
  grp_temp=$3

  FINALDIR=$FINALTEMPDIR/template${grp_temp}
  FINALWORKDIR=$FINALDIR/finalwork/

  # Apply the warp to each of the labels
  for sub in ${LABEL_IDS[*]}; do

    antsApplyTransforms -d 3 \
      -r $DATADIR/${id}_${side}_seg.nii.gz \
      -i $FINALWORKDIR/template_${sub}.nii.gz \
      -o $FINALDIR/labelwarp_todata/${id}_${side}_${sub}_tempfit.nii.gz \
      -t $FINALDIR/work/${id}_${side}_totempInverseWarp.nii.gz \
      -t [$FINALWORKDIR/${id}_${side}_totempAffine.txt,1] \
      -t $FINALWORKDIR/${id}_${side}_totempInverseWarp.nii.gz \
      -n BSpline

  done

  # Vote in original subject space to get the tempfit segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALDIR/labelwarp_todata/${id}_${side}_${sub}_tempfit.nii.gz; done) -vote -type ushort \
    -o $FINALDIR/labelwarp_todata/${id}_${side}_seg_tempfit.nii.gz

  for ((i=0;i<${#LABEL_IDS[*]};i++)); do
    rm $FINALDIR/labelwarp_todata/${id}_${side}_${LABEL_IDS[i]}_tempfit.nii.gz
  done
}

##############################################################################
function WarpFinalMesh()
{
  PREFIX=WFM${expid}
  for grp_temp in ${groups_temp[*]}; do
    FINALDIR=$FINALTEMPDIR/template${grp_temp}
    FINALWORKDIR=$FINALDIR/finalwork/

    mkdir -p $FINALDIR/meshwarp
    mkdir -p $FINALDIR/jacobian
    mkdir -p $FINALDIR/meshwarp/template

    ITDIR=$FINALWORKDIR/$(printf iter_%02d $((FINALITER-1)) )

    # Generate meshes for the individual subfields
    for sub in ${LABEL_FG[*]}; do

      vtklevelset \
        $ITDIR/template_${sub}.nii.gz \
        $FINALDIR/meshwarp/template_${sub}.vtk \
        $subfield_th

    done

    # Generate one mesh for all non-DG non-CS subfields
    for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

      vtklevelset \
        $ITDIR/template_${EXTRAMESHES[i]}.nii.gz \
        $FINALDIR/meshwarp/template_${EXTRAMESHES[i]}.vtk \
        $template_th

    done

    cp $FINALDIR/meshwarp/template*.vtk \
       $FINALDIR/meshwarp/template/

    IDS=$(cat $SUBJ_TXT)

    for id in $IDS; do
 
      id=$(ls $ASHSRUNDIR | grep $id)

      for side in left right; do 

        # Submit job for this subject
        qsub  -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -N "${PREFIX}_${id}_${side}_${grp_temp}" \
          $0 WarpFinalMesh_sub $id $side $grp_temp

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
     -q all.q,basic.q \
     -hold_jid "${PREFIX}_*" -sync y -b y \
     sleep 1
}

function WarpFinalMesh_sub()
{
  id=$1
  side=$2
  grp_temp=$3
  ALLSF="${LABEL_FG[*]} ${EXTRAMESHES[*]}"
  FINALDIR=$FINALTEMPDIR/template${grp_temp}
  FINALWORKDIR=$FINALDIR/finalwork/

  # add on transform
  if [[ $side == "left" ]]; then
    ADD_ON_TRANS=$IDENTITYFN
  else
    ADD_ON_TRANS=$DATADIR/tempreg/right_to_left_init_affine.txt
  fi

  # Compose the transformation between the template and the subject
  $ANTSPATH/ComposeMultiTransform 3 \
     $TMPDIR/compose.nii \
     -R $FINALWORKDIR/template_tse.nii.gz \
     $FINALWORKDIR/${id}_${side}_totempWarp.nii.gz \
     $FINALWORKDIR/${id}_${side}_totempAffine.txt \
     $FINALDIR/work/${id}_${side}_totempWarp.nii.gz \
     $ADD_ON_TRANS \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempWarp.nii \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
     $ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

  # Split the transformations into xyz
  c3d -mcs $TMPDIR/compose.nii -oo $TMPDIR/comp%02d.nii

  # Apply the warp to each of the meshes
  for sub in $ALLSF; do

    # Warp the subfield into subject space
    warpmesh -w ants -m ras \
      $FINALDIR/meshwarp/template_${sub}.vtk \
      $FINALDIR/meshwarp/${id}_${side}_${sub}_tempfit.vtk \
      $TMPDIR/comp00.nii $TMPDIR/comp01.nii $TMPDIR/comp02.nii

    # Extract the thickness of the subfield
    cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
      -T $FINALDIR/meshwarp/${id}_${side}_${sub}_thickmap.vtk \
      -p $thick_p -e $thick_e \
      $FINALDIR/meshwarp/${id}_${side}_${sub}_tempfit.vtk \
      $TMPDIR/meshwarp/skel_${id}_${side}_${sub}.vtk

  done

  # Compute the Jacobian map
  #ANTSJacobian 3 $TMPDIR/compose.nii \
  #  $FINALDIR/jacobian/${id}_${side}_totemp 1
}

##############################################################################
function EvalFinalTemp()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=EFT${expid}
  # Iterate over side and subject
  for grp_temp in ${groups_temp[*]}; do

    FINALEVALTMPDIR=$FINALTEMPDIR/template${grp_temp}/evaluation/tmp
    rm -rf $FINALEVALTMPDIR
    mkdir -p $FINALEVALTMPDIR

    for side in left right; do
      for id in $IDS; do

        id=$(ls $ASHSRUNDIR | grep $id)
        # Submit job for this subject
        qsub -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -N "${PREFIX}_${id}_${side}_${grp_temp}" \
          $0 EvalFinalTemp_sub \
          $id $side $grp_temp $FINALEVALTMPDIR
        sleep 0.1

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # Combine scores
  for grp_temp in ${groups_temp[*]}; do

    FINALEVALDIR=$FINALTEMPDIR/template${grp_temp}/evaluation
    FINALEVALTMPDIR=$FINALEVALDIR/tmp

    echo "ID ${EVALLABELS[*]} ${MESH_EVAL[*]} ${EVALLABELS[*]}" \
        > $FINALEVALDIR/overlap.txt

    echo "ID ${EVALLABELS[*]} ${MESH_EVAL[*]}" \
        > $FINALEVALDIR/overlap_refseg.txt

    for side in left right; do

      echo "ID ${EVALLABELS[*]} ${MESH_EVAL[*]} ${EVALLABELS[*]}" \
        > $FINALEVALDIR/overlap_${side}.txt

      echo "ID ${EVALLABELS[*]} ${MESH_EVAL[*]}" \
        > $FINALEVALDIR/overlap_refseg_${side}.txt

      for id in $IDS; do

        id=$(ls $ASHSRUNDIR | grep $id)
        cat $FINALEVALTMPDIR/${id}_${side}_overlap.txt \
            >> $FINALEVALDIR/overlap_${side}.txt

        cat $FINALEVALTMPDIR/${id}_${side}_overlap_combined.txt \
            >> $FINALEVALDIR/overlap.txt

        if [[ -f $FINALEVALTMPDIR/${id}_${side}_refseg_overlap.txt ]]; then
          cat $FINALEVALTMPDIR/${id}_${side}_refseg_overlap.txt \
            >> $FINALEVALDIR/overlap_refseg_${side}.txt

          cat $FINALEVALTMPDIR/${id}_${side}_refseg_overlap_combined.txt \
            >> $FINALEVALDIR/overlap_refseg.txt

        fi

      done
    done
  done
}

function EvalFinalTemp_sub()
{
  id=$1
  side=$2
  grp_temp=$3
  EVALTMPDIR=$4
  ALLOVL=""
  FINALDIR=$FINALTEMPDIR/template${grp_temp}
  FINALWORKDIR=$FINALDIR/finalwork/

  ###################################
  # Compute dice overlap in subject space
  FITTED=$FINALDIR/labelwarp/${id}_${side}_seg_tempfit.nii.gz
  ASHSSEG=$FINALDIR/labelwarp/${id}_${side}_seg_ashs.nii.gz
  do_pair $ASHSSEG $FITTED
  ALLOVL="$ALLOVL $FULLOVL"

  ###################################
  # Compute H distance
  MESH_DIST=""
  for ((i=0;i<${#MESH_EVAL[*]};i++)); do

    c3d $ASHSSEG \
      -replace $(for k in ${MESH_EVAL_RANGES[i]}; do echo $k 999; done) \
      -thresh 999 999 1 0 \
      -trim 2mm -resample 100x100x300% \
      -o $TMPDIR/mask_seg.nii.gz

    vtklevelset $TMPDIR/mask_seg.nii.gz $TMPDIR/mask_seg.vtk 0.5
    SIZE=$(cat $TMPDIR/mask_seg.vtk | wc -l)
    if [[ $SIZE -gt 6 ]]; then

      MESH_DIST_TMP="$(meshdiff $FINALDIR/meshwarp/${id}_${side}_${MESH_EVAL[i]}_tempfit.vtk $TMPDIR/mask_seg.vtk | grep RESULT | awk '{print $7}')"

      if [[ $MESH_DIST_TMP == "" ]]; then
        MESH_DIST="$MESH_DIST NA"
      else
        MESH_DIST="$MESH_DIST $MESH_DIST_TMP"
      fi

    else

      MESH_DIST="$MESH_DIST NA"

    fi

  done

  ALLOVL="$ALLOVL $MESH_DIST"

  ###################################
  # Compute dice overlap in template space
  do_pair $FINALWORKDIR/template_seg.nii.gz \
    $FINALWORKDIR/${id}_${side}_totemp_reslice_seg.nii.gz
  ALLOVL="$ALLOVL $FULLOVL"

  ###################################
  # save to temporary file
  echo $id $ALLOVL > $EVALTMPDIR/${id}_${side}_overlap.txt
  echo $id $side $ALLOVL > $EVALTMPDIR/${id}_${side}_overlap_combined.txt

  ###################################
  ALLOVL=""
  # Compute dice overlap between warp seg and manual seg
  FITTED=$FINALDIR/labelwarp/${id}_${side}_seg_tempfit_ref.nii.gz
  REFSEG=$FINALDIR/labelwarp/${id}_${side}_seg_refseg.nii.gz
  if [[ -f $REFSEG ]]; then
  # compute dice
  c3d $FINALDIR/labelwarp/${id}_${side}_seg_tempfit.nii.gz -replace 5 6 -o $FITTED
  do_pair $REFSEG $FITTED
  ALLOVL="$ALLOVL $FULLOVL"

  # Compute H distance
  MESH_DIST=""
  for ((i=0;i<${#MESH_EVAL[*]};i++)); do

    c3d $REFSEG \
      -replace $(for k in ${MESH_EVAL_RANGES[i]}; do echo $k 999; done) \
      -thresh 999 999 1 0 \
      -trim 2mm -resample 100x100x300% \
      -o $TMPDIR/mask_seg.nii.gz

    vtklevelset $TMPDIR/mask_seg.nii.gz $TMPDIR/mask_seg.vtk 0.5
    SIZE=$(cat $TMPDIR/mask_seg.vtk | wc -l)
    if [[ $SIZE -gt 6 ]]; then

      MESH_DIST_TMP="$(meshdiff $FINALDIR/meshwarp/${id}_${side}_${MESH_EVAL[i]}_tempfit.vtk $TMPDIR/mask_seg.vtk | grep RESULT | awk '{print $7}')"

      if [[ $MESH_DIST_TMP == "" ]]; then
        MESH_DIST="$MESH_DIST NA"
      else
        MESH_DIST="$MESH_DIST $MESH_DIST_TMP"
      fi

    else

      MESH_DIST="$MESH_DIST NA"

    fi

  done

  ALLOVL="$ALLOVL $MESH_DIST"

  # save to temporary file
  echo $id $ALLOVL > $EVALTMPDIR/${id}_${side}_refseg_overlap.txt
  echo $id $side $ALLOVL > $EVALTMPDIR/${id}_${side}_refseg_overlap_combined.txt
  fi

}

##############################################################################
function EvalFinalTempData()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=EFT${expid}
  # Iterate over side and subject
  for grp_temp in ${groups_temp[*]}; do

    FINALEVALTMPDIR=$FINALTEMPDIR/template${grp_temp}/evaluation_data/tmp
    rm -rf $FINALEVALTMPDIR
    mkdir -p $FINALEVALTMPDIR

    for side in left right; do
      for id in $IDS; do

        id=$(ls $ASHSRUNDIR | grep $id)
        # Submit job for this subject
        qsub -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -N "${PREFIX}_${id}_${side}_${grp_temp}" \
          $0 EvalFinalTempData_sub \
          $id $side $grp_temp $FINALEVALTMPDIR
        sleep 0.1

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # Combine scores
  for grp_temp in ${groups_temp[*]}; do

    FINALEVALDIR=$FINALTEMPDIR/template${grp_temp}/evaluation_data
    FINALEVALTMPDIR=$FINALEVALDIR/tmp

    echo "ID ${EVALLABELS[*]}" \
        > $FINALEVALDIR/overlap.txt

    for side in left right; do

      echo "ID ${EVALLABELS[*]}" \
        > $FINALEVALDIR/overlap_${side}.txt

      for id in $IDS; do

        id=$(ls $ASHSRUNDIR | grep $id)
        cat $FINALEVALTMPDIR/${id}_${side}_overlap.txt \
            >> $FINALEVALDIR/overlap_${side}.txt

        cat $FINALEVALTMPDIR/${id}_${side}_overlap_combined.txt \
            >> $FINALEVALDIR/overlap.txt

      done
    done
  done
}

function EvalFinalTempData_sub()
{
  id=$1
  side=$2
  grp_temp=$3
  EVALTMPDIR=$4
  ALLOVL=""
  FINALDIR=$FINALTEMPDIR/template${grp_temp}
  FINALWORKDIR=$FINALDIR/finalwork/

  ###################################
  # Compute dice overlap in subject space
  FITTED=$FINALDIR/labelwarp_todata/${id}_${side}_seg_tempfit.nii.gz
  ASHSSEG=$DATADIR/${id}_${side}_seg.nii.gz
  do_pair $ASHSSEG $FITTED
  ALLOVL="$ALLOVL $FULLOVL"

  echo $id $ALLOVL > $EVALTMPDIR/${id}_${side}_overlap.txt
  echo $id $side $ALLOVL > $EVALTMPDIR/${id}_${side}_overlap_combined.txt
}

##############################################################################
function DispFinalStat()
{
  ALLSF="${LABEL_FG[*]} ${EXTRAMESHES[*]}"
  PREFIX=DFS${expid}
  for grp_temp in ${groups_temp[*]}; do

    FINALDIR=$FINALTEMPDIR/template${grp_temp}/
    rm -rf $FINALDIR/meshwarp/analysis
    mkdir -p $FINALDIR/meshwarp/analysis

    for side in left right; do
      for sub in $ALLSF; do

        # Submit job for this subject
        qsub -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -N "${PREFIX}_${side}_${grp_temp}_${sub}" \
          $0 DispFinalStat_sub $side $grp_temp $sub

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function DispFinalStat_sub()
{
  side=$1
  grp_temp=$2
  sub=$3

  IDS=$(cat $SUBJ_TXT)
  FINALDIR=$FINALTEMPDIR/template${grp_temp}/

  # Displacement analysis
  MESHES=$(for id in $IDS; do \
    echo $(find $FINALDIR/meshwarp | grep $id | grep $side | grep ${sub}_tempfit ); done)

  meshdisp $MESHES \
    $FINALDIR/meshwarp/analysis/disp_${side}_${sub}.vtk

  # Thickness analysis
  MESHES=$(for id in $IDS; do \
    echo $(find $FINALDIR/meshwarp | grep $id | grep $side | grep ${sub}_thickmap ); done)

  mesh_merge_arrays \
    -r $FINALDIR/meshwarp/template_${sub}.vtk \
    $FINALDIR/meshwarp/analysis/thick_${side}_${sub}.vtk \
    Thickness $MESHES
}

##############################################################################
#ANGRP[0]="MRG DG"
#ANGRP[1]="CA DG SUB ERC BA35 BA36"

#GRPNM[0]="merged"
#GRPNM[1]="all"
function ThickFinalStat()
{
  PREFIX=TFS${expid}
  for grp_temp in ${groups_temp[*]}; do

    rm -rf $FINALTEMPDIR/template${grp_temp}/meshwarp/analysis/design*

    #for side in left right; do
     side="left"
     for design in $(ls $ALLSTATDIR | grep "design_.*txt" | grep $side ); do

        exp=$(echo $design | sed -e "s/^design_${side}_//" | sed -e "s/\.txt//")

        for ((igrp=0;igrp<${#ANGRP[*]};igrp++)); do

          # Submit job for this subject
          qsub -cwd -o $DUMPDIR -j y \
            -q all.q,basic.q \
            -N "${PREFIX}_${exp}_${grp_temp}_${GRPNM[igrp]}" \
            $0 ThickFinalStat_sub $exp $grp_temp $igrp

        done
      done
    #done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function ThickFinalStat_sub()
{
  exp=$1
  grp_temp=$2
  igrp=$3

  MYGRP=${ANGRP[igrp]}
  GNAME=${GRPNM[igrp]}

  FINALDIR=$FINALTEMPDIR/template${grp_temp}/

  # Create the work directory for this analysis
  WORK=$FINALDIR/meshwarp/analysis/design_${exp}_group_${GNAME}
  mkdir -p $WORK

  # both side
  for side in left right; do
    # Merge the meshes for this analysis
    for sub in $MYGRP; do

      # Get the list of subjects
      DESIGNTXT=$ALLSTATDIR/design_${side}_${exp}.txt
      SUBJ=$(cat $DESIGNTXT | awk '{print $1}')

      # Generate the design matrix for meshglm
      cat $DESIGNTXT | awk '{$1=""; print}' > $WORK/design_${side}.txt

      MESHES=$(for id in $SUBJ; do \
        echo $(find $FINALDIR/meshwarp | grep ${id} | grep ${side}_${sub}_thickmap.vtk); done)

      mesh_merge_arrays -r \
        $FINALDIR/meshwarp/template_${sub}.vtk \
        $WORK/thick_${side}_${grp_temp}_${sub}.vtk Thickness $MESHES
    done
  done

  # Go through the list of contrasts
  for con in $(ls $ALLSTATDIR | grep "contrast_${exp}_.*\.txt"); do

    # Get the suffix for this contrast
    suffix=$(echo $con | sed -e "s/^contrast_${exp}_//" -e "s/\.txt//")

    # Create the directory for this contrast
    CWORK=$WORK/contrast_${suffix}
    mkdir -p $CWORK

    FULLNM="design_${exp}_group_${GNAME}_con_${suffix}"

    # Copy the contrast
    cp $ALLSTATDIR/$con $CWORK/contrast.txt

    # Build the list of meshes to include
    for sub in $MYGRP; do

      # pool both hemisphere for multiple comparisons
      MESHPARAM=""
      for side in left right; do
        MESHPARAM="$MESHPARAM -m $WORK/thick_${side}_${grp_temp}_${sub}.vtk $CWORK/thickstat_${FULLNM}_${side}_${grp_temp}_${sub}.vtk"
      done

      meshglm $MESHPARAM \
        -g $WORK/design_${side}.txt $CWORK/contrast.txt \
        -a Thickness \
        -d 4 \
        -p 1000 \
        -s P \
        -t 0.01 \
        -e

    done
  done
}

##############################################################################
function FinalMeanThickness()
{
  PREFIX=FMT${expid}
  for grp_temp in ${groups_temp[*]}; do
    MESHGRPDIR=$FINALTEMPDIR/template${grp_temp}/meshwarp
    rm -rf $MESHGRPDIR/analysis/labels
    mkdir -p $MESHGRPDIR/analysis/labels
    for side in left right; do

      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=2.1G,s_vmem=2G \
           -N "${PREFIX}_${grp_temp}_${side}" \
           $0 FinalMeanThickness_sub $grp_temp $side
      sleep 0.1
      #FinalMeanThickness_sub $grp_temp $side

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # generate header
  header="ID"
  for side in left right; do
    for ((i=0;i<${#MESH_LABEL[*]};i++)); do
      sub=${MESH_LABEL[$i]}
      header="${header},${sub}_${side}"
    done
    header="${header},aBA35_${side},pBA35_${side}"
  done

  # generate csv file with summary measurements!
  IDS=($(cat $SUBJ_TXT))
  for grp_temp in ${groups_temp[*]}; do
    FINALTEMPGRPDIR=$FINALTEMPDIR/template${grp_temp}
    rm -rf $FINALTEMPGRPDIR/analysis
    mkdir -p $FINALTEMPGRPDIR/analysis

    # save header
    echo "$header" > $FINALTEMPGRPDIR/analysis/summary_thick_${grp_temp}.csv

    # combine files
    for ((i=0;i<${#IDS[*]};i++)); do

      # id
      id=${IDS[$i]}

      # load left side
      side=left
      THICK_LEFT=$(cat $FINALTEMPGRPDIR/meshwarp/analysis/labels/thick_${side}_${grp_temp}_MRG.txt | grep $id)

      # load right side
      side=right
      THICK_RIGHT=$(cat $FINALTEMPGRPDIR/meshwarp/analysis/labels/thick_${side}_${grp_temp}_MRG.txt | grep $id | cut -d ' ' -f 2)

      echo "${THICK_LEFT}${THICK_RIGHT}" >> \
        $FINALTEMPGRPDIR/analysis/summary_thick_${grp_temp}.csv

    done

    # perform statistical analysis
    Rscript $WORKDIR/scripts/summary_thickness_analysis.R \
      $WORKDIR/analysis_input/demog.csv \
      $FINALTEMPGRPDIR/analysis/summary_thick_${grp_temp}.csv \
      $FINALTEMPGRPDIR/analysis/mean_thickness_stats${grp_temp}.csv

  done
}

function FinalMeanThickness_sub()
{
  grp_temp=$1
  side=$2

  # sample on probability maps (don`t run BKG)
  FINALTEMPGRPDIR=$FINALTEMPDIR/template${grp_temp}
  MESHGRPDIR=$FINALTEMPGRPDIR/meshwarp
  WORKGRPDIR=$FINALTEMPGRPDIR/finalwork
  for ((ii=0;ii<${#EXTRAMESHES[*]};ii++)); do
    LABELDEF=(${EXTRAMESHESDEF[ii]})
    MESHES=""
    for ((i=0;i<${#LABEL_IDS[*]};i++)); do
      if [[ ${LABELDEF[$i]} == "1" ]]; then
        sub=${LABEL_IDS[$i]}
        mesh_image_sample \
          $MESHGRPDIR/template_${EXTRAMESHES[ii]}.vtk \
          $WORKGRPDIR/template_${sub}.nii.gz \
          $MESHGRPDIR/analysis/labels/template_${side}_${grp_temp}_${EXTRAMESHES[ii]}_${sub}.vtk \
          PROB
        MESHES="$MESHES $MESHGRPDIR/analysis/labels/template_${side}_${grp_temp}_${EXTRAMESHES[ii]}_${sub}.vtk"
      fi
    done

    # merge prob arrays
    mesh_merge_arrays \
      -r $MESHGRPDIR/analysis/thick_${side}_${EXTRAMESHES[ii]}.vtk \
      $MESHGRPDIR/analysis/labels/thick_${side}_${grp_temp}_${EXTRAMESHES[ii]}.vtk \
      PROB $MESHES

    # run matlab script to generate label
      $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
        addpath('/data/picsl/longxie/ASHS/thickness_newlabel/matlabcode');
        compute_mesh_label('$MESHGRPDIR/analysis/labels/thick_${side}_${grp_temp}_${EXTRAMESHES[ii]}.vtk');
MATCODE

  done

  # run matlab script to generate label
  #mkdir -p $TMPDIR/${id}
  source $PKGDIR/matlab_batch.sh \
    $TMPDIR/${id}/ \
    $MATLABCODEDIR \
    compute_mean_thickness_apBA35 \
    $MESHGRPDIR/analysis/labels/thick_${side}_${grp_temp}_MRG.vtk \
    2 \
    $WORKDIR/analysis_input/subj_vertical.txt
}

##############################################################################
function reset_dir()
{
  rm -rf $DUMPDIR/*
}

##############################################################################
if [[ $# -lt 2 ]]; then

  main

else

  cmd=$1
  shift
  $cmd $@

fi
