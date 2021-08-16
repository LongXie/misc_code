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
expid=800
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

# similarity type
SIM_TYPE="PRCCS_seg_dice"
SIM_DIR=$expdir/PWsim/sim_${SIM_TYPE}

# 1.2 path directory
ALLPATHDIR=$expdir/MST/paths

########################################
# 2. MST registration
MSTREGDIR=$expdir/MST/registration
MSTTEMPDIR=$expdir/MST/template

# 2.2 generate template
EXTRAMESHES=(PRC Anterior MRG)
#                BKG ERC BA35 BA36 PHC aCS pCS
EXTRAMESHESDEF=("-1  -1   1    1   -1  -1  -1" \
                "-1   1   1    1   -1  -1  -1" \
                "-1   1   1    1    1  -1  -1")



########################################
# 3. shooting
SHOOTDIR=$expdir/gshoot
GSHOOT_NITER=5







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
  #PreprocessLabels
  #copy_data

  #######################################
  # 1 perform MST
  # 1.1 pairwise registration
  #pairwise
  #completeness  

  # 1.2 find paths between all subjects and the root
  #PrintAdjacency

  #######################################
  # 2. build template using MST
  # 2.1 registration through MST
  #RegisterMST
  
  # 2.2 generate initial template and generate landmarks
  #GenerateInitTemp

  #######################################
  # 3. refine template using geodesic shooting
  # 3.1 perform shape correction
  #ShootingCorrection

  #######################################
  # 4. compute mean and mode
  #ShapeStats

  MakeMovies
    


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

############################################################
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
        #/data/picsl/pauly/wolk/ashs/ext/Linux/bin/ml_affine \
        #  $DATADIR/${fnside_fix}_seg.nii.gz \
        #  $DATADIR/${fnside_mov}_seg.nii.gz \
        #  $TMPDIR/${fnside_mov}_to_${fnside_fix}_mlaffine.txt

        # Convert that to ITK format
        #c3d_affine_tool \
        #  $TMPDIR/${fnside_mov}_to_${fnside_fix}_mlaffine.txt \
        #  -oitk \
        #  $TMPDIR/${fnside_mov}_to_${fnside_fix}_mlaffine_itk.txt

        # file names
        MAT_MOMENTS=$TMPDIR/${fnside_mov}_to_${fnside_fix}_moment.mat
        MAT_AFFINE=$TMPDIR/${fnside_mov}_to_${fnside_fix}_affine.mat
        WARP=$TMPDIR/${fnside_mov}_to_${fnside_fix}_warp.mat

        # greedy command
        CMD=""
        for sub in ${REG_LABELS_1[*]}; do
          CMD="$CMD -w 1 -i $DATADIR/${fnside_fix}_${sub}.nii.gz $DATADIR/${fnside_mov}_${sub}.nii.gz "
        done

        # Perform moments of intertia matching between the two masks
        greedy -d 3 -threads $NSLOTS \
          $CMD \
          -moments \
          -o $MAT_MOMENTS

        # Perform affine matching between the two masks
        greedy -d 3 -threads $NSLOTS \
          $CMD \
         -a -ia $MAT_MOMENTS \
         -n 100x100 \
         -o $MAT_AFFINE

        # Run greedy between these two images
        if [[ $ANTs_x_1 == "Y" ]]; then
          c3d $DATADIR/${fnside_fix}_seg.nii.gz -dup \
            $DATADIR/${fnside_mov}_seg.nii.gz \
            -int 0 -reslice-identity \
            -add -binarize -dilate 1 10x10x10vox \
            -o $TMPDIR/${fnside_mov}_to_${fnside_fix}_mask.nii.gz
          MASK="-gm $TMPDIR/${fnside_mov}_to_${fnside_fix}_mask.nii.gz"
        else
          MASK=""
        fi

        greedy -d 3 -threads $NSLOTS \
          $GREEDY_INPUTS \
          -it $MAT_AFFINE \
          $MASK \
          -n 50x50x20x0 \
          -s 2.0mm 0.1mm -e 0.5 \
          -o $WARP
 
        # Reslice the segmentations from raw space
        RM=""
        for sub in ${KINDS[*]}; do

          RM="$RM -rm $DATADIR/${fnside_mov}_${sub}.nii.gz $TMPDIR/${fnside_mov}_to_${fnside_fix}_reslice_${sub}.nii.gz"

        done

        greedy -d 3 \
          -rf $DATADIR/${fnside_fix}_seg.nii.gz \
          $RM \
          -r $WARP $MAT_AFFINE

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

##########################################################
function PrintAdjacency()
{
  IDS=$(cat $SUBJ_TXT)
  rm -rf $ALLPATHDIR
  mkdir -p $ALLPATHDIR

  # create adjacency matrix
  if [ ! -f $ALLPATHDIR/adj.txt ]; then
  idx_fix=1
  for side_fix in left right; do
  for id_fix in $IDS; do
  idx_mov=1
  for side_mov in left right; do
  for id_mov in $IDS; do

    idside_fix=${id_fix}_${side_fix}
    fn_fix=$(ls $ASHSRUNDIR | grep $id_fix)
    fnside_fix=${fn_fix}_${side_fix}
    idside_mov=${id_mov}_${side_mov}
    fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)
    fnside_mov=${fn_mov}_${side_mov}
    PWREGDIR=$PWDIR/${fnside_fix}/${fnside_mov}_to_${fnside_fix}
    SIMFILE=$PWREGDIR/${fnside_mov}_to_${fnside_fix}_sim.txt

    if [[ $idside_fix != $idside_mov ]]; then
      if [[ -f $SIMFILE ]]; then
        echo "$idx_fix $idx_mov $(cat $SIMFILE)" \
          >> $ALLPATHDIR/adj.txt
      fi
    fi

    idx_mov=$(($idx_mov+1))

  done
  done
  idx_fix=$(($idx_fix+1))
  echo "$fnside_fix" >> $ALLPATHDIR/IDSide.txt
  done
  done
  fi

  /share/apps/R/R-3.1.1/bin/Rscript compute_mst.R $ALLPATHDIR/adj.txt \
    > $ALLPATHDIR/mst_paths.txt
}

############################################################
function RegisterMST()
{
  # id
  IDS=$(cat $SUBJ_TXT)
  #IDS="DW124"

  # submit job for each subject
  PREFIX=RMST${expid}
  idx=1
  for side in left right; do
    for id in $IDS; do

      qsubp2 -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${id}_${side}" $0 \
           RegisterMST_sub $idx
      sleep 0.1
      idx=$((${idx}+1))

    done

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function RegisterMST_sub()
{
  idxmov=$1
  IDS=$(cat $SUBJ_TXT)
  fnside_mov=$(cat $ALLPATHDIR/IDSide.txt | head -n $idxmov | tail -n 1)

  # get path
  GPATH=($(cat $ALLPATHDIR/mst_paths.txt | head -n $idxmov | tail -n 1))

  # We start with the warp chain as an empty list
  WARPCHAIN=""
  
  # The id in whose space the moving image is - at the beginning this is the
  # same as the subject's ID
  IDMOV=$fnside_mov

  # Loop. The first element of GPATH is the image itself
  for ((j=1;j<${#GPATH[*]};++j)); do

    # The next image in the path
    idxref=${GPATH[j]}
    IDREF=$(cat $ALLPATHDIR/IDSide.txt | head -n $idxref | tail -n 1)

    # Create directory for this registration
    WORK=$MSTREGDIR/$fnside_mov/step_${j}_${IDMOV}_to_${IDREF}
    mkdir -p $WORK

    # file names
    PREFIX=$WORK/${IDMOV}_to_${IDREF}
    MAT_MOMENTS=${PREFIX}_moment.mat
    MAT_AFFINE=${PREFIX}_affine.mat
    WARP=${PREFIX}_warp.nii.gz

    if [[ ! -f $MAT_AFFINE || ! -f $WARP ]]; then

    # greedy command
    CMD=""
    for sub in ${REG_LABELS_1[*]}; do
      CMD="$CMD -w 1 -i $DATADIR/${IDREF}_${sub}.nii.gz $DATADIR/${IDMOV}_${sub}.nii.gz "
    done

    # Perform moments of intertia matching between the two masks
    greedy -d 3 -threads $NSLOTS \
      $CMD \
      -moments \
      -o $MAT_MOMENTS

    # Perform affine matching between the two masks
    greedy -d 3 -threads $NSLOTS \
      $CMD \
     -a -ia $MAT_MOMENTS \
     -n 100x100 \
     -o $MAT_AFFINE

    # Run greedy between these two images
    c3d $DATADIR/${IDREF}_seg.nii.gz -dup \
      $DATADIR/${IDMOV}_seg.nii.gz \
      -int 0 -reslice-identity \
      -add -binarize -dilate 1 10x10x10vox \
      -o $TMPDIR/${IDMOV}_to_${IDREF}_mask.nii.gz

    greedy -d 3 -threads $NSLOTS \
      $CMD \
      -it $MAT_AFFINE \
      -gm $TMPDIR/${IDMOV}_to_${IDREF}_mask.nii.gz \
      -n 50x40x20 -float \
      -s 0.6mm 0.1mm -e 0.5 \
      -o $WARP

    fi

    # Update the warp chain
    WARPCHAIN="$WARP $MAT_AFFINE $WARPCHAIN"
    IDMOV=${IDREF}
    
    # link the reference image
    for sub in ${KINDS[*]} seg; do
      ln -sf $DATADIR/${IDREF}_${sub}.nii.gz \
        $WORK/${IDREF}_${sub}.nii.gz
    done

    # Reslice the segmentations from raw space
    RM=""
    for sub in ${KINDS[*]}; do

      RM="$RM -rm $DATADIR/${fnside_mov}_${sub}.nii.gz $WORK/${fnside_mov}_to_${IDREF}_reslice_${sub}.nii.gz"

    done

    greedy -d 3 \
      -rf $DATADIR/${IDREF}_seg.nii.gz \
      $RM \
      -r $WARPCHAIN

    # Create seg
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $WORK/${fnside_mov}_to_${IDREF}_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $WORK/${fnside_mov}_to_${IDREF}_reslice_seg.nii.gz

  done

  # make link to the final directory
  FINALWORK=$MSTREGDIR/$fnside_mov/final
  mkdir -p $FINALWORK
  if [[ ${#GPATH[*]} == 1 ]]; then
    # for root subject, just link original data
    for sub in ${KINDS[*]} seg; do
      ln -sf $DATADIR/${fnside_mov}_${sub}.nii.gz \
        $FINALWORK/${fnside_mov}_to_MSTRoot_reslice_${sub}.nii.gz 
    done
  else
    for sub in ${KINDS[*]} seg; do
      ln -sf $WORK/${fnside_mov}_to_${IDREF}_reslice_${sub}.nii.gz \
        $FINALWORK/${fnside_mov}_to_MSTRoot_reslice_${sub}.nii.gz
    done
  fi
  echo $WARPCHAIN > $FINALWORK/chain_unwarp_to_final.txt
}

##########################################################
function GenerateInitTemp()
{
  mkdir -p $MSTTEMPDIR
 
  # create mean image for different labels
  PREFIX=GIT${expid}
  for sub in ${LABEL_IDS[*]}; do

    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX}_${sub}" $0 \
         GenerateTemp_sub $sub INIT
    sleep 0.1

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # generate seg
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $MSTTEMPDIR/template_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $MSTTEMPDIR/template_seg.nii.gz
  
  # generate meshes
  for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

    LABELDEF=(${EXTRAMESHESDEF[i]})
    echo $LABELDEF

    c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$MSTTEMPDIR/template_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
      -accum -add -endaccum \
      -o $MSTTEMPDIR/template_${EXTRAMESHES[i]}.nii.gz 

    vtklevelset $MSTTEMPDIR/template_${EXTRAMESHES[i]}.nii.gz \
      $MSTTEMPDIR/template_${EXTRAMESHES[i]}.stl 0.5

  done 
  
  # subsample mesh
  /data/picsl-build/pauly/vcg/gcc64rel/mesh_poisson_sample \
    $MSTTEMPDIR/template_MRG.stl \
    $MSTTEMPDIR/template_MRG_sampled.ply \
    2000

  # convert into vtk mesh
  TEMPLATE=$MSTTEMPDIR/template_MRG_sampled.vtk
  NV=$(cat $MSTTEMPDIR/template_MRG_sampled.ply | grep 'element vertex' | awk '{print $3}')

  # Write the header of the VTK file
  echo "# vtk DataFile Version 3.0" > $TEMPLATE
  echo "vtk output" >> $TEMPLATE
  echo "ASCII" >> $TEMPLATE
  echo "DATASET POLYDATA" >> $TEMPLATE
  echo "POINTS ${NV} float" >> $TEMPLATE

  # Write the point coordinates
  grep -A $NV end_header $MSTTEMPDIR/template_MRG_sampled.ply | tail -n $NV >> $TEMPLATE
}

function GenerateTemp_sub()
{
  sub=$1
  MODE=$2

  if [[ $MODE == "INIT" ]]; then
    c3d $MSTREGDIR/*/final/*reslice_${sub}.nii.gz \
      -mean -o $MSTTEMPDIR/template_${sub}.nii.gz
  else
    c3d $SHOOTDIR/*/iter_$((GSHOOT_NITER-1))/reslice*${sub}.nii.gz \
      -mean -o $SHOOTDIR/template/template_gshoot_${sub}.nii.gz
  fi
}

##########################################################
function ShootingCorrection()
{
  START_ITER=0
  END_ITER=$GSHOOT_NITER
  LANDMARKS=$MSTTEMPDIR/template_MRG_sampled.vtk
  IDS=$(cat $SUBJ_TXT)
  #IDS="DW101"
  mkdir -p $SHOOTDIR

  # Loop
  for ((iter=$START_ITER;iter<$END_ITER;iter++)); do

    # loop through all subjects
    PREFIX=SCI${expid}
    i=1
    for side in left right; do
    for id in $IDS; do
    
      # submit jobs
      qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=6.1G,s_vmem=6G \
         -N "${PREFIX}_${i}_${iter}" $0 \
         ShootingCorrectionSubj_sub $i $iter $LANDMARKS
      sleep 0.1
      i=$((i+1))

    done
    done

    # Wait for completion
    qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

    # perform averaging
    PREFIX=SCA${expid}
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=6.1G,s_vmem=6G \
         -N "${PREFIX}_shape_${iter}" $0 \
         ShootingCorrectionShapeAvg_sub $iter $LANDMARKS

    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=6.1G,s_vmem=6G \
         -N "${PREFIX}_volume_${iter}" $0 \
         ShootingCorrectionVolumeAvg_sub $iter

    # Wait for completion
    qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1


    # landmark for the next iteration
    LANDMARKS=$SHOOTDIR/shape_avg/iter_$((iter))/shavg_landmarks.vtk

  done
  
}

function ShootingCorrectionSubj_sub()
{
  # The index of the subject being registered using the MST
  idx=${1?}
  ID=$(cat $ALLPATHDIR/IDSide.txt | head -n $idx | tail -n 1)

  # The iteration - 0 is the initial iteration, involves extra work
  iter=${2?}

  # The path to the landmarks
  LANDMARKS=${3?}

  # directory
  WORK=$SHOOTDIR/$ID
  WITER=$WORK/iter_${iter}
  mkdir -p $WITER
  
  # Reference space (root node in cm-rep space)
  REFSPACE=$DATADIR/${ID}_tse.nii.gz

  # Result meshes
  TARGET=$WORK/shooting_target_native.vtk
  LM_PROCRUSTES=$WITER/shooting_target_procrustes.vtk
  LM_PROCRUSTES_MAT=$WITER/target_to_root_procrustes.mat
  SHOOTING_WARP=$WITER/shooting_warp.nii.gz

  # Mesh containing the momenta
  MOMENTA=$WITER/shooting_momenta.vtk

  # Target-related stuff in the WORK directory that is only done in the first iter
  if [[ $iter -eq 0 ]]; then
 
    # Get the warp chain from file
    WARPCHAIN=$(cat $MSTREGDIR/${ID}/final/chain_unwarp_to_final.txt) 

    # Apply the warp chain to the landmark mesh in template space, creating
    # the target locations for the geodesic shooting registration
    # -- this code works when the WARPCHAIN is empty (MST root)
    if [ ! -f $TARGET ]; then
    greedy -d 3 \
      -rf $REFSPACE \
      -rs $LANDMARKS $TARGET \
      -r $WARPCHAIN
    fi

  fi

  # Landmarks in reference space
  ln -sf $LANDMARKS $WITER/landmarks.vtk

  # Bring the target mesh back near the root mesh using procrustes alignment
  if [ ! -f $LM_PROCRUSTES_MAT ]; then
    vtkprocrustes $TARGET $LANDMARKS $LM_PROCRUSTES_MAT
  fi

  # Apply procrustes to the landmarks.
  #warpmesh $TARGET $LM_PROCRUSTES $LM_PROCRUSTES_MAT
  if [ ! -f $LM_PROCRUSTES ]; then
  greedy -d 3 \
    -rf $REFSPACE \
    -rs $TARGET $LM_PROCRUSTES \
    -r $LM_PROCRUSTES_MAT
  fi

  # Perform geodesic shooting between the procrustes landmarks and the
  # warped landmarks - this is going to allow us to interpolate the correspondence
  # found by the MST to the rest of the images
  if [[ ! -f $MOMENTA || ! -f $SHOOTING_WARP ]]; then

    time lmshoot -d 3 \
      -m $LANDMARKS $LM_PROCRUSTES \
      -s 2.0 -l 5000 -n 40 -i 200 0 -f \
      -o $MOMENTA \

    # Convert the shooting result into a warp
    lmtowarp -d 3 -n 40 -r $REFSPACE \
      -m $MOMENTA -o $SHOOTING_WARP \
      -s 2.0

  fi

  # Warp the native space image into the template
  if [ ! -f $WITER/reslice_${ID}_shooting_to_template_seg.nii.gz ]; then

  RM=""
  for sub in ${KINDS[*]}; do

    RM="$RM -rm $DATADIR/${ID}_${sub}.nii.gz $WITER/reslice_${ID}_shooting_to_template_${sub}.nii.gz"

  done

  greedy -d 3 \
    -rf $REFSPACE \
    $RM \
    -r $SHOOTING_WARP \
       $LM_PROCRUSTES_MAT,-1

  # Create seg
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $WITER/reslice_${ID}_shooting_to_template_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $WITER/reslice_${ID}_shooting_to_template_seg.nii.gz

  fi
}

function ShootingCorrectionShapeAvg_sub()
{
  # The iteration - 0 is the initial iteration, involves extra work
  iter=${1?}

  # The path to the landmarks
  SRC_LANDMARKS=${2?}

  # directory
  WORK=$SHOOTDIR/shape_avg/iter_${iter}
  mkdir -p $WORK

  # Reference space (root node in cm-rep space)
  REFSPACE=$MSTTEMPDIR/template_seg.nii.gz

  # The result landmarks - after shape averating
  SHAVG_LANDMARKS_NOPROC=$WORK/shavg_landmarks_noprocrustes.vtk
  SHAVG_LANDMARKS=$WORK/shavg_landmarks.vtk

  # Average the momentum maps from the previous iteration
  avgmesharr \
    $SHOOTDIR/*/iter_${iter}/shooting_momenta.vtk \
    InitialMomentum $SRC_LANDMARKS $WORK/average_momenta.vtk

  # Perform the shooting and generate warp
  if [ ! -f $WORK/average_momenta_warp.nii.gz ]; then
    lmtowarp \
      -d 3 -n 40 -r $REFSPACE \
      -m $WORK/average_momenta.vtk \
      -o $WORK/average_momenta_warp.nii.gz -s 2.0
  fi

  # Apply the warp to the landmarks to bring them to shape-averaging position
  greedy \
    -d 3 -rs $SRC_LANDMARKS $SHAVG_LANDMARKS_NOPROC \
    -rf $REFSPACE -r $WORK/average_momenta_warp.nii.gz

  # This transformation of the landmarks can cause shrinkage of the template. This
  # is not at all what we want in the template, we actually want the template to
  # keep its size during this iterative process. The way to correct this is to perform
  # procrustes between the source lanmarks and the new shape average
  vtkprocrustes $SRC_LANDMARKS $SHAVG_LANDMARKS_NOPROC $WORK/residual_procrustes.mat \
    | grep RMS_ | tee $WORK/procrustes_metric.txt

  # Applying the inverse of this procrustes to the SHAVG_LANDMARKS_NOPROC gives the new
  # template landmarks which are shape averaged but still the same size as the original
  # template.
  greedy \
    -d 3 -rs $SRC_LANDMARKS $SHAVG_LANDMARKS \
    -rf $REFSPACE \
    -r $WORK/average_momenta_warp.nii.gz $WORK/residual_procrustes.mat,-1
}

function ShootingCorrectionVolumeAvg_sub()
{
  # The iteration - 0 is the initial iteration, involves extra work
  iter=${1?}
  ITERDIR=$SHOOTDIR/template/iter_${iter}
  mkdir -p $ITERDIR

  # create mean image for different labels
  for sub in ${LABEL_IDS[*]}; do
    c3d $SHOOTDIR/*/iter_${iter}/reslice*${sub}.nii.gz \
      -mean -o $ITERDIR/template_gshoot_${sub}.nii.gz
  done

  # generate seg
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $ITERDIR/template_gshoot_${sub}.nii.gz; done) \
    -vote -type ushort \
    -o $ITERDIR/template_gshoot_seg.nii.gz

  # generate meshes
  for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

    LABELDEF=(${EXTRAMESHESDEF[i]})
    echo $LABELDEF

    c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$ITERDIR/template_gshoot_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
      -accum -add -endaccum \
      -o $ITERDIR/template_gshoot_${EXTRAMESHES[i]}.nii.gz

    vtklevelset \
      $ITERDIR/template_gshoot_${EXTRAMESHES[i]}.nii.gz \
      $ITERDIR/template_gshoot_${EXTRAMESHES[i]}.vtk 0.5

  done
}

##########################################################
function ShapeStats()
{
  # perform shape statistics at iteration ?
  iter=4

  PREFIX=SS${expid}
  # submit jobs
  #qsub -cwd -o $DUMPDIR -j y \
  #     -q all.q,basic.q \
  #     -l h_vmem=6.1G,s_vmem=6G \
  #     -N "${PREFIX}_${iter}" $0 \
       ShapeStats_sub $iter

  # Wait for completion
  #qsub -cwd -o $DUMPDIR -j y \
  #   -q all.q,basic.q \
  #   -hold_jid "${PREFIX}_*" -sync y -b y \
  #   sleep 1
}

function ShapeStats_sub()
{
  iter=$1
  ITERDIR=$SHOOTDIR/shape_avg/iter_${iter}
  PCADIR=$ITERDIR/pca
  rm -rf $PCADIR/*
  mkdir -p $PCADIR
  IDs=$(cat $ALLPATHDIR/IDSide.txt)

  # Generate the momentum data in a format that is readable by R script
  idx=1
  echo "ID" > $PCADIR/IDs.csv
  for id in $IDs; do
    local MOMENTS=$SHOOTDIR/$id/iter_${iter}/shooting_momenta.vtk
    echo $(echo $idx) $(dumpmeshattr $MOMENTS InitialMomentum) >> $PCADIR/initial_momenta.txt
    echo $idx >> $PCADIR/IDs.csv
    idx=$((idx+1))
  done

  # Generate an array of mesh coordinates for R
  dumpmeshpoints $ITERDIR/average_momenta.vtk | tail -n +3 > $PCADIR/template_points.txt

  # Run the R statistics notebook
  #Rscript run_pca_oneroot.R

}

##########################################################
function MakeMovies()
{
  # perform shape statistics at iteration ?
  iter=4

  PREFIX=MM${expid}
  for what in mode1 mode2 mode3 mode4 mode5; do

    # submit jobs
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=6.1G,s_vmem=6G \
         -N "${PREFIX}_${iter}_${what}" $0 \
         MakeMovies_sub $iter $what

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
     -q all.q,basic.q \
     -hold_jid "${PREFIX}_*" -sync y -b y \
     sleep 1  
}

function MakeMovies_sub()
{
  iter=$1
  what=$2
  ITERDIR=$SHOOTDIR/shape_avg/iter_${iter}
  VOLITERDIR=$SHOOTDIR/template/iter_${iter}
  PCADIR=$ITERDIR/pca
  MOVIEDIR=$PCADIR/movie_${what}
  mkdir -p $MOVIEDIR

  # Template
  local TEMPLATE=$VOLITERDIR/template_gshoot_seg.nii.gz
  local TEMPMESH=$VOLITERDIR/template_gshoot_MRG.vtk

  # Flow the mode vector forward and backward
  /data/picsl-build/pauly/cmrep/gcc64rel/lmtowarp \
    -m $PCADIR/${what}_vector_neg.vtk -s 2.0 -n 40 -d 3 -a 4 \
    -M $TEMPMESH $TMPDIR/neg_flow_%d.vtk
 
  /data/picsl-build/pauly/cmrep/gcc64rel/lmtowarp \
    -m $PCADIR/${what}_vector_pos.vtk -s 2.0 -n 40 -d 3 -a 4 \
    -M $TEMPMESH $TMPDIR/pos_flow_%d.vtk

  # Combine the flows into a single movie
  merge_neg_pos_movie $MOVIEDIR/movie_${what}_MRG_%02d.vtk $TEMPMESH
}

function merge_neg_pos_movie()
{
  local OUT_PATTERN=${1?}
  local CENTER=${2?}

  for ((i=4;i<=40;i+=4)); do
    local kneg=$((40-i))
    local kpos=$((40+i))
    cp $TMPDIR/neg_flow_${i}.vtk $(printf $OUT_PATTERN $kneg)
    cp $TMPDIR/pos_flow_${i}.vtk $(printf $OUT_PATTERN $kpos)
  done

  cp $CENTER $(printf $OUT_PATTERN 40)
}

##########################################################
function reset_dir()
{
  rm -rf $DUMPDIR/*
}

##########################################################
if [[ $# -lt 2 ]]; then

  main

else

  cmd=$1
  shift
  $cmd $@

fi


