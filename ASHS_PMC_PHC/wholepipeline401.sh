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

##############################################################################
# Parameters needs to be specify
# Experiment number
expid=401
expdir=$WORKDIR/exp/exp${expid}
DUMPDIR=${expdir}/dump
mkdir -p ${expdir}/dump

########################################
# 0. copy data
LABEL_IDS_ALL=(BKG      CA      DG  SUB ERC  BA35 BA36 PHC  CS   OTS)
LABEL_MRG_ALL=("0 7 16" "1 2 4" "3" "8" "10" "11" "12" "13" "14" "16")
# Labels to get segmentation
LABEL_IDS=(BKG      CA      DG  SUB ERC  BA35 BA36 PHC  CS)
LABEL_MRG=("0 7 16" "1 2 4" "3" "8" "10" "11" "12" "13" "14")
LABEL_NEW=(0        1       2   3   4    5    6    7    8)
KINDS="tse mprage ${LABEL_IDS[*]}"
DATADIR=${expdir}/data

########################################
# 1 pairwise registration and clustering
# 1.1 ANTs parameters
WGT_1=1
REG_LABELS_1=(${LABEL_IDS[*]})
ANTs_t_1="SyN[0.25]"
ANTs_r_1="Gauss[0.5,0]"
ANTs_i_1="15x8x0"
ANTs_x_1="Y"
ANTs_all_metrics_1="--use-all-metrics-for-convergence"
PWDIR=${expdir}/pairwise/

# 1.2 similarity type
#SIM_TYPE="PRCCS_seg_dice"
SIM1_TYPE="CS_HFdist"
SIM1_DIR=$PWDIR/sim_${SIM1_TYPE}
SIM2_TYPE="PRCCS_seg_dice"
SIM2_DIR=$PWDIR/sim_${SIM2_TYPE}

# 1.3 clustering
INITGRP_DIR=$expdir/clustering/InitGroup/

# 1.4 prepare statitics documents
groups=(1 2 3)
GPSTATDIR=$expdir/stat/groups/
ALLSTATDIR=$expdir/stat/all/

#########################################
# 2
# 2.1 assume only one parameter in each type of graph
GRAPHTYPES=(NNeighbor Epsilon)
GRAPHPARAMMIN=(0 -0.5)
IDS=($(cat $SUBJ_TXT))
N_all=${#IDS[*]}
N_max=$((N_all-1))
GRAPHPARAMMAX=(0 0.5)
GRAPHPARAMINC=(2 0.1)

# directories
GSEARCHDIR=$expdir/graphsearch
GSEARCHPARAMDIR=$GSEARCHDIR/params
GSEARCHPATHDIR=$GSEARCHDIR/paths
GSEARCHEVLDIR=$GSEARCHDIR/evaluation

# 2.2
FINALTEMPDIR=$expdir/FinalTemp
FINALGRAPHDIR=$FINALTEMPDIR/graph
FINALINITDIR=$FINALTEMPDIR/init
FINALWORKDIR=$FINALTEMPDIR/work

# 2.5
FINALITER=8
ANTs_start_final=3
WGT_3=1
REG_LABELS_3=(${LABEL_IDS[*]})
ANTs_t_3="SyN[0.25]"
ANTs_r_3="Gauss[0.5,0]"
ANTs_i_3="80x80x20"
ANTs_x_3="Y"
ANTs_all_metrics_3="--use-all-metrics-for-convergence"

# 2.6
# Relevant labels
LABEL_FG=(        CA  DG  SUB ERC BA35 BA36       PHC CS)
EXTRAMESHES=(HIPPO PRC ExtHippo MRG)
#                BKG CA DG SUB ERC BA35 BA36 PHC CS
EXTRAMESHESDEF=("-1   1  1  -1  -1  -1   -1   -1 -1" \
                "-1  -1 -1  -1  -1   1    1   -1 -1" \
                "-1  -1 -1   1   1   1    1    1 -1" \
                "-1   1 -1   1   1   1    1    1 -1")

# 2.7
FINALLABELDIR=$FINALTEMPDIR/labelwarp

# 2.8 warp meshes
FINALMESHDIR=$FINALTEMPDIR/meshwarp
subfield_th=0.5
template_th=0.0
thick_p=1.2
thick_e=6

# 2.9
EVALLABELS=(CA DG SUB ERC BA35 BA36 PRC   PHC CS HIPP     EXPHIPP   ALL)
RANGES=(    1  2  3   4   5    6    "5 6" 7   8  "1 2 3"  "4 5 6 7" "1 2 3 4 5 6 7")
MESH_EVAL=(       CA  DG  SUB ERC BA35 BA36 PRC   PHC CS)
MESH_EVAL_RANGES=(1   2   3   4   5    6    "5 6" 7   8)
FINALEVALDIR=$FINALTEMPDIR/evaluation

##############################################################################
function main()
{
  #reset_dir

  #######################################
  # 0. preparation
  #copy_data

  #######################################
  # 1 cluster subjects into three groups
  # 1.1 pairwise registration
  #pairwise
  #completeness

  # 1.2 compute similarity between subjects
  #similarity (need to be changed to have both!!)
  
  # 1.3 dimension reduction and perform clustering
  #clustering

  # 1.4 prepare the design and contrast documents
  #PrepStat

  ########################################
  # 2 build graph and generate super template
  # 2.1 compute paths with different parameters, apply warps, compute dice overlap
  #SearchGraphParam

  # 2.2 choose the best graph and the center subject (MATLAB)
  #ChooseBestGraph

  # 2.3 warp to center subject space 
  #WarpSuperTemplate

  # 2.4 average
  #AverageSuperTemplate

  # 2.5 register to the final super template
  #RegisterSuperTemplate

  # 2.6 make final images
  #MakeFinalImages

  # 2.7 warp super template labels
  WarpFinalLabel

  # 2.8 warp super template meshes
  WarpFinalMesh
  
  # 2.9 evaluate super template
  EvalFinalTemp

  # 2.10 displacement
  #DispFinalStat

  # 2.11 thickness analysis
  #ThickFinalStat


  ########################################
  # 5 statistical analyses for the whole group


}

#######################################################################
# Copy data
function copy_data()
{
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
           -l h_vmem=2.1G,s_vmem=2G \
           -N "${PREFIX}_${id}_${side}" \
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
           $DATADIR/${fn}_${side}_tse.nii.gz
    ln -sf $ASHSRUNDIR/$fn/mprage_to_chunktemp_${side}.nii.gz \
           $DATADIR/${fn}_${side}_mprage.nii.gz
  fi

  # Generate a binary image for the label
  for ((i=0;i<${#LABEL_IDS_ALL[*]};i++)); do

    c3d $SEG -replace $(for k in ${LABEL_MRG_ALL[i]}; do echo $k 999; done) -thresh 999 999 1 0 \
      -o $TMPDIR/binary_${LABEL_IDS_ALL[i]}_${fn}_${side}.nii.gz

    WarpImageMultiTransform 3  \
      $TMPDIR/binary_${LABEL_IDS_ALL[i]}_${fn}_${side}.nii.gz \
      $DATADIR/${fn}_${side}_${LABEL_IDS_ALL[i]}.nii.gz \
      -R $DATADIR/${fn}_${side}_tse.nii.gz \
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
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${id_fix}_${side}" $0 \
           pairwise_sub $id_fix $side
      sleep 0.1

    done

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function pairwise_sub()
{
  id_fix=$1
  side=$2

  IDS=$(cat $SUBJ_TXT)
  fn_fix=$(ls $ASHSRUNDIR | grep $id_fix)
  for id_mov in $IDS; do

    if [[ $id_mov != $id_fix ]]; then

      fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)
      OUTDIR=$PWDIR/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}
      mkdir -p $OUTDIR

      if [[ -f $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz ]]; then

        echo "Seg file exists."

      else

      # Use ml_affine for nice affine alignment
      /data/picsl/pauly/wolk/ashs/ext/Linux/bin/ml_affine \
        $DATADIR/${fn_fix}_${side}_seg.nii.gz \
        $DATADIR/${fn_mov}_${side}_seg.nii.gz \
        $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine.txt

      # Convert that to ITK format
      c3d_affine_tool \
        $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine.txt \
        -oitk $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine_itk.txt


      CMD=""
      for sub in ${REG_LABELS_1[*]}; do
        CMD="$CMD -m MSQ[$DATADIR/${fn_fix}_${side}_${sub}.nii.gz,$DATADIR/${fn_mov}_${side}_${sub}.nii.gz,$WGT_1]"
      done

      if [[ $ANTs_x_1 == "Y" ]]; then
        c3d $DATADIR/${fn_fix}_${side}_seg.nii.gz -dup \
          $DATADIR/${fn_mov}_${side}_seg.nii.gz \
          -int 0 -reslice-identity \
          -add -binarize -dilate 1 10x10x10vox \
          -o $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mask.nii.gz
        ANTs_mask_1="-x $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mask.nii.gz"
      else
        ANTs_mask_1=""
      fi

      # Perform ANTs registration
      ANTS 3 $CMD \
           -t $ANTs_t_1 \
           -r $ANTs_r_1 \
           -i $ANTs_i_1 \
           $ANTs_mask_1 \
           -a $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_mlaffine_itk.txt \
           --continue-affine 0 \
           $ANTs_all_metrics_1 \
           -o $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_pairwise.nii.gz \
           | tee $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_antsoutput.txt

      for sub in ${KINDS[*]}; do

        WarpImageMultiTransform 3 \
          $DATADIR/${fn_mov}_${side}_${sub}.nii.gz \
          $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_${sub}.nii.gz \
          -R $DATADIR/${fn_fix}_${side}_tse.nii.gz \
          $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_pairwiseWarp.nii.gz \
          $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_pairwiseAffine.txt

      done

      # Create seg
      c3d $(for sub in ${LABEL_IDS[*]}; do echo $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_${sub}.nii.gz; done) \
        -vote -type ushort \
        -o $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz

      fi
    fi

  done
}

##############################################################################
function completeness
{
  rm -rf $PWDIR/check_pairwise.txt

  for side in left right; do

    IDS=$(ls $DATADIR | grep ${side}_tse | sed -e "s/_${side}_tse.nii.gz//")

    for id_fix in $IDS; do

      for id_mov in $IDS; do

      if [[ $id_fix != $id_mov ]]; then

      # Check whether all the files exist
      OUTDIR=$PWDIR/${id_fix}_${side}/${id_mov}_to_${id_fix}
      if [[ -f $OUTDIR/${id_mov}_to_${id_fix}_${side}_reslice_seg.nii.gz ]]; then
          echo "ok"
      else
        echo "${id_mov} to ${id_fix} of ${side} is missing" \
          >> $PWDIR/check_pairwise.txt
      fi

      fi
      done
    done
  done

  if [ -f $PWDIR/check_pairwise.txt ]; then
    echo "pairwise registration has missing data"
    exit
  fi
}

##############################################################################
function similarity()
{
  rm -rf ${SIM_DIR}
  IDS=$(cat $SUBJ_TXT)

  PREFIX=SIM${expid}
  for side in left right; do

    mkdir -p ${SIM_DIR}/${side}

    for id_fix in $IDS; do

      qsub -cwd -o $DUMPDIR -j y \
           -l h_vmem=2.1G,s_vmem=2G \
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
  fn_out=$SIM_DIR/${side}/${id_fix}.txt
  fn_seg_fix=$TMPDIR/${fn_fix}_${side}_seg_tmp.nii.gz
  fn_seg_mesh_fix=$TMPDIR/${fn_fix}_${side}_seg_tmp.vtk
  OVL=""

  if [ ${SIM_TYPE} == "PRC_seg_dice" ]; then

    c3d $DATADIR/${fn_fix}_${side}_seg.nii.gz \
      -thresh 5 6 1 0 \
      -o $fn_seg_fix

  elif [ ${SIM_TYPE} == "CS_seg_dice" ]; then

    c3d $DATADIR/${fn_fix}_${side}_seg.nii.gz \
      -thresh 8 8 1 0 \
      -o $fn_seg_fix

  elif [ ${SIM_TYPE} == "PRCCS_seg_dice" ]; then

    c3d $DATADIR/${fn_fix}_${side}_seg.nii.gz \
      -replace 1 0 2 0 3 0 4 0 7 0 \
      -o $fn_seg_fix

  elif [ ${SIM_TYPE} = "CS_HFdist" ]; then

    c3d $DATADIR/${fn_fix}_${side}_seg.nii.gz \
      -replace 8 99 -thresh 99 99 1 0 \
      -o $fn_seg_fix
    vtklevelset $fn_seg_fix $fn_seg_mesh_fix 0.5

  fi

  # Go through all the other subjects
  IDS=$(cat $SUBJ_TXT)
  for id_mov in $IDS; do

    fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)
    OVL_tmp="NA"

    if [[ $id_mov == $id_fix ]]; then

      if [ ${SIM_TYPE} = "PRC_HF_dist" ]; then

        OVL="$OVL 0"

      else

        OVL="$OVL 1"

      fi

    else

      if [ ${SIM_TYPE} == "PRC_seg_dice" ]; then

        OVL_tmp=$(c3d $PWDIR/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz \
          -thresh 5 6 1 0 \
          $fn_seg_fix -overlap 1 \
          | grep OVL | awk -F '[ ,]+' '{print $6}' )

      elif [ ${SIM_TYPE} == "CS_seg_dice" ]; then

        OVL_tmp=$(c3d $PWDIR/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz \
          -thresh 8 8 1 0 \
          $fn_seg_fix -overlap 1 \
          | grep OVL | awk -F '[ ,]+' '{print $6}' )

      elif [ ${SIM_TYPE} == "PRCCS_seg_dice" ]; then

        OVL_tmp=$(c3d $PWDIR/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz \
          -replace 1 0 2 0 3 0 4 0 7 0 \
          $fn_seg_fix -label-overlap \
          | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}' )

      elif [ ${SIM_TYPE} = "CS_HFdist" ]; then

        c3d $PWDIR/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz \
          -replace 8 99 -thresh 99 99 1 0 \
          -o $TMPDIR/${fn_mov}_to_${fn_fix}_seg_tmp.nii.gz
        vtklevelset $TMPDIR/${fn_mov}_to_${fn_fix}_seg_tmp.nii.gz \
          $TMPDIR/${fn_mov}_to_${fn_fix}_seg_tmp.vtk 0.5

        # meshdiff
        OVL_tmp1=($(meshdiff $fn_seg_mesh_fix $TMPDIR/${fn_mov}_to_${fn_fix}_seg_tmp.vtk | grep RESULT | awk '{print $9} {print $10}'))
        
        if [[ $OVL_tmp1 == "" ]]; then
          OVL_tmp='NA'
        else
          st=$(echo "${OVL_tmp1[0]} > ${OVL_tmp1[1]}" | bc)
          if [ $st -gt 0 ]; then
            OVL_tmp=${OVL_tmp1[0]}
          else
            OVL_tmp=${OVL_tmp1[1]}
          fi
        fi

      fi

      OVL="$OVL $OVL_tmp"

    fi
  done

  echo $OVL > ${fn_out}

}

##############################################################################
function clustering()
{
  mkdir -p $INITGRP_DIR

  PREFIX=CLU${expid}
  for side in left right; do

      qsubp5 -cwd -o $DUMPDIR -j y \
           -l h_vmem=8.1G,s_vmem=8G \
           -N "${PREFIX}_${side}" $0 \
           clustering_sub $side
      sleep 0.1

  done

  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function clustering_sub()
{
  side=$1

  # Run matlab to denoise T1 image
  mkdir -p $TMPDIR/${id}
  source $PKGDIR/matlab_batch.sh \
     $TMPDIR/${id}/ \
     $MATLABCODEDIR \
     clustering_ProbPCA_2sim \
     ${SIM1_DIR}/sim_${SIM1_TYPE}_${side}.txt \
     ${SIM2_DIR}/sim_${SIM2_TYPE}_${side}.txt \
     $GTGROUPDIR/group_${side}_200.txt \
     $SUBJ_TXT \
     $side \
     $INITGRP_DIR
}

#############################################################################
function PrepStat()
{
  rm -rf $GPSTATDIR $ALLSTATDIR
  mkdir -p $GPSTATDIR $ALLSTATDIR
  
  # get IDs
  IDS=($(cat $SUBJ_TXT))
  MCIGRP=($(cat $WORKDIR/analysis_input/CNMCI.txt))  
  AGE=($(cat $WORKDIR/analysis_input/age.txt))
  ICV=($(cat $WORKDIR/analysis_input/icv.txt))

  for side in left right; do

    # group info
    GRP=($(cat $INITGRP_DIR/group_${side}.txt))
    
    for ((i=0;i<${#GRP[*]};i++)); do

      id=${IDS[i]}
      grp=${GRP[i]}
      age=${AGE[i]}
      icv=${ICV[i]}
      mcigrp=${MCIGRP[i]}

      if [[ $mcigrp -eq 1 ]]; then
        str="$id	1	0	$age	$icv"
      else
        str="$id	0       1       $age    $icv"
      fi
      
      echo $str | sed "s/\r//g" >> \
        $GPSTATDIR/design_${side}_${grp}_group-age-icvcr.txt
      echo $str | sed "s/\r//g" >> \
        $ALLSTATDIR/design_${side}_group-age-icvcr.txt

    done

  done

  echo "1 -1 0 0" > \
    $GPSTATDIR/contrast_group-age-icvcr_nc-mci.txt
  echo "-1 1 0 0" > \
    $GPSTATDIR/contrast_group-age-icvcr_mci-nc.txt
  echo "1 -1 0 0" > \
    $ALLSTATDIR/contrast_group-age-icvcr_nc-mci.txt
  echo "-1 1 0 0" > \
    $ALLSTATDIR/contrast_group-age-icvcr_mci-nc.txt
}

##############################################################################
function SearchGraphParam()
{
  # make directory
  mkdir -p $GSEARCHPARAMDIR
  rm -rf $GSEARCHPARAMDIR/*

  # IDS
  IDS=($(cat $SUBJ_TXT))

  # for each graph type
  for ((i=0;i<${#GRAPHTYPES[*]};i++)); do

    graphtype=${GRAPHTYPES[i]}
    parammin=${GRAPHPARAMMIN[i]}
    parammax=${GRAPHPARAMMAX[i]}
    paraminc=${GRAPHPARAMINC[i]}

    # for each parameter
    count=0
    while true; do

      # evaluation tem
      GSEARCHEVLTMPDIR=$GSEARCHEVLDIR/tmp
      rm -rf $GSEARCHEVLTMPDIR $GSEARCHPATHDIR
      mkdir -p $GSEARCHEVLTMPDIR $GSEARCHPATHDIR

      # compute parameter
      param=$(echo "$parammin + $paraminc * $count" | bc)

      # check whether parameter is bigger than maximum
      st=$(echo "$param > $parammax" | bc)
      if [[ $st -eq 1 ]]; then
        break
      else
        count=$((count+1))
      fi

      # save the searching parameter
      echo "$param" >> $GSEARCHPARAMDIR/${graphtype}.txt

      ###########################################
      # get paths using MATLAB script
      PREFIX="GP${expid}"
      for side in left right; do

        # if exist continue
        if [ -f $GSEARCHEVLDIR/${graphtype}_${param}_${side}.txt ]; then
          continue
        fi

        qsub -cwd -o $DUMPDIR -j y -N \
           "${PREFIX}_${side}_${graphtype}_${param}" $0 \
           GraphPaths_sub $side $graphtype $param
      done

      # Wait for completion
      qsub -cwd -o $DUMPDIR -j y \
        -hold_jid "${PREFIX}_*" -sync y -b y \
        sleep 1

      ###########################################
      # compute pairwise similarity
      PREFIX="PS${expid}"
      for side in left right; do

        # if exist continue
        if [ -f $GSEARCHEVLDIR/${graphtype}_${param}_${side}.txt ]; then
          continue
        fi

        for ((j=0;j<${#IDS[*]};j++)); do
          qsub -cwd -o $DUMPDIR -j y -N \
            "${PREFIX}_${j}_${side}" $0 \
            PWSim_sub $side $j $GSEARCHEVLTMPDIR

        done
      done

      # Wait for completion
      qsub -cwd -o $DUMPDIR -j y \
        -hold_jid "${PREFIX}_*" -sync y -b y \
        sleep 1

      ############################################
      # combine the evaluation result
      for side in left right; do

        # if exist continue
        if [ -f $GSEARCHEVLDIR/${graphtype}_${param}_${side}.txt ]; then
          continue
        fi        

        for id_fix in ${IDS[*]}; do
          cat $GSEARCHEVLTMPDIR/${id_fix}_${side}.txt \
              >> $GSEARCHEVLDIR/${graphtype}_${param}_${side}.txt
        done
      done

    done
  done
}

function GraphPaths_sub()
{
  side=$1
  graphtype=$2
  param=$3

  # Run matlab to denoise T1 image
  mkdir -p $TMPDIR/${side}_${graphtype}_${param}
  source $PKGDIR/matlab_batch.sh \
     $TMPDIR/${side}_${graphtype}_${param} \
     $MATLABCODEDIR \
     construct_one_graph \
     $INITGRP_DIR/allinfo_${side}.mat \
     $GSEARCHPATHDIR \
     $side \
     $graphtype \
     $param
}

function PWSim_sub()
{
  side=$1
  idx_fix=$2
  EVLTMPDIR=$3
  IDS=($(cat $SUBJ_TXT))
  id_fix=${IDS[$idx_fix]}
  fn_fix=$(ls $ASHSRUNDIR | grep $id_fix)
  OVL=""

  c3d $DATADIR/${fn_fix}_${side}_seg.nii.gz \
      -replace 1 0 2 0 3 0 4 0 7 0 \
      -o $TMPDIR/${fn_fix}_${side}_seg_tmp.nii.gz

  for id_mov in ${IDS[*]}; do

    fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)

    if [[ $id_mov == $id_fix ]]; then
      OVL_tmp="1"
    else

      # get the path
      paths=($(cat $GSEARCHPATHDIR/${id_mov}_${side}.txt))
      cur_path=${paths[$((idx_fix))]}

      # loop
      idx=1
      warps=""
      while true; do

        set +e
        # from subject
        idx_from=$(echo $cur_path | cut -d , -f $idx)
        id_from=${IDS[$((idx_from-1))]}
        id_from=$(ls $ASHSRUNDIR | grep $id_from)

        # to subject
        idx=$((idx+1))
        idx_to=$(echo $cur_path | cut -d , -f $idx)
        id_to=${IDS[$((idx_to-1))]}
        id_to=$(ls $ASHSRUNDIR | grep $id_to)
        set -e

        if [[ $idx_to == "" ]] || [[ $idx_from == $idx_to ]]; then
          break
        fi

        # get warps
        warps="$PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_pairwiseWarp.nii.gz $PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_mlaffine_itk.txt $warps"

      done

      $ANTSPATH/ComposeMultiTransform 3 \
        $TMPDIR/${fn_mov}_to_${fn_fix}_${side}_Warp.nii.gz \
        -R $DATADIR/${fn_fix}_${side}_tse.nii.gz \
        $warps

      # apply init warps
      for sub in $KINDS; do
        WarpImageMultiTransform 3 \
          $DATADIR/${fn_mov}_${side}_${sub}.nii.gz \
          $TMPDIR/${fn_mov}_to_${fn_fix}_${side}_${sub}.nii.gz \
          -R $DATADIR/${fn_fix}_${side}_tse.nii.gz \
          $TMPDIR/${fn_mov}_to_${fn_fix}_${side}_Warp.nii.gz
      done

      # Create the segmentation for the template
      c3d $(for sub in ${LABEL_IDS[*]}; do echo $TMPDIR/${fn_mov}_to_${fn_fix}_${side}_${sub}.nii.gz; done) \
        -vote -type ushort \
        -o $TMPDIR/${fn_mov}_to_${fn_fix}_${side}_seg.nii.gz

      # compute overlap
      OVL_tmp=$(c3d $TMPDIR/${fn_mov}_to_${fn_fix}_${side}_seg.nii.gz \
          -replace 1 0 2 0 3 0 4 0 7 0 \
          $TMPDIR/${fn_fix}_${side}_seg_tmp.nii.gz -label-overlap \
          | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}' )
    fi
    OVL="$OVL $OVL_tmp"

  done

  echo $OVL > $EVLTMPDIR/${id_fix}_${side}.txt
}

##############################################################################
function ChooseBestGraph()
{
  mkdir -p $INITGRP_DIR

  PREFIX=CBG${expid}
  for side in left right; do

      qsub -cwd -o $DUMPDIR -j y \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${side}" $0 \
           ChooseBestGraph_sub $side
      sleep 0.1

  done

  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function ChooseBestGraph_sub()
{
  side=$1

  # Run matlab to denoise T1 image
  mkdir -p $TMPDIR/${side}
  source $PKGDIR/matlab_batch.sh \
     $TMPDIR/${side}/ \
     $MATLABCODEDIR \
     construct_one_graph_2sim \
     $INITGRP_DIR/allinfo_${side}.mat \
     $FINALGRAPHDIR/paths \
     $side \
     Epsilon \
     3 \
     combine \
     
  if [[ $side == "left" ]]; then
    echo "65" > $FINALGRAPHDIR/centerIdx_${side}.txt
  else
    echo "4" > $FINALGRAPHDIR/centerIdx_${side}.txt
  fi

     #choose_graph \
     #$INITGRP_DIR/allinfo_${side}.mat \
     #$GSEARCHPARAMDIR \
     #$GSEARCHEVLDIR \
     #$side \
     #$FINALGRAPHDIR
}

##############################################################################
function WarpSuperTemplate()
{
  # IDS
  IDS=($(cat $SUBJ_TXT))
  mkdir -p $FINALINITDIR

  PREFIX=WST${expid}
  for side in left right; do

    # load the fix idx
    idx_fix=$(cat $FINALGRAPHDIR/centerIdx_${side}.txt)

    # loop through each subject
    for ((i=0;i<${#IDS[*]};i++)); do

      id_mov=${IDS[i]}

      # submit jobs
      qsub -cwd -o $DUMPDIR -j y \
        -N "${PREFIX}_${side}_${id_mov}_${idx_fix}" $0 \
        WarpSuperTemplate_sub $side $id_mov $idx_fix
      sleep 0.1

    done
  done

  # wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function WarpSuperTemplate_sub()
{
  side=$1
  id_mov=$2
  idx_fix=$3

  # IDS
  IDS=($(cat $SUBJ_TXT))
  id_fix=${IDS[$((idx_fix-1))]}
  fn_fix=$(ls $ASHSRUNDIR | grep $id_fix)
  fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)

  if [[ $id_fix == $id_mov ]]; then

    for sub in $KINDS; do

      cp $DATADIR/${fn_mov}_${side}_${sub}.nii.gz \
        $FINALINITDIR/${fn_mov}_${side}_totemp_reslice_${sub}.nii.gz

    done

  else

    # get the path
    paths=($(cat $FINALGRAPHDIR/paths/${id_mov}_${side}.txt))
    cur_path=${paths[$((idx_fix-1))]}

    # loop
    idx=1
    warps=""
    invwarps=""
    while true; do

      set +e
      # from subject
      idx_from=$(echo $cur_path | cut -d , -f $idx)
      id_from=${IDS[$((idx_from-1))]}
      id_from=$(ls $ASHSRUNDIR | grep $id_from)

      # to subject
      idx=$((idx+1))
      idx_to=$(echo $cur_path | cut -d , -f $idx)
      id_to=${IDS[$((idx_to-1))]}
      id_to=$(ls $ASHSRUNDIR | grep $id_to)
      set -e

      if [[ $idx_to == "" ]] || [[ $idx_from == $idx_to ]]; then
        break
      fi

      # get warps
      warps="$PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_pairwiseWarp.nii.gz $PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_mlaffine_itk.txt $warps"

      invwarps="$invwarps -i $PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_mlaffine_itk.txt $PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_pairwiseInverseWarp.nii.gz"

    done

    # compose warps
    $ANTSPATH/ComposeMultiTransform 3 \
       $FINALINITDIR/${fn_mov}_${side}_totempWarp.nii.gz \
       -R $DATADIR/${fn_fix}_${side}_tse.nii.gz \
       $warps
    $ANTSPATH/ComposeMultiTransform 3 \
       $FINALINITDIR/${fn_mov}_${side}_totempInverseWarp.nii.gz \
       -R $DATADIR/${fn_mov}_${side}_tse.nii.gz \
       $invwarps

    for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        $DATADIR/${fn_mov}_${side}_${sub}.nii.gz \
        $FINALINITDIR/${fn_mov}_${side}_totemp_reslice_${sub}.nii.gz \
        -R $DATADIR/${fn_fix}_${side}_tse.nii.gz \
        $FINALINITDIR/${fn_mov}_${side}_totempWarp.nii.gz

    done

  fi

  # Create the segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALINITDIR/${fn_mov}_${side}_totemp_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $FINALINITDIR/${fn_mov}_${side}_totemp_reslice_seg.nii.gz
}

##############################################################################
function AverageSuperTemplate()
{
  mkdir -p $FINALWORKDIR

  # generate the average template
  PREFIX=AST${expid}  
  for side in left right; do
    # Compute average images
    for kind in $KINDS; do

      qsub -cwd -o $DUMPDIR -j y -N \
        "${PREFIX}_${kind}_${side}" $0 \
        AverageSuperTemplate_sub $side $kind

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
     -hold_jid "${PREFIX}_*" -sync y -b y \
     sleep 1

  for side in left right; do

    # Create the segmentation for the template
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALWORKDIR/template_${side}_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $FINALWORKDIR/template_${side}_seg.nii.gz

    # generate mask
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALWORKDIR/template_${side}_${sub}.nii.gz; done | grep -v BKG) \
      -mean -thresh 1e-5 inf 1 0 \
      -o $FINALWORKDIR/template_${side}_mask.nii.gz

  done
}

function AverageSuperTemplate_sub()
{
  side=$1
  kind=$2
  
  if [[ $kind = "tse" ]]; then NORM=1; else NORM=0; fi
    AverageImages 3 \
      $FINALWORKDIR/template_${side}_${kind}.nii.gz \
      $NORM \
      $FINALINITDIR/*_${side}_totemp_reslice_${kind}.nii.gz
}

##############################################################################
function RegisterSuperTemplate()
{
  # choose the best template to be the super template
  PREFIX=RST${expid}
  for side in left right; do

    for ((iter=0;iter<$FINALITER;iter++)); do

      # Back up template
      ITDIR=$FINALWORKDIR/$(printf iter_%s_%02d $side $iter)
      mkdir -p $ITDIR
      cp -a $FINALWORKDIR/template_${side}_*.nii.gz $ITDIR/

      # Do ants?
      if [[ $iter -lt $ANTs_start_final ]]; then doants=0; else doants=1; fi

      IDS=$(cat $SUBJ_TXT)

      # Run ANTS for each image
      for id in $IDS; do

        id=$(ls $ASHSRUNDIR | grep $id)

        #Submit ANTS job
        qsub -cwd -o $DUMPDIR -j y -N \
           "${PREFIX}_${id}_${side}" $0 \
           RegisterSuperTemplate_sub $id $side $doants
        sleep 0.1

      done

      # Wait for completion
      qsub -cwd -o $DUMPDIR -j y \
        -hold_jid "${PREFIX}_*" -sync y -b y \
        sleep 1

      # If this is the last iteration, we don't want to recompute the template
      if [[ $iter -lt $((FINALITER-1)) ]]; then

        for grp in ${groups[*]}; do

          # Compute average images
          for kind in $KINDS; do

            if [[ $kind = "tse" ]]; then NORM=1; else NORM=0; fi
            AverageImages 3 $FINALWORKDIR/template_${side}_${kind}.nii.gz $NORM $FINALWORKDIR/*_${side}_totemp_reslice_${kind}.nii.gz

          done

          # Create the segmentation for the template
          c3d \
            $(for sub in ${LABEL_IDS[*]}; do \
              echo $FINALWORKDIR/template_${side}_${sub}.nii.gz; \
            done) \
            -vote -type ushort \
            -o $FINALWORKDIR/template_${side}_seg.nii.gz

          # Perform shape averaging
          if [[ $doants -eq 1 ]]; then
            shape_update_to_template $side $grp
          fi

        done

      fi

    done
  done
}

function RegisterSuperTemplate_sub()
{
  id=$1
  side=$2
  doants=$3

  if [[ ! -f $FINALINITDIR/${id}_${side}_totemp_reslice_seg.nii.gz ]]; then
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALINITDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $FINALINITDIR/${id}_${side}_totemp_reslice_seg.nii.gz
  fi

  # Before we vote, use ml_affine for nice affine alignment
  ~pauly/wolk/ashs/ext/Linux/bin/ml_affine \
    $FINALWORKDIR/template_${side}_seg.nii.gz \
    $FINALINITDIR/${id}_${side}_totemp_reslice_seg.nii.gz \
    $FINALWORKDIR/${id}_${side}_mlaffine.txt

  # Do we do ANTS or not?
  if [[ $doants -eq 0 ]]; then

    for sub in $KINDS; do

      c3d \
        $FINALWORKDIR/template_${side}_tse.nii.gz \
        $FINALINITDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz \
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
      CMD="$CMD -m MSQ[$FINALWORKDIR/template_${side}_${sub}.nii.gz,$FINALINITDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz,$WGT_3]"
    done

    if [[ $ANTs_x_3 == "Y" ]]; then
       ANTs_mask_3="-x $FINALWORKDIR/template_${side}_mask.nii.gz"
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
      -o $FINALWORKDIR/${id}_${side}_totemp.nii.gz \
      | tee $FINALWORKDIR/${id}_${side}_antsoutput.txt

    for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        $FINALINITDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz \
        $FINALWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz \
        -R $FINALWORKDIR/template_${side}_tse.nii.gz \
        $FINALWORKDIR/${id}_${side}_totempWarp.nii.gz \
        $FINALWORKDIR/${id}_${side}_totempAffine.txt

    done

    c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $FINALWORKDIR/${id}_${side}_totemp_reslice_seg.nii.gz

  fi
}

function shape_update_to_template()
{
  side=$1
  grp=$2

  # Borrowed from ANTS buildtemplateparallel

  # Average the warp fields
  local TEMPWARP=$FINALWORKDIR/template_${side}_warp.nii.gz
  $NEWANTSDIR/AverageImages 3 $TEMPWARP 0 \
    $FINALWORKDIR/*_${side}_totempWarp.nii.gz

  # Scale by -0.25 (gradient step)
  MultiplyImages 3 $TEMPWARP -0.25 $TEMPWARP

  # Create a new template affine
  local TEMPAFF=$FINALWORKDIR/template_${side}_Affine.txt
  if [[ -f $TEMPAFF ]]; then rm -rf $TEMPAFF; fi

  cat $FINALWORKDIR/*_${side}_totempAffine.txt \
    | grep '^Parameters:' | awk '\
    BEGIN { for(i=0;i<12;i++) x[i]=0} \
    { for(i=0;i<12;i++) x[i]=x[i]+$(i+2) } \
    END { \
      printf "Transform: MatrixOffsetTransformBase_double_3_3\nParameters: "; \
      for(i=0;i<12;i++) printf "%f ",x[i]/NR; \
      printf "\nFixedParameters: 0 0 0\n";}' > $TEMPAFF

  # Compose the warps
  WarpImageMultiTransform 3 $TEMPWARP $TEMPWARP -i $TEMPAFF \
    -R $FINALWORKDIR/template_${side}_tse.nii.gz

  TEMPWARPFULL=$FINALWORKDIR/template_${side}_fullwarp.nii.gz
  $ANTSPATH/ComposeMultiTransform 3 \
    $TEMPWARPFULL -R $FINALWORKDIR/template_${side}_tse.nii.gz \
    -i $TEMPAFF $TEMPWARP $TEMPWARP $TEMPWARP $TEMPWARP

  # Apply this warp to all the template derivatives
  for kind in $KINDS; do

    WarpImageMultiTransform 3 \
      $FINALWORKDIR/template_${side}_${kind}.nii.gz \
      $FINALWORKDIR/template_${side}_${kind}.nii.gz \
      $TEMPWARPFULL

  done
}

##############################################################################
function MakeFinalImages()
{
  for side in left right; do

    # Generate one compound label for all other subfields
    ITDIR=$FINALWORKDIR/$(printf iter_%s_%02d $side $((FINALITER-1)) )

    for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

      LABELDEF=(${EXTRAMESHESDEF[i]})
      echo $LABELDEF

      c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$ITDIR/template_${side}_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
        -accum -add -endaccum \
        -o $ITDIR/template_${side}_${EXTRAMESHES[i]}.nii.gz

    done
  done
}

##############################################################################
function WarpFinalLabel()
{
  mkdir -p $FINALLABELDIR

  PREFIX=WFL${expid}
  # Iterate over side and subject
  for side in left right; do

    IDS=$(cat $SUBJ_TXT)

    for id in $IDS; do

      id=$(ls $ASHSRUNDIR | grep $id)

      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
        -N "${PREFIX}_${id}_${side}" \
        $0 WarpFinalLabel_sub $id $side

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
    -hold_jid "${PREFIX}_*" -sync y -b y \
    sleep 1
}

function WarpFinalLabel_sub()
{
  id=$1
  side=$2

  # check whether it is the center subject
  IDS=($(cat $SUBJ_TXT))
  center_idx=$(cat $FINALGRAPHDIR/centerIdx_${side}.txt)
  center_id=${IDS[$((center_idx-1))]}
  center_id=$(ls $ASHSRUNDIR | grep $center_id)

  # Apply the warp to each of the labels
  for sub in ${LABEL_IDS[*]}; do

    if [[ $center_id == $id ]]; then

    antsApplyTransforms -d 3 \
      -r $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray.nii.gz \
      -i $FINALWORKDIR/template_${side}_${sub}.nii.gz \
      -o $FINALLABELDIR/${id}_${side}_${sub}_tempfit.nii.gz \
      -t [$ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt,1] \
      -t [$ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt,1] \
      -t $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii \
      -t [$FINALWORKDIR/${id}_${side}_totempAffine.txt,1] \
      -t $FINALWORKDIR/${id}_${side}_totempInverseWarp.nii.gz \
      -n BSpline

    else

    antsApplyTransforms -d 3 \
      -r $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray.nii.gz \
      -i $FINALWORKDIR/template_${side}_${sub}.nii.gz \
      -o $FINALLABELDIR/${id}_${side}_${sub}_tempfit.nii.gz \
      -t [$ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt,1] \
      -t [$ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt,1] \
      -t $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii \
      -t $FINALINITDIR/${id}_${side}_totempInverseWarp.nii.gz \
      -t [$FINALWORKDIR/${id}_${side}_totempAffine.txt,1] \
      -t $FINALWORKDIR/${id}_${side}_totempInverseWarp.nii.gz \
      -n BSpline

    fi

  done

  # Vote in original subject space to get the tempfit segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALLABELDIR/${id}_${side}_${sub}_tempfit.nii.gz; done) -vote -type ushort \
    -o $FINALLABELDIR/${id}_${side}_seg_tempfit.nii.gz


  # change ASHS labels to the current formate
  for ((i=0;i<${#LABEL_IDS[*]};i++)); do
    c3d $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray.nii.gz \
      -replace $(for k in ${LABEL_MRG[i]}; do echo $k 999; done) \
      -thresh 999 999 1 0 \
      -o $FINALLABELDIR/${id}_${side}_${LABEL_IDS[i]}_ashs.nii.gz
  done

  c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALLABELDIR/${id}_${side}_${sub}_ashs.nii.gz; done) -vote -type ushort \
    -o $FINALLABELDIR/${id}_${side}_seg_ashs.nii.gz
}

##############################################################################
function WarpFinalMesh()
{
  PREFIX=WFM${expid}

  mkdir -p $FINALMESHDIR
  mkdir -p $FINALTEMPDIR/jacobian
  mkdir -p $FINALMESHDIR/template

  # Iterate over side and subject
  for side in left right; do

    ITDIR=$FINALWORKDIR/$(printf iter_%s_%02d $side $((FINALITER-1)) )

    # Generate meshes for the individual subfields
    for sub in ${LABEL_FG[*]}; do

      vtklevelset \
        $ITDIR/template_${side}_${sub}.nii.gz \
        $FINALMESHDIR/template_${side}_${sub}.vtk \
        $subfield_th

    done

    # Generate one mesh for all non-DG non-CS subfields
    for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

      vtklevelset \
        $ITDIR/template_${side}_${EXTRAMESHES[i]}.nii.gz \
        $FINALMESHDIR/template_${side}_${EXTRAMESHES[i]}.vtk \
        $template_th

    done

    IDS=$(cat $SUBJ_TXT)

    for id in $IDS; do

      id=$(ls $ASHSRUNDIR | grep $id)

      # Submit job for this subject
      qsub  -cwd -o $DUMPDIR -j y \
        -N "${PREFIX}_${id}_${side}" \
        $0 WarpFinalMesh_sub $id $side

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
     -hold_jid "${PREFIX}_*" -sync y -b y \
     sleep 1

  cp $FINALMESHDIR/template*.vtk \
     $FINALMESHDIR/template/
}

function WarpFinalMesh_sub()
{
  id=$1
  side=$2
  ALLSF="${LABEL_FG[*]} ${EXTRAMESHES[*]}"

  # check whether it is the center subject
  IDS=($(cat $SUBJ_TXT))
  center_idx=$(cat $FINALGRAPHDIR/centerIdx_${side}.txt)
  center_id=${IDS[$((center_idx-1))]}
  center_id=$(ls $ASHSRUNDIR | grep $center_id)

  if [[ $center_id == $id ]]; then
  # Compose the transformation between the template and the subject
  $ANTSPATH/ComposeMultiTransform 3 \
     $TMPDIR/compose.nii \
     -R $FINALWORKDIR/template_${side}_tse.nii.gz \
     $FINALWORKDIR/${id}_${side}_totempWarp.nii.gz \
     $FINALWORKDIR/${id}_${side}_totempAffine.txt \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempWarp.nii \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
     $ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt
  else
  # Compose the transformation between the template and the subject
  $ANTSPATH/ComposeMultiTransform 3 \
     $TMPDIR/compose.nii \
     -R $FINALWORKDIR/template_${side}_tse.nii.gz \
     $FINALWORKDIR/${id}_${side}_totempWarp.nii.gz \
     $FINALWORKDIR/${id}_${side}_totempAffine.txt \
     $FINALINITDIR/${id}_${side}_totempWarp.nii.gz \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempWarp.nii \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
     $ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt
  fi

  # Split the transformations into xyz
  c3d -mcs $TMPDIR/compose.nii -oo $TMPDIR/comp%02d.nii

  # Apply the warp to each of the meshes
  for sub in $ALLSF; do

    # Warp the subfield into subject space
    warpmesh -w ants -m ras \
      $FINALMESHDIR/template_${side}_${sub}.vtk \
      $FINALMESHDIR/${id}_${side}_${sub}_tempfit.vtk \
      $TMPDIR/comp00.nii $TMPDIR/comp01.nii $TMPDIR/comp02.nii

    # Extract the thickness of the subfield
    cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
      -T $FINALMESHDIR/${id}_${side}_${sub}_thickmap.vtk \
      -p $thick_p -e $thick_e \
      $FINALMESHDIR/${id}_${side}_${sub}_tempfit.vtk \
      $FINALMESHDIR/skel_${id}_${side}_${sub}.vtk

  done

  # Compute the Jacobian map
  ANTSJacobian 3 $TMPDIR/compose.nii \
    $FINALTEMPDIR/jacobian/${id}_${side}_totemp 1
}

##############################################################################
function EvalFinalTemp()
{
  FINALEVALTMPDIR=$FINALEVALDIR/tmp
  rm -rf $FINALEVALTMPDIR
  mkdir -p $FINALEVALTMPDIR
  IDS=$(cat $SUBJ_TXT)

  PREFIX=EFT${expid}
  # Iterate over side and subject
  for side in left right; do
    for id in $IDS; do

      id=$(ls $ASHSRUNDIR | grep $id)

      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
        -l h_vmem=4.1G,s_vmem=4G \
        -N "${PREFIX}_${id}_${side}" \
        $0 EvalFinalTemp_sub $id $side $FINALEVALTMPDIR
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # Combine scores
  for side in left right; do

    echo "ID ${EVALLABELS[*]} ${MESH_EVAL[*]} ${EVALLABELS[*]}" \
      > $FINALEVALDIR/overlap_${side}.txt

    for id in $IDS; do

      id=$(ls $ASHSRUNDIR | grep $id)
      cat $FINALEVALTMPDIR/${id}_${side}_overlap.txt \
          >> $FINALEVALDIR/overlap_${side}.txt

    done
  done
}

function EvalFinalTemp_sub()
{
  id=$1
  side=$2
  EVALTMPDIR=$3
  ALLOVL=""

  ###################################
  # Compute dice overlap in subject space
  FITTED=$FINALLABELDIR/${id}_${side}_seg_tempfit.nii.gz
  ASHSSEG=$FINALLABELDIR/${id}_${side}_seg_ashs.nii.gz
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

      MESH_DIST_TMP="$(meshdiff $FINALMESHDIR/${id}_${side}_${MESH_EVAL[i]}_tempfit.vtk $TMPDIR/mask_seg.vtk | grep RESULT | awk '{print $7}')"

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
  do_pair $FINALWORKDIR/template_${side}_seg.nii.gz \
    $FINALWORKDIR/${id}_${side}_totemp_reslice_seg.nii.gz
  ALLOVL="$ALLOVL $FULLOVL"

  ###################################
  # save to temporary file
  echo $id $ALLOVL > $EVALTMPDIR/${id}_${side}_overlap.txt
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

  # Save this temporary dice file
  #echo $id $FULLOVL > $out_dice_file
}

##############################################################################
function DispFinalTemp()
{
  #rm -rf $FINALMESHDIR/analysis
  mkdir -p $FINALMESHDIR/analysis
  ALLSF="${LABEL_FG[*]} ${EXTRAMESHES[*]}"

  PREFIX=DFL${expid}
  # Iterate over side and subject
  for side in left right; do
    for sub in $ALLSF; do

      id=$(ls $ASHSRUNDIR | grep $id)

      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
        -N "${PREFIX}_${side}_${sub}" \
        $0 DispFinalTemp_sub $side $sub

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
    -hold_jid "${PREFIX}_*" -sync y -b y \
    sleep 1
}

function DispFinalTemp_sub()
{
  side=$1
  sub=$2
  IDS=$(cat $SUBJ_TXT)

  # Displacement analysis
  MESHES=$(for id in $IDS; do \
    echo $(find $FINALMESHDIR | grep $id | grep $side | grep ${sub}_tempfit ); done)

  meshdisp $MESHES \
    $FINALMESHDIR/analysis/disp_${side}_${sub}.vtk

  # Thickness analysis
  MESHES=$(for id in $IDS; do \
    echo $(find $FINALMESHDIR | grep $id | grep $side | grep ${sub}_thickmap ); done)

  mesh_merge_arrays \
    -r $FINALMESHDIR/template_${side}_${sub}.vtk \
    $FINALMESHDIR/analysis/thick_${side}_${sub}.vtk \
    Thickness $MESHES
}

##############################################################################
ANGRP[0]="MRG DG"
ANGRP[1]="CA DG SUB ERC BA35 BA36"

GRPNM[0]="merged"
GRPNM[1]="all"

function ThickFinalTemp()
{
  mkdir -p $FINALMESHDIR/analysis

  PREFIX=THICK${expid}
  for side in left right; do

    for design in $(ls $ALLSTATDIR | grep "design_.*txt" | grep $side ); do

      exp=$(echo $design | sed -e "s/^design_${side}_//" | sed -e "s/\.txt//")

      for ((igrp=0;igrp<${#ANGRP[*]};igrp++)); do

        # Submit job for this subject
        qsub -cwd -o $DUMPDIR -j y \
          -N "${PREFIX}_${exp}_${side}_${GRPNM[igrp]}" \
          $0 ThickFinalTemp_sub $exp $side $igrp

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function ThickFinalTemp_sub()
{
  exp=$1
  side=$2
  igrp=$3

  MYGRP=${ANGRP[igrp]}
  GNAME=${GRPNM[igrp]}

  # Create the work directory for this analysis
  WORK=$FINALMESHDIR/analysis/design_${exp}_group_${GNAME}
  mkdir -p $WORK
  rm -rf $WORK/*

  # Merge the meshes for this analysis
  for sub in $MYGRP; do

    # Get the list of subjects
    DESIGNTXT=$ALLSTATDIR/design_${side}_${exp}.txt
    SUBJ=$(cat $DESIGNTXT | awk '{print $1}')

    # Generate the design matrix for meshglm
    cat $DESIGNTXT | awk '{$1=""; print}' > $WORK/design_${side}.txt

    MESHES=$(for id in $SUBJ; do \
      echo $(find $FINALMESHDIR | grep ${id} | grep ${side}_${sub}_thickmap.vtk); done)

    mesh_merge_arrays \
      -r $FINALMESHDIR/template_${side}_${grp}_${sub}.vtk \
      $WORK/thick_${side}_${sub}.vtk Thickness $MESHES

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
    MESHPARAM=""
    for sub in $MYGRP; do
      MESHPARAM="$MESHPARAM -m $WORK/thick_${side}_${sub}.vtk $CWORK/thickstat_${FULLNM}_${side}_${sub}.vtk"

      meshglm $MESHPARAM \
        -g $WORK/design_${side}.txt $CWORK/contrast.txt \
        -a Thickness

    done
  done
}

##############################################################################
function reset_dir()
{
  rm -rf $DUMPDIR/*
}

##############################################################################
if [[ $1 == "copy_subject" ]]; then

  copy_subject $2 $3

elif [[ $1 == "pairwise_sub" ]]; then

  pairwise_sub $2 $3

elif [[ $1 == "similarity_sub" ]]; then

  similarity_sub $2 $3

elif [[ $1 == "clustering_sub" ]]; then

  clustering_sub $2

elif [[ $1 == "average_subfield" ]]; then

  average_subfield $2 $3 $4

elif [[ $1 == "ants_iter" ]]; then

  ants_iter $2 $3 $4 $5

elif [[ $1 == "warp_labels_subj" ]]; then

  warp_labels_subj $2 $3 $4

elif [[ $1 == "warp_meshes_subj" ]]; then

  warp_meshes_subj $2 $3 $4

elif [[ $1 == "eval_subj" ]]; then

  eval_subj $2 $3 $4 $5

elif [[ $1 == "disp_stats_sub" ]]; then

  disp_stats_sub $2 $3 $4

elif [[ $1 == "thick_stats_sub" ]]; then

  thick_stats_sub $2 $3 $4 $5

elif [[ $1 == "Insertion_sub" ]]; then

  Insertion_sub $2 $3 $4

elif [[ $1 == "TempSim_sub" ]]; then

  TempSim_sub $2 $3

elif [[ $1 == "InsertGraph_sub" ]]; then

  InsertGraph_sub $2

elif [[ $1 == "PropagateTemplates_sub" ]]; then

  PropagateTemplates_sub $2 $3 $4

elif [[ $1 == "WarpSuperTemplate_sub" ]]; then

  WarpSuperTemplate_sub $2 $3 $4 $5

elif [[ $1 == "RegisterSuperTemplate_sub" ]]; then

  RegisterSuperTemplate_sub $2 $3 $4 $5

elif [[ $1 == "WarpFinalLabel_sub" ]]; then

  WarpFinalLabel_sub $2 $3

elif [[ $1 == "WarpFinalMesh_sub" ]]; then

  WarpFinalMesh_sub $2 $3

elif [[ $1 == "EvalFinalTemp_sub" ]]; then

  EvalFinalTemp_sub $2 $3 $4

elif [[ $1 == "GraphPaths_sub" ]]; then

  GraphPaths_sub $2 $3 $4

elif [[ $1 == "PWSim_sub" ]]; then

  PWSim_sub $2 $3 $4

elif [[ $1 == "ChooseBestGraph_sub" ]]; then

  ChooseBestGraph_sub $2

elif [[ $1 == "AverageSuperTemplate_sub" ]]; then

  AverageSuperTemplate_sub $2 $3

elif [[ $1 == "DispFinalTemp" ]]; then

  DispFinalTemp $2 $3

elif [[ $1 == "ThickFinalTemp_sub" ]]; then

  ThickFinalTemp_sub $2 $3 $4

else

  main

  exit

fi


