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
expid=004
expdir=$WORKDIR/SCexp/exp${expid}
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
ANTs_i_1="20x10x2"
ANTs_all_metrics_1="--use-all-metrics-for-convergence"
PWDIR=${expdir}/pairwise/

# 1.2 similarity type
SIM_TYPE="PRCCS_seg_dice"
SIM_DIR=$PWDIR/sim_${SIM_TYPE}

# 1.3 clustering
INITGRP_DIR=$expdir/clustering/InitGroup/

# 1.4 prepare statitics documents
GPSTATDIR=$expdir/stat/groups/
ALLSTATDIR=$expdir/stat/all/

#########################################
# 2 generate three templates
# Relevant labels
LABEL_FG=(        CA  DG  SUB ERC BA35 BA36       PHC CS)
EXTRAMESHES=(HIPPO PRC ExtHippo MRG)
#                BKG CA DG SUB ERC BA35 BA36 PHC CS
EXTRAMESHESDEF=("-1   1  1  -1  -1  -1   -1   -1 -1" \
                "-1  -1 -1  -1  -1   1    1   -1 -1" \
                "-1  -1 -1   1   1   1    1    1 -1" \
                "-1   1 -1   1   1   1    1    1 -1")

# 2.1 initilization
groups=(1 2 3)
INITTEMP_DIR=$expdir/InitTemps/

# 2.2 ANTs
ITER=8
ANTs_start=3
WGT_2=1
REG_LABELS_2=(${LABEL_IDS[*]})
ANTs_t_2="SyN[0.25]"
ANTs_r_2="Gauss[0.5,0]"
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
EVALLABELS=(CA DG SUB ERC BA35 BA36 PRC   PHC CS HIPP     EXPHIPP   ALL)
RANGES=(    1  2  3   4   5    6    "5 6" 7   8  "1 2 3"  "4 5 6 7" "1 2 3 4 5 6 7")
MESH_EVAL=(       CA  DG  SUB ERC BA35 BA36 PRC   PHC CS)
MESH_EVAL_RANGES=(1   2   3   4   5    6    "5 6" 7   8)
INITEVALDIR=$INITTEMP_DIR/evaluation

# 2.6 

# 2.7

#########################################
# 3 insert three tempaltes into the graph

# 3.3 
FINALGRAPHDIR=$expdir/clustering/finalgraph/



#########################################
# 4 build graph and generate super template
FINALTEMPDIR=$expdir/FinalTemp

# 4.1 
ALLPATHSDIR=$FINALTEMPDIR/paths/allpaths

# 4.3 
BESTPATHSDIR=$FINALTEMPDIR/paths/bestpaths

# 4.6
FINALITER=8
ANTs_start_final=3
WGT_3=1
REG_LABELS_3=(${LABEL_IDS[*]})
ANTs_t_3="SyN[0.25]"
ANTs_r_3="Gauss[0.5,0]"
ANTs_i_3="80x80x20"
ANTs_x_3="Y"
ANTs_all_metrics_3="--use-all-metrics-for-convergence"


FINALEVALDIR=$FINALTEMPDIR/evaluation

##############################################################################
function main()
{
  reset_dir

  #######################################
  # 0. preparation
  #copy_data

  #######################################
  # 1 cluster subjects into three groups
  # 1.1 pairwise registration
  #pairwise
  #completeness

  # 1.2 compute similarity between subjects
  #similarity
  
  # 1.3 dimension reduction and perform clustering
  #clustering

  # 1.4 prepare the design and contrast documents
  #PrepStat

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

  ########################################
  # 3 insert three templates into the graph
  # 3.1 perform pairwise registration
  #Insertion

  # 3.2 measure similarity between templates and subjects
  #TempSim

  # 3.3 insert the three templates into the graph 
  #InsertGraph

  ########################################
  # 4 build graph and generate super template
  # 4.1 compute paths with different parameters


  # 4.2 apply warps and compute dice overlap


  # 4.3 choose the best graph and generate the pathes between templates
 

  # 4.4 propagate three templates to build super template
  #PropagateTemplates  

  # 4.5 warp all the subjects to the super template space
  #WarpSuperTemplate

  # 4.6 register to the final super template
  #RegisterSuperTemplate

  # 4.6 make final images
  #MakeFinalImages

  # 4.7 warp super template labels
  #WarpFinalLabel

  # 4.8 warp super template meshes
  #WarpFinalMesh
  
  # 4.9 evaluate super template
  #EvalFinalTemp

  # 4.10 displacement
  #DispFinalStat

  # 4.11 thickness analysis
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

      # Perform ANTs registration
      ANTS 3 $CMD \
           -t $ANTs_t_1 \
           -r $ANTs_r_1 \
           -i $ANTs_i_1 \
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
     clustering_ProbPCA\
     ${SIM_DIR}/sim_${SIM_TYPE}_${side}.txt \
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

#############################################################################
# Take the one that is most similar to others as the initial template
function average_subfield()
{
  side=$1
  grp=$2
  kind=$3

  GRP_IDS=($(cat $INITGRP_DIR/subjID_${side}_${grp}.txt))  

  AverageImages 3 $INITTEMP_DIR/group${grp}/work/template_fullchunk_${side}_${grp}_${kind}.nii.gz 0 $(for id in ${GRP_IDS[*]}; do echo $DATADIR/*${id}*${side}_${kind}.nii.gz; done)
}

function initialization()
{
  # Initialization for all the groups
  for grp in ${groups[*]}; do

    GROUPDIR=$INITTEMP_DIR/group${grp}
    mkdir -p ${GROUPDIR}/work

    # Compute initial average for each subfield mask to generate mask
    PREFIX=AVG${expid}
    for side in left right; do

      # Average all the input images
      for kind in $KINDS; do
        qsub -cwd -o $DUMPDIR -j y -N \
             "${PREFIX}_${kind}_${side}_${grp}" $0 \
             average_subfield $side $grp $kind
      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # Compute the mask for all the groups and cp init template
  for grp in ${groups[*]}; do

    i=$((grp-1))
    for side in left right; do

      INIT_ID=($(cat $INITGRP_DIR/AutoCenters_${side}.txt)) 
      INIT=${INIT_ID[i]}

      id=$(ls $ASHSRUNDIR | grep $INIT)

      # Compute the initial mask for the segmentation
      c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITTEMP_DIR/group${grp}/work/template_fullchunk_${side}_${grp}_${sub}.nii.gz; done | grep -v BKG) \
      -mean -thresh 1e-5 inf 1 0 -trim 5vox \
      -o $INITTEMP_DIR/group${grp}/work/template_mask_${side}_${grp}.nii.gz

      # Trim every template component using the mask
      for kind in $KINDS; do

        c3d $INITTEMP_DIR/group${grp}/work/template_mask_${side}_${grp}.nii.gz \
          $DATADIR/${id}_${side}_${kind}.nii.gz \
          -reslice-identity \
          -o $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_${kind}.nii.gz

      done
    done
  done
}

##############################################################################
# Build template using ANTs
function shape_update_to_template()
{

  side=$1
  grp=$2

  # Borrowed from ANTS buildtemplateparallel

  # Average the warp fields
  local TEMPWARP=$INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_warp.nii.gz
  $NEWANTSDIR/AverageImages 3 $TEMPWARP 0 \
    $INITTEMP_DIR/group${grp}/work/*_${side}_${grp}_totempWarp.nii.gz

  # Scale by -0.25 (gradient step)
  MultiplyImages 3 $TEMPWARP -0.25 $TEMPWARP

  # Create a new template affine
  local TEMPAFF=$INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_Affine.txt
  if [[ -f $TEMPAFF ]]; then rm -rf $TEMPAFF; fi

  cat $INITTEMP_DIR/group${grp}/work/*_${side}_${grp}_totempAffine.txt \
    | grep '^Parameters:' | awk '\
    BEGIN { for(i=0;i<12;i++) x[i]=0} \
    { for(i=0;i<12;i++) x[i]=x[i]+$(i+2) } \
    END { \
      printf "Transform: MatrixOffsetTransformBase_double_3_3\nParameters: "; \
      for(i=0;i<12;i++) printf "%f ",x[i]/NR; \
      printf "\nFixedParameters: 0 0 0\n";}' > $TEMPAFF

  # Compose the warps
  WarpImageMultiTransform 3 $TEMPWARP $TEMPWARP -i $TEMPAFF \
    -R $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_tse.nii.gz

  TEMPWARPFULL=$INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_fullwarp.nii.gz
  $ANTSPATH/ComposeMultiTransform 3 \
    $TEMPWARPFULL -R $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_tse.nii.gz \
    -i $TEMPAFF $TEMPWARP $TEMPWARP $TEMPWARP $TEMPWARP

  # Apply this warp to all the template derivatives
  for kind in $KINDS; do

    WarpImageMultiTransform 3 \
      $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_${kind}.nii.gz \
      $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_${kind}.nii.gz \
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
    $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_seg.nii.gz \
    $DATADIR/${id}_${side}_seg.nii.gz \
    $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_mlaffine.txt

  # Do we do ANTS or not?
  if [[ $doants -eq 0 ]]; then

    for sub in $KINDS; do

      c3d \
        $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_tse.nii.gz \
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
      CMD="$CMD -m MSQ[$INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_${sub}.nii.gz,$DATADIR/${id}_${side}_${sub}.nii.gz,$WGT_2]"
    done

    if [[ $ANTs_x_2 == "Y" ]]; then
       ANTs_mask_2="-x $INITTEMP_DIR/group${grp}/work/template_mask_${side}_${grp}.nii.gz"
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
      -o $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp.nii.gz \
      | tee $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_antsoutput.txt

    for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        $DATADIR/${id}_${side}_${sub}.nii.gz \
        $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp_reslice_${sub}.nii.gz \
        -R $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_tse.nii.gz \
        $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempWarp.nii.gz \
        $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempAffine.txt

    done

  fi
}

function main_loop()
{
  # Main iteration loop
  PREFIX=ANTs${expid}
  for side in left right; do

    for ((iter=0;iter<$ITER;iter++)); do

      for grp in ${groups[*]}; do

        # Create the segmentation for the template
        c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_${sub}.nii.gz; done) \
          -vote -type ushort \
          -o $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_seg.nii.gz

        # Back up template
        ITDIR=$INITTEMP_DIR/group${grp}/work/$(printf iter_%s_%02d $side $iter)
        mkdir -p $ITDIR
        cp -a $INITTEMP_DIR/group${grp}/work/template_${side}_*.nii.gz $ITDIR/

        # Do ants?
        if [[ $iter -lt $ANTs_start ]]; then doants=0; else doants=1; fi

        #IDS=$(ls ${expdir}/data | grep ${side}_${grp}_tse | sed -e "s/_${side}_${grp}_tse.nii.gz//")
        GRP_IDS=($(cat $INITGRP_DIR/subjID_${side}_${grp}.txt))

        # Run ANTS for each image
        for grpid in ${GRP_IDS[*]}; do

          id=$(ls $ASHSRUNDIR | grep $grpid)

          #Submit ANTS job
          qsub -cwd -o $DUMPDIR -j y -N \
               "${PREFIX}_${id}_${side}_${grp}" $0 \
               ants_iter $id $side $grp $doants

        done

      done

      # Wait for completion
      qsub -cwd -o $DUMPDIR -j y \
        -hold_jid "${PREFIX}_*" -sync y -b y \
        sleep 1

      # If this is the last iteration, we don't want to recompute the template
      if [[ $iter -lt $((ITER-1)) ]]; then

        for grp in ${groups[*]}; do

          # Compute average images
          for kind in $KINDS; do

            if [[ $kind = "tse" ]]; then NORM=1; else NORM=0; fi
            AverageImages 3 $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_${kind}.nii.gz $NORM $INITTEMP_DIR/group${grp}/work/*_${side}_${grp}_totemp_reslice_${kind}.nii.gz

          done

          # Perform shape averaging
          if [[ $doants -eq 1 ]]; then
            shape_update_to_template $side $grp
          fi

        done

      fi

    done
  done
}

##############################################################################
# Make images for further analysis
function make_images()
{
  for grp in ${groups[*]}; do

    for side in left right; do

      # Generate one compound label for all other subfields
      ITDIR=$INITTEMP_DIR/group${grp}/work/$(printf iter_%s_%02d $side $((ITER-1)) )
      for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

        LABELDEF=(${EXTRAMESHESDEF[i]})
        echo $LABELDEF

        c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
          -accum -add -endaccum \
          -o $ITDIR/template_${side}_${grp}_${EXTRAMESHES[i]}.nii.gz

      done

      # for each subject
      #IDS=$(ls ${expdir}/data | grep ${side}_${grp}_tse \
      #  | sed -e "s/_${side}_${grp}_tse.nii.gz//")
      GRP_IDS=($(cat $INITGRP_DIR/subjID_${side}_${grp}.txt))

      for grpid in ${GRP_IDS[*]}; do

        id=$(ls $ASHSRUNDIR | grep $grpid)

        # Create seg for each subject
        c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp_reslice_${sub}.nii.gz; done) -vote -type ushort \
          -o $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp_reslice_seg.nii.gz

      done
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

    # Iterate over side and subject
    for side in left right; do

      #IDS=$(ls ${expdir}/data | grep ${side}_${grp}_tse | sed -e "s/_${side}_${grp}_tse.nii.gz//")
      GRP_IDS=($(cat $INITGRP_DIR/subjID_${side}_${grp}.txt))

      for grpid in ${GRP_IDS[*]}; do

        id=$(ls $ASHSRUNDIR | grep $grpid)

        # Submit job for this subject
        qsub -cwd -o $DUMPDIR -j y \
          -N "${PREFIX}_${id}_${side}_${grp}" \
          $0 warp_labels_subj $id $side $grp

      done
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
  grp=$3

  # Apply the warp to each of the labels
  for sub in ${LABEL_IDS[*]}; do

    antsApplyTransforms -d 3 \
      -r $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray.nii.gz \
      -i $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_${sub}.nii.gz \
      -o $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${sub}_tempfit.nii.gz \
      -t [$ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt,1] \
      -t [$ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt,1] \
      -t $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii \
      -t [$INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempAffine.txt,1] \
      -t $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempInverseWarp.nii.gz \
      -n BSpline

  done

  # Vote in original subject space to get the tempfit segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${sub}_tempfit.nii.gz; done) -vote -type ushort \
    -o $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_tempfit.nii.gz


  # change ASHS labels to the current formate
  for ((i=0;i<${#LABEL_IDS[*]};i++)); do
    c3d $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray.nii.gz \
      -replace $(for k in ${LABEL_MRG[i]}; do echo $k 999; done) \
      -thresh 999 999 1 0 \
      -o $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${LABEL_IDS[i]}_ashs.nii.gz
  done

  c3d $(for sub in ${LABEL_IDS[*]}; do echo $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_${sub}_ashs.nii.gz; done) -vote -type ushort \
    -o $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_ashs.nii.gz
}

##############################################################################
# Extract meshes for the subfields and apply warps to these meshes. This allows us to# perform statistical analysis on the mesh boundaries
function warp_meshes()
{
  PREFIX=WM${expid}
  for grp in ${groups[*]}; do

    mkdir -p $INITTEMP_DIR/group${grp}/meshwarp
    mkdir -p $INITTEMP_DIR/group${grp}/jacobian
    mkdir -p $INITTEMP_DIR/group${grp}/meshwarp/template

    # Iterate over side and subject
    for side in left right; do

      # Generate meshes for the individual subfields
      for sub in ${LABEL_FG[*]}; do

        vtklevelset \
          $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_${sub}.nii.gz \
          $INITTEMP_DIR/group${grp}/meshwarp/template_${side}_${grp}_${sub}.vtk \
          $subfield_th

      done

      # Generate one mesh for all non-DG non-CS subfields
      ITDIR=$INITTEMP_DIR/group${grp}/work/$(printf iter_%s_%02d $side $((ITER-1)) )
      for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

        vtklevelset \
          $ITDIR/template_${side}_${grp}_${EXTRAMESHES[i]}.nii.gz \
          $INITTEMP_DIR/group${grp}/meshwarp/template_${side}_${grp}_${EXTRAMESHES[i]}.vtk \
          $template_th

      done

      #IDS=$(ls ${expdir}/data | grep ${side}_${grp}_tse | sed -e "s/_${side}_${grp}_tse.nii.gz//")
      GRP_IDS=($(cat $INITGRP_DIR/subjID_${side}_${grp}.txt))

      for grpid in ${GRP_IDS[*]}; do

        id=$(ls $ASHSRUNDIR | grep $grpid)

        # Submit job for this subject
        qsub  -cwd -o $DUMPDIR -j y \
          -N "${PREFIX}_${id}_${side}_${grp}" \
          $0 warp_meshes_subj $id $side $grp

      done

    done

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
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

  # Compose the transformation between the template and the subject
  $ANTSPATH/ComposeMultiTransform 3 \
     $TMPDIR/compose.nii \
     -R $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_tse.nii.gz \
     $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempWarp.nii.gz \
     $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totempAffine.txt \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempWarp.nii \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
     $ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

  # Split the transformations into xyz
  c3d -mcs $TMPDIR/compose.nii -oo $TMPDIR/comp%02d.nii

  # Apply the warp to each of the meshes
  for sub in $ALLSF; do

    # Warp the subfield into subject space
    warpmesh -w ants -m ras \
      $INITTEMP_DIR/group${grp}/meshwarp/template_${side}_${grp}_${sub}.vtk \
      $INITTEMP_DIR/group${grp}/meshwarp/${id}_${side}_${grp}_${sub}_tempfit.vtk \
      $TMPDIR/comp00.nii $TMPDIR/comp01.nii $TMPDIR/comp02.nii

    # Extract the thickness of the subfield
    cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
      -T $INITTEMP_DIR/group${grp}/meshwarp/${id}_${side}_${grp}_${sub}_thickmap.vtk \
      -p $thick_p -e $thick_e \
      $INITTEMP_DIR/group${grp}/meshwarp/${id}_${side}_${grp}_${sub}_tempfit.vtk \
      $INITTEMP_DIR/group${grp}/meshwarp/skel_${id}_${side}_${grp}_${sub}.vtk

  done

  # Compute the Jacobian map
  ANTSJacobian 3 $TMPDIR/compose.nii \
    $INITTEMP_DIR/group${grp}/jacobian/${id}_${side}_${grp}_totemp 1
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

    # Iterate over side and subject
    for side in left right; do

      #IDS=$(ls ${expdir}/data | grep ${side}_${grp}_tse | sed -e "s/_${side}_${grp}_tse.nii.gz//")
      GRP_IDS=($(cat $INITGRP_DIR/subjID_${side}_${grp}.txt))

      for grpid in ${GRP_IDS[*]}; do

        id=$(ls $ASHSRUNDIR | grep $grpid)

        # Submit job for this subject
        qsub -cwd -o $DUMPDIR -j y \
          -N "${PREFIX}_${id}_${side}_${grp}" \
          $0 eval_subj $id $side $grp $INITEVALTMPDIR
        sleep 0.1

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # Combine scores
  IDs=($(cat $SUBJ_TXT))
  for side in left right; do

    echo "ID ${EVALLABELS[*]} ${MESH_EVAL[*]} ${EVALLABELS[*]}" \
      > $INITEVALDIR/overlap_${side}.txt

    for ids in ${IDs[*]}; do

      id=$(ls $ASHSRUNDIR | grep $ids)
      cat $INITEVALTMPDIR/${id}_${side}_overlap.txt \
          >> $INITEVALDIR/overlap_${side}.txt

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
  do_pair $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_seg.nii.gz \
    $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp_reslice_seg.nii.gz
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
}

##############################################################################
# Display statistics
ALLSF="${LABEL_FG[*]} MRG"
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
          -N "${PREFIX}_${side}_${grp}_${sub}" \
          $0 disp_stats_sub $side $grp $sub

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function disp_stats_sub()
{
  side=$1
  grp=$2
  sub=$3

  #IDS=$(ls ${expdir}/data | grep ${side}_${grp}_tse | sed -e "s/_${side}_${grp}_tse.nii.gz//")
  IDS=$(cat $INITGRP_DIR/subjID_${side}_${grp}.txt)

  # Displacement analysis
  MESHES=$(for id in $IDS; do \
    echo $(find $INITTEMP_DIR/group${grp}/meshwarp | grep $id | grep $side | grep ${sub}_tempfit ); done)

  meshdisp $MESHES \
    $INITTEMP_DIR/group${grp}/meshwarp/analysis/disp_${side}_${sub}.vtk

  # Thickness analysis
  MESHES=$(for id in $IDS; do \
    echo $(find $INITTEMP_DIR/group${grp}/meshwarp | grep $id | grep $side | grep ${sub}_thickmap ); done)

  mesh_merge_arrays \
    -r $INITTEMP_DIR/group${grp}/meshwarp/template_${side}_${grp}_${sub}.vtk \
    $INITTEMP_DIR/group${grp}/meshwarp/analysis/thick_${side}_${sub}.vtk \
    Thickness $MESHES
}

##############################################################################
# For thickness analysis, we can run analysis on multiple groups of meshes at once
ANGRP[0]="MRG DG"
ANGRP[1]="CA DG SUB ERC BA35 BA36"

GRPNM[0]="merged"
GRPNM[1]="all"

function thick_stats()
{
  PREFIX=THICK${expid}
  for grp in ${groups[*]}; do

    for side in left right; do

     for design in $(ls $GPSTATDIR | grep "design_.*txt" | grep $side | grep $grp ); do

        exp=$(echo $design | sed -e "s/^design_${side}_${grp}_//" | sed -e "s/\.txt//")

        for ((igrp=0;igrp<${#ANGRP[*]};igrp++)); do

          # Submit job for this subject
          qsub -cwd -o $DUMPDIR -j y \
            -N "${PREFIX}_${exp}_${side}_${grp}_${GRPNM[igrp]}" \
            $0 thick_stats_sub $exp $side $grp $igrp

        done
      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
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
  rm -rf $WORK/*

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
      $INITTEMP_DIR/group${grp}/meshwarp/template_${side}_${grp}_${sub}.vtk \
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

      meshglm $MESHPARAM \
        -g $WORK/design_${side}.txt $CWORK/contrast.txt \
        -a Thickness

    done

  done
}

##############################################################################
function Insertion()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=IS${expid}
  for side in left right; do
    for grp in ${groups[*]}; do      

      # iter directory
      ITDIR=$INITTEMP_DIR/group${grp}/work/$(printf iter_%s_%02d $side $((ITER-1)) )

      # link templates to data folder
      for kind in $KINDS seg; do
        
        ln -sf $ITDIR/template_${side}_${grp}_${kind}.nii.gz \
          $DATADIR/template${grp}_${side}_${kind}.nii.gz

      done

      for id in $IDS; do

        fn_id=$(ls $ASHSRUNDIR | grep $id)
        # perform registration from templates to subjects
        qsub -cwd -o $DUMPDIR -j y \
           -N "${PREFIX}_${id}_${grp}_${side}" $0 \
           Insertion_sub template${grp} $fn_id $side
        sleep 0.1

        # perform registration from subjects to templates
        qsub -cwd -o $DUMPDIR -j y \
           -N "${PREFIX}_${grp}_${id}_${side}" $0 \
           Insertion_sub $fn_id template${grp} $side
        sleep 0.1

      done
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
  side=$3
  OUTDIR=$PWDIR/${id_fix}_${side}/${id_mov}_to_${id_fix}
  mkdir -p $OUTDIR

  if [[ -f $OUTDIR/${id_mov}_to_${id_fix}_${side}_reslice_seg.nii.gz ]]; then

    echo "Seg file exists."

  else

    # Use ml_affine for nice affine alignment
    /data/picsl/pauly/wolk/ashs/ext/Linux/bin/ml_affine \
      $DATADIR/${id_fix}_${side}_seg.nii.gz \
      $DATADIR/${id_mov}_${side}_seg.nii.gz \
      $OUTDIR/${id_mov}_to_${id_fix}_${side}_mlaffine.txt

    # Convert that to ITK format
    c3d_affine_tool \
      $OUTDIR/${id_mov}_to_${id_fix}_${side}_mlaffine.txt \
      -oitk $OUTDIR/${id_mov}_to_${id_fix}_${side}_mlaffine_itk.txt

    CMD=""
    for sub in ${REG_LABELS_1[*]}; do
      CMD="$CMD -m MSQ[$DATADIR/${id_fix}_${side}_${sub}.nii.gz,$DATADIR/${id_mov}_${side}_${sub}.nii.gz,$WGT_1]"
    done

    # Perform ANTs registration
    ANTS 3 $CMD \
         -t $ANTs_t_1 \
         -r $ANTs_r_1 \
         -i $ANTs_i_1 \
         -a $OUTDIR/${id_mov}_to_${id_fix}_${side}_mlaffine_itk.txt \
         --continue-affine 0 \
         $ANTs_all_metrics_1 \
         -o $OUTDIR/${id_mov}_to_${id_fix}_${side}_pairwise.nii.gz \
         | tee $OUTDIR/${id_mov}_to_${id_fix}_${side}_antsoutput.txt

    for sub in ${KINDS[*]}; do

      WarpImageMultiTransform 3 \
        $DATADIR/${id_mov}_${side}_${sub}.nii.gz \
        $OUTDIR/${id_mov}_to_${id_fix}_${side}_reslice_${sub}.nii.gz \
        -R $DATADIR/${id_fix}_${side}_tse.nii.gz \
        $OUTDIR/${id_mov}_to_${id_fix}_${side}_pairwiseWarp.nii.gz \
        $OUTDIR/${id_mov}_to_${id_fix}_${side}_pairwiseAffine.txt

    done

    # Create seg
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $OUTDIR/${id_mov}_to_${id_fix}_${side}_reslice_${sub}.nii.gz; done) \
        -vote -type ushort \
        -o $OUTDIR/${id_mov}_to_${id_fix}_${side}_reslice_seg.nii.gz

  fi
}

##############################################################################
function TempSim()
{
  mkdir -p $SIM_DIR

  PREFIX=TS${expid}
  for side in left right; do
    for grp in ${groups[*]}; do

      qsub -cwd -o $DUMPDIR -j y \
           -l h_vmem=2.1G,s_vmem=2G \
           -N "${PREFIX}_${grp}_${side}" $0 \
           TempSim_sub $grp $side
      sleep 0.1

    done
  done

  # wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function TempSim_sub()
{
  grp=$1
  side=$2
  fn_fix=template${grp}
  fn_out=$SIM_DIR/${fn_fix}_${side}.txt
  fn_seg_fix=$TMPDIR/${fn_fix}_${side}_seg_tmp.nii.gz
  rm -rf $fn_out

  # preprocess the template seg
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

  fi

  # go through all the subjects
  OVL_1=""
  OVL_2=""
  IDS=$(cat $SUBJ_TXT)
  for id_mov in $IDS; do

    fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)

    # subject to template overlap
    OVL_tmp="NA"
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

    fi

    OVL_1="$OVL_1 $OVL_tmp"

    # template to subject overlap
    OVL_tmp="NA"
    if [ ${SIM_TYPE} == "PRC_seg_dice" ]; then

      OVL_tmp=$(c3d $PWDIR/${fn_mov}_${side}/${fn_fix}_to_${fn_mov}/${fn_fix}_to_${fn_mov}_${side}_reslice_seg.nii.gz \
          -thresh 5 6 1 0 \
          $DATADIR/${fn_mov}_${side}_seg.nii.gz \
          -replace 5 6 1 0 \ 
          -overlap 1 \
          | grep OVL | awk -F '[ ,]+' '{print $6}' )

    elif [ ${SIM_TYPE} == "CS_seg_dice" ]; then

      OVL_tmp=$(c3d $PWDIR/${fn_mov}_${side}/${fn_fix}_to_${fn_mov}/${fn_fix}_to_${fn_mov}_${side}_reslice_seg.nii.gz \
          -thresh 8 8 1 0 \
          $DATADIR/${fn_mov}_${side}_seg.nii.gz \
          -replace 8 8 1 0 \
          -overlap 1 \
          | grep OVL | awk -F '[ ,]+' '{print $6}' )

    elif [ ${SIM_TYPE} == "PRCCS_seg_dice" ]; then

      OVL_tmp=$(c3d $PWDIR/${fn_mov}_${side}/${fn_fix}_to_${fn_mov}/${fn_fix}_to_${fn_mov}_${side}_reslice_seg.nii.gz \
          -replace 1 0 2 0 3 0 4 0 7 0 \
          $DATADIR/${fn_mov}_${side}_seg.nii.gz \
          -replace 1 0 2 0 3 0 4 0 7 0 \
          -label-overlap \
          | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}' )

    fi

    OVL_2="$OVL_2 $OVL_tmp"

  done

  # output similarity
  echo $OVL_1 > $fn_out
  echo $OVL_2 >> $fn_out
}

##############################################################################
function InsertGraph()
{
  mkdir -p $FINALGRAPHDIR

  PREFIX=IG${expid}
  for side in left right; do

      qsubp5 -cwd -o $DUMPDIR -j y \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${side}" $0 \
           InsertGraph_sub $side
      sleep 0.1

  done

  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function InsertGraph_sub()
{
  side=$1

  # Run matlab to denoise T1 image
  mkdir -p $TMPDIR/${id}
  source $PKGDIR/matlab_batch.sh \
     $TMPDIR/${id}/ \
     $MATLABCODEDIR \
     insert_template \
     $INITGRP_DIR/allinfo_${side}.mat \
     $SIM_DIR/ \
     $side \
     $FINALGRAPHDIR/graph_${side}.mat
}

##############################################################################
function PropagateTemplates()
{
  PREFIX=PT${expid}
  for side in left right; do
    for grp_fix in ${groups[*]}; do
      for grp_mov in ${groups[*]}; do

      if [[ $grp_fix == $grp_mov ]]; then
        continue
      fi

      # submit jobs
      qsub -cwd -o $DUMPDIR -j y \
        -l h_vmem=4.1G,s_vmem=4G \
        -N "${PREFIX}_${side}_${grp_fix}_${grp_mov}" $0 \
        PropagateTemplates_sub $side $grp_fix $grp_mov
      sleep 0.1
  
      done
    done
  done

  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function PropagateTemplates_sub()
{
  side=$1
  grp_fix=$2
  grp_mov=$3
  IDS=($(cat $SUBJ_TXT))

  # get the path
  paths=($(cat $BESTPATHSDIR/$side/template${grp_mov}_${side}.txt))
  N=${#IDS[*]}
  N_grps=${#groups[*]}
  N_all=${#paths[*]}
  idx=$((N_all-N_grps-1+grp_fix))
  cur_path=${paths[$idx]}

  # generate IDS array with templates
  for ((i=0;i<$N;i++)); do
    id=${IDS[i]}
    IDS[$i]=$(ls $ASHSRUNDIR | grep $id)
  done

  for ((i=0;i<$N_grps;i++)); do
    k=$((i+1))
    IDS[$((N+i))]="template${k}"
  done

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
  TRANSDIR=$FINALTEMPDIR/template${grp_fix}/tempwarps
  mkdir -p $TRANSDIR
  $ANTSPATH/ComposeMultiTransform 3 \
    $TRANSDIR/template${grp_mov}_${side}_InitWarp.nii.gz \
    -R $DATADIR/template${grp_fix}_${side}_tse.nii.gz \
    $warps

  $ANTSPATH/ComposeMultiTransform 3 \
    $TRANSDIR/template${grp_mov}_${side}_InitInverseWarp.nii.gz \
    -R $DATADIR/template${grp_mov}_${side}_tse.nii.gz \
    $invwarps

  # apply init warps
  for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        $DATADIR/template${grp_mov}_${side}_${sub}.nii.gz \
        $TRANSDIR/template${grp_mov}_${side}_inittotemp_reslice_${sub}.nii.gz \
        -R $DATADIR/template${grp_fix}_${side}_tse.nii.gz \
        $TRANSDIR/template${grp_mov}_${side}_InitWarp.nii.gz

  done

  # Create the segmentation for the template
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $TRANSDIR/template${grp_mov}_${side}_inittotemp_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $TRANSDIR/template${grp_mov}_${side}_inittotemp_reslice_seg.nii.gz
}

##############################################################################
function WarpSuperTemplate()
{
  # IDS
  IDS=($(cat $SUBJ_TXT))

  PREFIX=WST${expid}
  for side in left right; do

    # load groups
    GRPS=($(cat $INITGRP_DIR/group_${side}.txt))

    # for all three templates
    for grp_temp in ${groups[*]}; do

      # loop through each subject
      for ((i=0;i<${#IDS[*]};i++)); do

      id=${IDS[i]}
      grp_subj=${GRPS[i]}

      # submit jobs
      #qsub -cwd -o $DUMPDIR -j y \
      #  -N "${PREFIX}_${side}_${grp_temp}_${id}_${grp_subj}" $0 \
      #  WarpSuperTemplate_sub $side $grp_temp $id $grp_subj
      #sleep 0.1

      done
    done
  done

  # wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # generate the average template
  for side in left right; do
    for grp_temp in ${groups[*]}; do

      FINALWORKDIR=$FINALTEMPDIR/template${grp_temp}/finalwork/
      mkdir -p $FINALWORKDIR

      # Compute average images
      for kind in $KINDS; do

        if [[ $kind = "tse" ]]; then NORM=1; else NORM=0; fi
        AverageImages 3 $FINALWORKDIR/template_${side}_${kind}.nii.gz $NORM $FINALTEMPDIR/template${grp_temp}/work/*_${side}_totemp_reslice_${kind}.nii.gz

      done

      # Create the segmentation for the template
      c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALWORKDIR/template_${side}_${sub}.nii.gz; done) \
        -vote -type ushort \
        -o $FINALWORKDIR/template_${side}_seg.nii.gz

      # generate mask
      c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALWORKDIR/template_${side}_${sub}.nii.gz; done | grep -v BKG) \
        -mean -thresh 1e-5 inf 1 0 -trim 5vox \
        -o $FINALWORKDIR/template_${side}_mask.nii.gz

    done
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
      warps="$FINALTEMPDIR/template${grp_temp}/tempwarps/template${grp_subj}_${side}_InitWarp.nii.gz $warps"

      invwarps="$invwarps $FINALTEMPDIR/template${grp_temp}/tempwarps/template${grp_subj}_${side}_InitInverseWarp.nii.gz"
    fi

    # compose warps
    $ANTSPATH/ComposeMultiTransform 3 \
      $FINALTEMPLATEWORKDIR/${fn_id}_${side}_totempWarp.nii.gz \
      -R $DATADIR/template${grp_temp}_${side}_tse.nii.gz \
      $warps  

    # compose inverse warps
    $ANTSPATH/ComposeMultiTransform 3 \
      $FINALTEMPLATEWORKDIR/${fn_id}_${side}_totempInverseWarp.nii.gz \
      -R $DATADIR/${fn_id}_${side}_tse.nii.gz \
      $invwarps

    WarpImageMultiTransform 3 \
      $DATADIR/${fn_id}_${side}_${sub}.nii.gz \
      $FINALTEMPLATEWORKDIR/${fn_id}_${side}_totemp_reslice_${sub}.nii.gz \
      -R $DATADIR/template${grp_temp}_${side}_tse.nii.gz \
      $FINALTEMPLATEWORKDIR/${fn_id}_${side}_totempWarp.nii.gz

    done
}

##############################################################################
function RegisterSuperTemplate()
{
  # choose the best template to be the super template
  PREFIX=RST${expid}
  for side in left right; do

    besttemplate=""
    maxdice=0
    for grp_fix in ${groups[*]}; do

      curdice=0
      for grp_mov in ${groups[*]}; do

        if [[ $grp_mov -eq $grp_fix ]]; then
          continue
        fi

        tmpdice=$(c3d \
          $DATADIR/template${grp_fix}_${side}_seg.nii.gz \
          -replace 1 0 2 0 3 0 4 0 7 0 \
          $FINALTEMPDIR/template${grp_fix}/tempwarps/template${grp_mov}_${side}_inittotemp_reslice_seg.nii.gz \
          -replace 1 0 2 0 3 0 4 0 7 0 \
          -label-overlap \
          | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}' )

        curdice=$(echo "$curdice+$tmpdice" | bc)

      done

      st=$(echo "$curdice > $maxdice" | bc)
      if [[ $st -eq 1 ]]; then
        maxdice=$curdice
        besttemplate=$grp_fix
      fi

    done

    # link to the best template
    ln -sf $FINALTEMPDIR/template${besttemplate} \
      $FINALTEMPDIR/template_${side}
    FINALWORKDIR=$FINALTEMPDIR/template_${side}/finalwork

    # write to file
    echo $besttemplate > $FINALGRAPHDIR/besttemplate_${side}.txt
    echo $besttemplate > $FINALTEMPDIR/template_${side}/${side}_best_is_${besttemplate}.txt

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
          c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALWORKDIR/template_${side}_${sub}.nii.gz; done) \
            -vote -type ushort \
            -o $FINALWORKDIR/template_${side}_seg.nii.gz

          # Perform shape averaging
          #if [[ $doants -eq 1 ]]; then
          #  shape_update_to_template $side $grp
          #fi

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

  besttemplate=$(cat $FINALGRAPHDIR/besttemplate_${side}.txt)
  SRCWORKDIR=$FINALTEMPDIR/template${besttemplate}/work
  FINALWORKDIR=$FINALTEMPDIR/template${besttemplate}/finalwork

  if [[ ! -f $SRCWORKDIR/${id}_${side}_totemp_reslice_seg.nii.gz ]]; then
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $SRCWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $SRCWORKDIR/${id}_${side}_totemp_reslice_seg.nii.gz
  fi

  # Before we vote, use ml_affine for nice affine alignment
  ~pauly/wolk/ashs/ext/Linux/bin/ml_affine \
    $FINALWORKDIR/template_${side}_seg.nii.gz \
    $SRCWORKDIR/${id}_${side}_totemp_reslice_seg.nii.gz \
    $FINALWORKDIR/${id}_${side}_mlaffine.txt

  # Do we do ANTS or not?
  if [[ $doants -eq 0 ]]; then

    for sub in $KINDS; do

      c3d \
        $FINALWORKDIR/template_${side}_tse.nii.gz \
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
      CMD="$CMD -m MSQ[$FINALWORKDIR/template_${side}_${sub}.nii.gz,$SRCWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz,$WGT_3]"
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
        $SRCWORKDIR/${id}_${side}_totemp_reslice_${sub}.nii.gz \
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

##############################################################################
function MakeFinalImages()
{
  for side in left right; do

    # Generate one compound label for all other subfields
    ITDIR=$FINALTEMPDIR/template_${side}/finalwork/$(printf iter_%s_%02d $side $((FINALITER-1)) )

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
  mkdir -p $FINALTEMPDIR/labelwarp

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

  # Apply the warp to each of the labels
  for sub in ${LABEL_IDS[*]}; do

    antsApplyTransforms -d 3 \
      -r $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray.nii.gz \
      -i $FINALTEMPDIR/template_${side}/finalwork/template_${side}_${sub}.nii.gz \
      -o $FINALTEMPDIR/labelwarp/${id}_${side}_${sub}_tempfit.nii.gz \
      -t [$ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt,1] \
      -t [$ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt,1] \
      -t $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii \
      -t $FINALTEMPDIR/template_${side}/work/${id}_${side}_totempInverseWarp.nii.gz \
      -t [$FINALTEMPDIR/template_${side}/finalwork/${id}_${side}_totempAffine.txt,1] \
      -t $FINALTEMPDIR/template_${side}/finalwork/${id}_${side}_totempInverseWarp.nii.gz \
      -n BSpline

  done

  # Vote in original subject space to get the tempfit segmentation
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALTEMPDIR/labelwarp/${id}_${side}_${sub}_tempfit.nii.gz; done) -vote -type ushort \
    -o $FINALTEMPDIR/labelwarp/${id}_${side}_seg_tempfit.nii.gz


  # change ASHS labels to the current formate
  for ((i=0;i<${#LABEL_IDS[*]};i++)); do
    c3d $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray.nii.gz \
      -replace $(for k in ${LABEL_MRG[i]}; do echo $k 999; done) \
      -thresh 999 999 1 0 \
      -o $FINALTEMPDIR/labelwarp/${id}_${side}_${LABEL_IDS[i]}_ashs.nii.gz
  done

  c3d $(for sub in ${LABEL_IDS[*]}; do echo $FINALTEMPDIR/labelwarp/${id}_${side}_${sub}_ashs.nii.gz; done) -vote -type ushort \
    -o $FINALTEMPDIR/labelwarp/${id}_${side}_seg_ashs.nii.gz
}

##############################################################################
function WarpFinalMesh()
{
  PREFIX=WFM${expid}

  mkdir -p $FINALTEMPDIR/meshwarp
  mkdir -p $FINALTEMPDIR/jacobian
  mkdir -p $FINALTEMPDIR/meshwarp/template

  # Iterate over side and subject
  for side in left right; do

    ITDIR=$FINALTEMPDIR/template_${side}/finalwork/$(printf iter_%s_%02d $side $((FINALITER-1)) )

    # Generate meshes for the individual subfields
    for sub in ${LABEL_FG[*]}; do

      vtklevelset \
        $ITDIR/template_${side}_${sub}.nii.gz \
        $FINALTEMPDIR/meshwarp/template_${side}_${sub}.vtk \
        $subfield_th

    done

    # Generate one mesh for all non-DG non-CS subfields
    for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

      vtklevelset \
        $ITDIR/template_${side}_${EXTRAMESHES[i]}.nii.gz \
        $FINALTEMPDIR/meshwarp/template_${side}_${EXTRAMESHES[i]}.vtk \
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

  cp $FINALTEMPDIR/meshwarp/template*.vtk \
     $FINALTEMPDIR/meshwarp/template/
}

function WarpFinalMesh_sub()
{
  id=$1
  side=$2
  ALLSF="${LABEL_FG[*]} ${EXTRAMESHES[*]}"

  # Compose the transformation between the template and the subject
  $ANTSPATH/ComposeMultiTransform 3 \
     $TMPDIR/compose.nii \
     -R $FINALTEMPDIR/template_${side}/finalwork/template_${side}_tse.nii.gz \
     $FINALTEMPDIR/template_${side}/finalwork/${id}_${side}_totempWarp.nii.gz \
     $FINALTEMPDIR/template_${side}/finalwork/${id}_${side}_totempAffine.txt \
     $FINALTEMPDIR/template_${side}/work/${id}_${side}_totempWarp.nii.gz \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempWarp.nii \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
     $ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

  # Split the transformations into xyz
  c3d -mcs $TMPDIR/compose.nii -oo $TMPDIR/comp%02d.nii

  # Apply the warp to each of the meshes
  for sub in $ALLSF; do

    # Warp the subfield into subject space
    warpmesh -w ants -m ras \
      $FINALTEMPDIR/meshwarp/template_${side}_${sub}.vtk \
      $FINALTEMPDIR/meshwarp/${id}_${side}_${sub}_tempfit.vtk \
      $TMPDIR/comp00.nii $TMPDIR/comp01.nii $TMPDIR/comp02.nii

    # Extract the thickness of the subfield
    cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
      -T $FINALTEMPDIR/meshwarp/${id}_${side}_${sub}_thickmap.vtk \
      -p $thick_p -e $thick_e \
      $FINALTEMPDIR/meshwarp/${id}_${side}_${sub}_tempfit.vtk \
      $FINALTEMPDIR/meshwarp/skel_${id}_${side}_${sub}.vtk

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
  FITTED=$FINALTEMPDIR/labelwarp/${id}_${side}_seg_tempfit.nii.gz
  ASHSSEG=$FINALTEMPDIR/labelwarp/${id}_${side}_seg_ashs.nii.gz
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

      MESH_DIST_TMP="$(meshdiff $FINALTEMPDIR/meshwarp/${id}_${side}_${MESH_EVAL[i]}_tempfit.vtk $TMPDIR/mask_seg.vtk | grep RESULT | awk '{print $7}')"

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
  do_pair $FINALTEMPDIR/template_${side}/finalwork/template_${side}_seg.nii.gz \
    $FINALTEMPDIR/template_${side}/finalwork/${id}_${side}_totemp_reslice_seg.nii.gz
  ALLOVL="$ALLOVL $FULLOVL"

  ###################################
  # save to temporary file
  echo $id $ALLOVL > $EVALTMPDIR/${id}_${side}_overlap.txt
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

else

  main

  exit

fi


