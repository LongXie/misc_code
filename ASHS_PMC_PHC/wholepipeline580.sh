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
MATLAB_BIN=/share/apps/matlab/R2016b/bin/matlab
#MATLAB_BIN=/share/apps/matlab/R2013a/bin/matlab
export PATH=$ANTSPATH:$C3DPATH:$PATH

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

# Directories
ROOT=/data/jet/longxie/ASHS_PHC
WORKDIR=$ROOT/thickness_newlabel
ASHSRUNDIR=$ROOT/SuperResolution/ASHSexp/exp003/fullset/ashs
MATLABCODEDIR=$WORKDIR/matlabcode/
GTGROUPDIR=$WORKDIR/group/
SUBJ_TXT=$WORKDIR/analysis_input/subj.txt
SUBJ_WITHTEMP_TXT=$WORKDIR/analysis_input/subj_withtemp.txt
LMMATCHINGDIR=/home/longxie/LMTemplateMatching/PointSetTemplateMatching/build_release
LMTOWARPDIR=/home/longxie/LMTemplateMatching/PointSetUtilities/build_release

##############################################################################
# Parameters needs to be specify
# Experiment number
expid=580
expdir=$WORKDIR/exp/exp${expid}
DUMPDIR=${expdir}/dump
mkdir -p ${expdir}/dump

########################################
# 0. copy data
#LABEL_IDS_ALL=(BKG      CA      DG  SUB ERC  BA35 BA36 PHC  CS   OTS)
#LABEL_MRG_ALL=("0 7 16" "1 2 4" "3" "8" "10" "11" "12" "13" "14" "16")
# Labels to get segmentation
#LABEL_IDS=(BKG      CA      DG  SUB ERC  BA35 BA36 PHC  CS)
#LABEL_MRG=("0 7 16" "1 2 4" "3" "8" "10" "11" "12" "13" "14")
#LABEL_NEW=(0        1       2   3   4    5    6    7    8)


LABEL_IDS_ALL=(BKG      CA      DG  SUB ERC  BA35 BA36 PHC  aCS  pCS)
LABEL_MRG_ALL=("0 7 16" "1 2 4" "3" "8" "10" "11" "12" "13" "15" "14")
# Labels to get segmentation
LABEL_IDS=(BKG      CA      DG  SUB  ERC  BA35 BA36 PHC  aCS  pCS)
LABEL_MRG=("0 7 16" "1 2 4" "3" "8" "10" "11" "12" "13" "15" "14")
LABEL_NEW=(0        1       2   3    4    5    6    7    8    9)
MESH_LABEL=(CA DG SUB ERC  BA35 BA36 PHC)
KINDS="tse mprage ${LABEL_IDS[*]}"
DATADIR=${expdir}/data

########################################
# 1 pairwise registration and clustering
# 1.1 ANTs parameters
REUSE="N"
REUSEID="506"
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
LABEL_FG=(           CA  DG  SUB ERC BA35 BA36 PHC aCS pCS)
EXTRAMESHES=(HIPPO ExtHippo  MRG)
#                BKG CA  DG   SUB  ERC BA35 BA36 PHC aCS pCS
EXTRAMESHESDEF=("-1   1  -1    1   -1   -1   -1   -1  -1  -1" \
                "-1  -1  -1   -1    1    1    1    1  -1  -1" \
                "-1   1  -1    1    1    1    1    1  -1  -1")

# 2.1 initilization
groups=(1 2 3)
INITTEMP_DIR=$expdir/MultiTemps/

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
EVALLABELS=(CA DG SUB ERC BA35 BA36 PRC   PHC aCS pCS HIPP     EXPHIPP   ALL)
RANGES=(    1  2  3   4   5    6    "5 6" 7   8   9   "1 2 3"  "4 5 6 7" "1 2 3 4 5 6 7")
MESH_EVAL=(       CA  DG  SUB ERC BA35 BA36 PRC   PHC aCS pCS)
MESH_EVAL_RANGES=(1   2   3   4   5    6    "5 6" 7   8   9)
INITEVALDIR=$INITTEMP_DIR/evaluation

# 2.6 

# 2.7 thickness analysis
#ANGRP[0]="MRG Anterior"
ANGRP[0]="MRG"
ANGRP[1]="ERC BA35 BA36 PHC"

GRPNM[0]="merged"
GRPNM[1]="all"

# 2.8 Mesh labels

GSMESHES=(ERC BA35 BA36 PHC CS MRG)
#             BKG ERC BA35 BA36 PHC CS 
GSMESHESDEF=("-1   1  -1   -1   -1  -1" \
             "-1  -1   1   -1   -1  -1" \
             "-1  -1  -1    1   -1  -1" \
             "-1  -1  -1   -1    1  -1" \
             "-1  -1  -1   -1   -1   1" \
             "-1   1   1    1    1  -1")
GSSAMPLEPOINT=(250 300 750 500 500 2000)
#GSSAMPLEPOINT=(500 600 1500 1000 1000 4000)
GSCOMBINEMESHES=(ERC BA35 BA36 PHC CS)

GSLABEL_IDS=(BKG                ERC  BA35 BA36 PHC  CS)
GSLABEL_MRG=("0 1 2 3 4 7 8 16" "10" "11" "12" "13" "14")

GSEVALLABELS=(ERC BA35 BA36 PRC   PHC CS ALL)
GSRANGES=(    1   2    3    "2 3" 4   5  "1 2 3 4")
VTLMEVALDIR=$INITTEMP_DIR/LMMEval

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

groups_temp=(3)

# 4.6
FINALITER=6
ANTs_start_final=2
WGT_3=1
REG_LABELS_3=(${LABEL_IDS[*]})
ANTs_t_3="SyN[0.25]"
ANTs_r_3="Gauss[0.5,0]"
ANTs_i_3="60x30x10"
ANTs_x_3="Y"
ANTs_all_metrics_3="--use-all-metrics-for-convergence"

# exp1: ANTs_i_3="60x30x10"
# exp2: ANTs_i_3="200x150x100"
# exp3: ANTs_i_3="40x20x5", 8, 3

# 5. UT geodesic shoot 

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

  # 2.9 generate labels for each template mesh
  #mean_thickness

  # 2.10 subsample template meshes for geodesic shooting
  #VTMeshSubSample  

  # 2.11 geodesic shooting to generate shape prior
  #VTGeodesicShoot

  # 2.12 shape average
  #VTGShootAvg

  # 2.13 compute mode and anamation
  #VTShapeStats
  #VTMakeMovies

  # 2.14 perform landmark template fitting
  #VTLMMatching

  # 2.15 evaluate landmark shooting
  #EvalLMMatching
 
  ########################################
  # 3 insert three templates into the graph
  # 3.1 perform pairwise registration
  #Insertion

  # 3.2 measure similarity between templates and subjects
  #TempSim

  ########################################
  # 4 build graph and generate super template

  #########
  # currently these steps are manually done in MATLAB
  # 4.1 compute paths with different parameters
  # 4.2 apply warps and compute dice overlap
  # 4.3 choose the best graph and generate the pathes between templates

  # in the new pipeline, we don`t necessary need this now!
  # 1. one way is to choose the parameter so that the three templates can connect with each other
  # 2. the other way is make the graph connected, then choose the parameter so that the overlap is the best!

  # anyway, for now, just choose one parameter based on how it looks in the matlab script!

  #########

  # 4.4 propagate three templates to build super template
  #PropagateTemplates
  #PropagateTemplatesNew

  # 4.5 warp all the subjects to the super template space
  #WarpSuperTemplate
  #AverageFinalTemplate

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

  # 4.12 generate labels for each template mesh
  FinalMeanThickness

  ##################### 
  # 5.0 subsample the template mesh using Paul's Poisson sample
  #UTMeshSubSample
 
  # 5.1 geodesic shooting to generate shape prior
  #UTGeodesicShoot

  # 5.2 shape average
  #UTGShootAvg
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
  SEG=$ASHSRUNDIR/${fn}/final/${fn}_${side}_lfseg_corr_usegray_dividedCS.nii.gz

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
      mkdir -p $PWDIR/${fn_fix}_${side}

      # check if we need to reuse
      if [[ $REUSE == "Y" ]]; then
        REUSEPWDIR=$REUSEEXPDIR/pairwise/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}
        if [[ -f $REUSEPWDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz ]]; then
          ln -sf $REUSEPWDIR $OUTDIR
          continue
        fi
      fi

       # perform registration
      if [[ -f $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz ]]; then

        echo "Seg file exists."

      else

      mkdir -p $OUTDIR

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
           $ANTs_G_1 \
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

      for sub in ${KINDS[*]}; do

        rm -rf $OUTDIR/${fn_mov}_to_${fn_fix}_${side}_reslice_${sub}.nii.gz

      done

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
           similarity_sub $id_fix $side $SUBJ_TXT
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
  IDS=$(cat $3)

  # take care of template situation
  if [[ $id_fix == "template1" ]] || \
     [[ $id_fix == "template2" ]] || \
     [[ $id_fix == "template3" ]]; then
    fn_fix=$id_fix
  else
    fn_fix=$(ls $ASHSRUNDIR | grep $id_fix)
  fi
  fn_out=$SIM_DIR/${side}/${id_fix}.txt
  fn_seg_fix=$TMPDIR/${fn_fix}_${side}_seg_tmp.nii.gz
  fn_seg_mesh_fix=$TMPDIR/${fn_fix}_${side}_seg_tmp.vtk
  OVL=""

  if [ ${SIM_TYPE} == "PRCCS_seg_dice" ]; then

    c3d $DATADIR/${fn_fix}_${side}_seg.nii.gz \
      -replace 1 0 2 0 3 0 4 0 7 0 9 0 \
      -o $fn_seg_fix

  elif [ ${SIM_TYPE} = "CS_HFdist" ]; then

    c3d $DATADIR/${fn_fix}_${side}_seg.nii.gz \
      -replace 5 99 -thresh 99 99 1 0 \
      -o $fn_seg_fix
    vtklevelset $fn_seg_fix $fn_seg_mesh_fix 0.5

  fi

  # Go through all the other subjects
  for id_mov in $IDS; do

    if [[ $id_mov == "template1" ]] || \
       [[ $id_mov == "template2" ]] || \
       [[ $id_mov == "template3" ]]; then
      fn_mov=$id_mov
    else
      fn_mov=$(ls $ASHSRUNDIR | grep $id_mov)
    fi
    OVL_tmp="NA"

    if [[ $id_mov == $id_fix ]]; then

      if [ ${SIM_TYPE} = "CS_HFdist" ]; then

        OVL="$OVL 0"

      else

        OVL="$OVL 1"

      fi

    else

      if [ ${SIM_TYPE} == "PRCCS_seg_dice" ]; then

        OVL_tmp=$(c3d $PWDIR/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz \
          -replace 1 0 2 0 3 0 4 0 7 0 9 0 \
          $fn_seg_fix -label-overlap \
          | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}' )

      elif [ ${SIM_TYPE} = "CS_HFdist" ]; then
      
        c3d $PWDIR/${fn_fix}_${side}/${fn_mov}_to_${fn_fix}/${fn_mov}_to_${fn_fix}_${side}_reslice_seg.nii.gz \
          -replace 5 99 -thresh 99 99 1 0 \
          -o $TMPDIR/${fn_mov}_to_${fn_fix}_seg_tmp.nii.gz
        vtklevelset $TMPDIR/${fn_mov}_to_${fn_fix}_seg_tmp.nii.gz \
          $TMPDIR/${fn_mov}_to_${fn_fix}_seg_tmp.vtk 0.5


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
     clustering_ProbPCA_${expid} \
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
  MCIGRP=($(cat $WORKDIR/analysis_input/CNMCI.txt | sed "s/\r//g"))
  AGE=($(cat $WORKDIR/analysis_input/age.txt | sed "s/\r//g"))
  ICV=($(cat $WORKDIR/analysis_input/icv.txt | sed "s/\r//g"))

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
        str1="$id        1       0       $age    $icv	$grp"
        str2="$id        1       0       $age"
        str3="$id	1	0	$age	$grp"
      else
        str="$id	0       1       $age    $icv"
        str1="$id        0       1       $age    $icv	$grp"
        str2="$id        0       1       $age"
        str3="$id	0	1	$age	$grp"
      fi
      
      echo $str | sed "s/\r//g" >> \
        $GPSTATDIR/design_${side}_${grp}_group-age-icvcr.txt
      #echo $str | sed "s/\r//g" >> \
      #  $ALLSTATDIR/design_${side}_group-age-icvcr.txt
      #echo $str1 | sed "s/\r//g" >> \
      #  $ALLSTATDIR/design_${side}_group-age-icvcr-temp.txt
      #echo $str2 | sed "s/\r//g" >> \
      #  $ALLSTATDIR/design_${side}_group-age.txt
      #echo $str3 | sed "s/\r//g" >> \
      #  $ALLSTATDIR/design_${side}_group-age-temp.txt

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
  #echo "1 -1 0 0 0" > \
  #  $ALLSTATDIR/contrast_group-age-icvcr-temp_nc-mci.txt
  #echo "-1 1 0 0 0" > \
  #  $ALLSTATDIR/contrast_group-age-icvcr-temp_mci-nc.txt
  #echo "1 -1 0" > \
  #  $ALLSTATDIR/contrast_group-age_nc-mci.txt
  #echo "-1 1 0" > \
  #  $ALLSTATDIR/contrast_group-age_mci-nc.txt
  #echo "1 -1 0 0" > \
  #  $ALLSTATDIR/contrast_group-age-temp_nc-mci.txt
  #echo "-1 1 0 0" > \
  #  $ALLSTATDIR/contrast_group-age-temp_mci-nc.txt
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
      -r $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray_dividedCS.nii.gz \
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
    c3d $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray_dividedCS.nii.gz \
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

    echo "ID ${EVALLABELS[*]} ${MESH_EVAL[*]}" \
      > $INITEVALDIR/overlap_refseg_${side}.txt

    for ids in ${IDs[*]}; do

      id=$(ls $ASHSRUNDIR | grep $ids)
      cat $INITEVALTMPDIR/${id}_${side}_overlap.txt \
          >> $INITEVALDIR/overlap_${side}.txt

      if [[ -f $INITEVALTMPDIR/${id}_${side}_refseg_overlap.txt ]]; then
        cat $INITEVALTMPDIR/${id}_${side}_refseg_overlap.txt \
          >> $INITEVALDIR/overlap_refseg_${side}.txt
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
  do_pair $INITTEMP_DIR/group${grp}/work/template_${side}_${grp}_seg.nii.gz \
    $INITTEMP_DIR/group${grp}/work/${id}_${side}_${grp}_totemp_reslice_seg.nii.gz
  ALLOVL="$ALLOVL $FULLOVL"

  ###################################
  # save to temporary file
  echo $id $ALLOVL > $EVALTMPDIR/${id}_${side}_overlap.txt

  ###################################
  ALLOVL=""
  # Compute dice overlap between warp seg and manual seg
  FITTED=$INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_tempfit_ref.nii.gz
  REFSEG=$INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_refseg.nii.gz
  if [[ -f $REFSEG ]]; then
  # compute dice
  c3d $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_tempfit.nii.gz -replace 8 9 -o $FITTED
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
            -N "${PREFIX}_${exp}_${side}_${grp}_${GRPNM[igrp]}" \
            $0 thick_stats_sub $exp $side $grp $igrp
          #thick_stats_sub $exp $side $grp $igrp

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

    done

    meshglm $MESHPARAM \
      -g $WORK/design_${side}.txt $CWORK/contrast.txt \
      -a Thickness \
      -d 6 \
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
           -l h_vmem=2.1G,s_vmem=2G \
           -N "${PREFIX}_${grp}_${side}" \
           $0 mean_thickness_sub $grp $side
      sleep 0.1
      #mean_thickness_sub $grp $side

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
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
      if [[ $sub != "DG" ]]; then
        header="${header},${sub}_${side}"
      fi
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
          $MESHGRPDIR/template_${side}_${grp}_${EXTRAMESHES[ii]}.vtk \
          $WORKGRPDIR/template_${side}_${grp}_${sub}.nii.gz \
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

  # run matlab script to generate label
  mkdir -p $TMPDIR/${id}
  source $PKGDIR/matlab_batch.sh \
    $TMPDIR/${id}/ \
    $MATLABCODEDIR \
    compute_mean_thickness_apBA35 \
    $MESHGRPDIR/analysis/labels/thick_${side}_${grp}_MRG.vtk \
    4 \
    $INITGRP_DIR/subjID_${side}_${grp}.txt
}

##################################################
function VTMeshSubSample()
{
  PREFIX=VSS${expid}
  # Iterate over side and subject
  for grp in ${groups[*]}; do
    for side in left right; do

      # combine aCS and pCS first
      WORKGRPDIR=$INITTEMP_DIR/group${grp}/work
      ITDIR=$WORKGRPDIR/$(printf iter_%s_%02d $side $((ITER-1)) )
      c3d $ITDIR/template_${side}_${grp}_aCS.nii.gz \
        $ITDIR/template_${side}_${grp}_pCS.nii.gz \
        -add -o $ITDIR/template_${side}_${grp}_CS.nii.gz
 
      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=6.1G,s_vmem=6G \
           -N "${PREFIX}_${side}_${grp}" \
           $0 VTMeshSubSample_sub \
           $side $grp
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function VTMeshSubSample_sub()
{
  side=$1
  grp=$2
  MESHGRPDIR=$INITTEMP_DIR/group${grp}/meshwarp
  WORKGRPDIR=$INITTEMP_DIR/group${grp}/work
  LABELWARPDIR=$INITTEMP_DIR/group${grp}/labelwarp
  ITDIR=$WORKGRPDIR/$(printf iter_%s_%02d $side $((ITER-1)) )

  # Generate one mesh for all non-DG non-CS subfields
  GSDIR=$MESHGRPDIR/template/geoshoot
  mkdir -p $GSDIR
  for ((i=0;i<${#GSMESHES[*]};i++)); do

    GSLABELDEF=(${GSMESHESDEF[i]})
    echo $LABELDEF

    c3d $(for ((j=0;j<${#GSLABEL_IDS[*]};j++)); do echo "$ITDIR/template_${side}_${grp}_${GSLABEL_IDS[j]}.nii.gz -scale ${GSLABELDEF[j]}"; done) \
      -accum -add -endaccum \
      -o $GSDIR/template_${side}_${GSMESHES[i]}.nii.gz

    vtklevelset \
      $GSDIR/template_${side}_${GSMESHES[i]}.nii.gz \
      $GSDIR/template_${side}_${GSMESHES[i]}.stl \
      0

    vtklevelset \
      $GSDIR/template_${side}_${GSMESHES[i]}.nii.gz \
      $GSDIR/template_${side}_${GSMESHES[i]}.vtk \
      0

    # perform subsampling
    SAM=$GSDIR/template_${side}_${grp}_${GSMESHES[i]}_sampled.ply
    /data/picsl-build/pauly/vcg/gcc64rel/mesh_poisson_sample \
      $GSDIR/template_${side}_${GSMESHES[i]}.stl \
      $SAM ${GSSAMPLEPOINT[i]}

  done

  # combine and convert into vtk mesh
  TEMPLATE=$GSDIR/template_${side}_${grp}_MRGcombined_sampled.vtk
  NV=0
  for ((i=0;i<${#GSCOMBINEMESHES[*]};i++)); do
    SAM=$GSDIR/template_${side}_${grp}_${GSCOMBINEMESHES[i]}_sampled.ply
    NV=$(($NV+$(cat $SAM | grep 'element vertex' | awk '{print $3}')))
  done

  # Write the header of the VTK file
  echo "# vtk DataFile Version 3.0" > $TEMPLATE
  echo "vtk output" >> $TEMPLATE
  echo "ASCII" >> $TEMPLATE
  echo "DATASET POLYDATA" >> $TEMPLATE
  echo "POINTS ${NV} float" >> $TEMPLATE

  # Write the point coordinates
  for ((i=0;i<${#GSCOMBINEMESHES[*]};i++)); do
    SAM=$GSDIR/template_${side}_${grp}_${GSCOMBINEMESHES[i]}_sampled.ply
    NV=$(cat $SAM | grep 'element vertex' | awk '{print $3}')
    grep -A $NV end_header $SAM | tail -n $NV >> $TEMPLATE
  done

  # sample on probability maps
  TEMPLATESAMPLE=$GSDIR/template_${side}_${grp}_MRGcombined_sampled_labelprobs.vtk
  cp $TEMPLATE $TEMPLATESAMPLE
  for ((i=0;i<${#GSLABEL_IDS[*]};i++)); do
    mesh_image_sample $TEMPLATESAMPLE \
      $ITDIR/template_${side}_${grp}_${GSLABEL_IDS[i]}.nii.gz \
      $TEMPLATESAMPLE Prob_${GSLABEL_IDS[i]}
  done
}

function VTMeshSubSample_sub_bk()
{
  side=$1
  grp=$2
  MESHGRPDIR=$INITTEMP_DIR/group${grp}/meshwarp
  WORKGRPDIR=$INITTEMP_DIR/group${grp}/work
  LABELWARPDIR=$INITTEMP_DIR/group${grp}/labelwarp
  ITDIR=$WORKGRPDIR/$(printf iter_%s_%02d $side $((ITER-1)) )

  # Generate one mesh for all non-DG non-CS subfields
  for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

    vtklevelset \
      $ITDIR/template_${side}_${grp}_${EXTRAMESHES[i]}.nii.gz \
      $MESHGRPDIR/template_${side}_${grp}_${EXTRAMESHES[i]}.stl \
      $template_th

  done

  cp $MESHGRPDIR/template_${side}*.stl \
     $MESHGRPDIR/template/

  # perform subsampling
  SAM=$MESHGRPDIR/template/template_${side}_${grp}_MRG_sampled.ply
  /data/picsl-build/pauly/vcg/gcc64rel/mesh_poisson_sample \
    $MESHGRPDIR/template/template_${side}_${grp}_MRG.stl \
    $SAM 2000

  # comvert into vtk mesh
  TEMPLATE=$MESHGRPDIR/template/template_${side}_${grp}_MRG_sampled.vtk
  NV=$(cat $SAM | grep 'element vertex' | awk '{print $3}')

  # Write the header of the VTK file
  echo "# vtk DataFile Version 3.0" > $TEMPLATE
  echo "vtk output" >> $TEMPLATE
  echo "ASCII" >> $TEMPLATE
  echo "DATASET POLYDATA" >> $TEMPLATE
  echo "POINTS ${NV} float" >> $TEMPLATE

  # Write the point coordinates
  grep -A $NV end_header $SAM | tail -n $NV >> $TEMPLATE
}

##################################################
function VTGeodesicShoot()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=VTGS${expid}
  # Iterate over side and subject
  for grp in ${groups[*]}; do
  for side in left right; do
    GRP_IDS=$(cat $INITGRP_DIR/subjID_${side}_${grp}.txt)
    for id in $GRP_IDS; do

      id=$(ls $ASHSRUNDIR | grep $id)
      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
           -q basic.q,all.q \
           -l h_vmem=6.1G,s_vmem=6G \
           -N "${PREFIX}_${id}_${side}_${grp}" \
           $0 VTGeodesicShoot_sub \
           $id $side $grp
      sleep 0.1

    done
  done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function VTGeodesicShoot_sub()
{
  id=$1
  side=$2
  grp=$3
  MESHGRPDIR=$INITTEMP_DIR/group${grp}/meshwarp
  GSDIR=$MESHGRPDIR/template/geoshoot
  WORKGRPDIR=$INITTEMP_DIR/group${grp}/work
  LABELWARPDIR=$INITTEMP_DIR/group${grp}/labelwarp
  VTGEOSHOOTDIR=$INITTEMP_DIR/group${grp}/VTGeoShoot
  ITDIR=$WORKGRPDIR/$(printf iter_%s_%02d $side $((ITER-1)) )
  SUBJDIR=$VTGEOSHOOTDIR/$id
  mkdir -p $SUBJDIR

  #if [ -f $SUBJDIR/${id}_${side}_momentum.vtk ]; then
  #  exit
  #fi

  # warp decimated template mesh back to subject space
  # Compose the transformation between the template and the subject
  $ANTSPATH/ComposeMultiTransform 3 \
     $TMPDIR/compose.nii \
     -R $WORKGRPDIR/template_${side}_${grp}_tse.nii.gz \
     $WORKGRPDIR/${id}_${side}_${grp}_totempWarp.nii.gz \
     $WORKGRPDIR/${id}_${side}_${grp}_totempAffine.txt \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempWarp.nii \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
     $ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

  # Apply the warp to each of the meshes
  # Warp the subfield into subject space
  greedy -d 3 \
    -rf $ASHSRUNDIR/$id/mprage.nii.gz \
    -rs $GSDIR/template_${side}_${grp}_MRGcombined_sampled.vtk \
        $SUBJDIR/${id}_${side}_${grp}_MRGcombined_tempfit.vtk \
    -rs $MESHGRPDIR/template/template_${side}_${grp}_MRG.vtk \
        $SUBJDIR/${id}_${side}_${grp}_MRG_tempfit.vtk \
    -r $TMPDIR/compose.nii

  # use procrustes to align the two meshes
  vtkprocrustes $SUBJDIR/${id}_${side}_${grp}_MRGcombined_tempfit.vtk \
    $GSDIR/template_${side}_${grp}_MRGcombined_sampled.vtk \
    $SUBJDIR/${id}_to_template_${side}_${grp}_meshrigid.mat
    #$SUBJDIR/template_to_${id}_${side}_${grp}_rigid.mat
  c3d_affine_tool $SUBJDIR/${id}_to_template_${side}_${grp}_meshrigid.mat \
    -inv -o $SUBJDIR/template_to_${id}_${side}_${grp}_meshrigid.mat
  greedy -d 3 \
    -rf $ASHSRUNDIR/$id/mprage.nii.gz \
    -rs $GSDIR/template_${side}_${grp}_MRGcombined_sampled.vtk \
        $SUBJDIR/${id}_${side}_${grp}_MRGcombined_tempfit_rigid.vtk \
    -rs $MESHGRPDIR/template/template_${side}_${grp}_MRG.vtk \
        $SUBJDIR/${id}_${side}_${grp}_MRG_tempfit_rigid.vtk \
    -r $SUBJDIR/template_to_${id}_${side}_${grp}_meshrigid.mat

  greedy -d 3 \
    -rf /home/longxie/ASHS_PHC/SuperResolution/ASHSexp/exp003/headtailatlas/final/template/template.nii.gz \
    -rs $SUBJDIR/${id}_${side}_${grp}_MRGcombined_tempfit.vtk \
        $SUBJDIR/${id}_${side}_${grp}_totemp_MRGcombined_sampled.vtk \
    -rs $SUBJDIR/${id}_${side}_${grp}_MRG_tempfit.vtk \
        $SUBJDIR/${id}_${side}_${grp}_totemp_MRG_rigid.vtk \
    -r $SUBJDIR/${id}_to_template_${side}_${grp}_meshrigid.mat

  # perform geodesic shooting
  #/home/cwang/pcmrep/PointSetGeodesicShooting_CUDA_c00/lmshoot_cuda -d 3 \
  #  -m $GSDIR/template_${grp}_MRGcombined_sampled.vtk \
  #     $SUBJDIR/${id}_${side}_${grp}_MRGcombined_tempfit_rigid.vtk \
  #  -o $SUBJDIR/${id}_${side}_${grp}_momentum.vtk \
  #  -s 2.0 -l 5000 \
  #  -n 40 -i 200 0 \
  #  -r 1

  if [ ! -f $SUBJDIR/${id}_${side}_${grp}_momentum.vtk ]; then
    lmshoot -d 3 \
      -m $GSDIR/template_${side}_${grp}_MRGcombined_sampled.vtk \
         $SUBJDIR/${id}_${side}_${grp}_totemp_MRGcombined_sampled.vtk \
      -o $SUBJDIR/${id}_${side}_${grp}_momentum.vtk \
      -s 2.0 -l 5000 \
      -n 40 -i 200 0 \
      -f -r 1
  fi

  #if [ ! -f $SUBJDIR/template_to_${id}_${side}_${grp}_GSShoot_meshwarp.nii.gz ]; then
    $LMTOWARPDIR/lmtowarp -d 3 -n 40 -s 2.0 \
      -r $ASHSRUNDIR/$id/mprage.nii.gz \
      -m $SUBJDIR/${id}_${side}_${grp}_momentum.vtk \
      -o $SUBJDIR/template_to_${id}_${side}_${grp}_GSShoot_meshwarp.nii.gz \
      -oinv $SUBJDIR/template_to_${id}_${side}_${grp}_GSShoot_invmeshwarp.nii.gz
  #fi

  greedy -d 3 \
    -rf $ASHSRUNDIR/$id/mprage.nii.gz \
    -rs $MESHGRPDIR/template/template_${side}_${grp}_MRG.vtk \
        $SUBJDIR/${id}_${side}_${grp}_MRG_tempfit_GSShoot.vtk \
    -r $SUBJDIR/template_to_${id}_${side}_${grp}_GSShoot_meshwarp.nii.gz \
       $SUBJDIR/template_to_${id}_${side}_${grp}_meshrigid.mat

  # inverse warp field
  #greedy -d 3 \
  #  -iw $SUBJDIR/template_to_${id}_${side}_${grp}_GSShoot_meshwarp.nii.gz \
  #      $SUBJDIR/template_to_${id}_${side}_${grp}_GSShoot_invmeshwarp.nii.gz \
  #  -invexp 4

  # warp template segmentation to subject space
  TRANS=""
   for ((i=0;i<${#GSLABEL_IDS[*]};i++)); do
     TRANS="$TRANS -rm $ITDIR/template_${side}_${grp}_${GSLABEL_IDS[i]}.nii.gz $SUBJDIR/template_to_${id}_${side}_${grp}_GSShoot_${GSLABEL_IDS[i]}.nii.gz"
   done

  greedy -d 3 \
    -rf $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray_dividedCS.nii.gz \
    $TRANS \
    -r $SUBJDIR/${id}_to_template_${side}_${grp}_meshrigid.mat \
       $SUBJDIR/template_to_${id}_${side}_${grp}_GSShoot_invmeshwarp.nii.gz

  # Vote in original subject space to get the tempfit segmentation
  c3d $(for sub in ${GSLABEL_IDS[*]}; do echo $SUBJDIR/template_to_${id}_${side}_${grp}_GSShoot_${sub}.nii.gz; done) -vote -type ushort \
    -o $SUBJDIR/template_to_${id}_${side}_${grp}_GSShoot_seg.nii.gz
}

##################################################
function VTGShootAvg()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=VTGSA${expid}
  # Iterate over side and subject
  for grp in ${groups[*]}; do
    for side in left right; do

      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=20.1G,s_vmem=20G \
           -N "${PREFIX}_${grp}_${side}" $0 \
           VTGShootAvg_sub $side $grp
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function VTGShootAvg_sub()
{
  side=$1
  grp=$2
  MESHGRPDIR=$INITTEMP_DIR/group${grp}/meshwarp
  GSDIR=$MESHGRPDIR/template/geoshoot
  WORKGRPDIR=$INITTEMP_DIR/group${grp}/work
  LABELWARPDIR=$INITTEMP_DIR/group${grp}/labelwarp
  VTGEOSHOOTDIR=$INITTEMP_DIR/group${grp}/VTGeoShoot
  mkdir -p $VTGEOSHOOTDIR/shavg

  # Average the momentum maps
  avgmesharr \
    $VTGEOSHOOTDIR/*/*_${side}_${grp}_momentum.vtk \
    InitialMomentum \
    $GSDIR/template_${side}_${grp}_MRGcombined_sampled.vtk \
    $VTGEOSHOOTDIR/shavg/template_${side}_${grp}_MRGcombined_average_momenta.vtk

  # Perform the shooting and generate warp
  lmtowarp -d 3 -n 40 \
    -r $WORKGRPDIR/template_${side}_${grp}_seg.nii.gz \
    -m $VTGEOSHOOTDIR/shavg/template_${side}_${grp}_MRGcombined_average_momenta.vtk \
    -o $VTGEOSHOOTDIR/shavg/template_${side}_${grp}_MRGcombined_average_momenta_warp.nii.gz \
    -s 2.0

  # Apply the warp to the template
  greedy -d 3 \
    -rf $WORKGRPDIR/template_${side}_${grp}_seg.nii.gz \
    -rs $MESHGRPDIR/template_${side}_${grp}_MRG.vtk \
        $VTGEOSHOOTDIR/shavg/template_${side}_${grp}_MRG_Gshoot.vtk \
    -r $VTGEOSHOOTDIR/shavg/template_${side}_${grp}_MRGcombined_average_momenta_warp.nii.gz
}

##################################################
function VTShapeStats()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=VTSS${expid}
  # Iterate over side and subject
  for grp in ${groups[*]}; do
    for side in left right; do

      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${grp}_${side}" $0 \
           VTShapeStats_sub $side $grp
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function VTShapeStats_sub()
{
  side=$1
  grp=$2
  MESHGRPDIR=$INITTEMP_DIR/group${grp}/meshwarp
  GSDIR=$MESHGRPDIR/template/geoshoot
  WORKGRPDIR=$INITTEMP_DIR/group${grp}/work
  LABELWARPDIR=$INITTEMP_DIR/group${grp}/labelwarp
  VTGEOSHOOTDIR=$INITTEMP_DIR/group${grp}/VTGeoShoot
  PCADIR=$VTGEOSHOOTDIR/shavg/pca/${side}
  rm -rf $PCADIR
  mkdir -p $PCADIR
  GRP_IDS=$(cat $INITGRP_DIR/subjID_${side}_${grp}.txt)

  # Generate the momentum data in a format that is readable by R script
  idx=1
  echo "ID" > $PCADIR/IDs_${side}_${grp}.csv
  for id in $GRP_IDS; do
    id=$(ls $ASHSRUNDIR | grep $id)
    local MOMENTS=$VTGEOSHOOTDIR/$id/${id}_${side}_${grp}_momentum.vtk
    echo $(echo $idx) $(dumpmeshattr $MOMENTS InitialMomentum) >> $PCADIR/initial_momenta_${side}_${grp}.txt
    echo $idx >> $PCADIR/IDs_${side}_${grp}.csv
    idx=$((idx+1))
  done

  # Generate the label probability data in the template space in a format that is readable by R script
  local LP=$GSDIR/template_${side}_${grp}_MRGcombined_sampled_labelprobs.vtk
  for ((i=0;i<${#GSLABEL_IDS[*]};i++)); do
    echo Prob_${GSLABEL_IDS[i]} $(dumpmeshattr $LP Prob_${GSLABEL_IDS[i]}) >> $PCADIR/label_probs_${side}_${grp}.txt
  done

  # Generate an array of mesh coordinates for R
  dumpmeshpoints $VTGEOSHOOTDIR/shavg/template_${side}_${grp}_MRGcombined_average_momenta.vtk | tail -n +3 > $PCADIR/template_points_${side}_${grp}.txt

  # Run the R statistics notebook
  #Rscript run_pca.R \
  #  $PCADIR/IDs_${side}_${grp}.csv \
  #  $PCADIR/initial_momenta_${side}_${grp}.txt \
  #  $PCADIR/template_points_${side}_${grp}.txt \
  #  $PCADIR/label_probs_${side}_${grp}.txt \
  #  $PCADIR
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('/data/picsl/longxie/ASHS/thickness_newlabel/matlabcode');
    runPCA('$PCADIR/initial_momenta_${side}_${grp}.txt', ...
           '$PCADIR/label_probs_${side}_${grp}.txt', ...
           '$PCADIR/template_points_${side}_${grp}.txt', ...
           '$PCADIR', '$PCADIR/template_allinfo_${side}_${grp}.vtk');
MATCODE
}

##################################################
function VTMakeMovies()
{
  PREFIX=VTMM${expid}
  # Iterate over side and subject
  for grp in ${groups[*]}; do
    for side in left right; do
      for what in mode1 mode2 mode3 mode4 mode5; do
      
      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${grp}_${side}_${what}" $0 \
           VTMakeMovies_sub $side $grp $what
      sleep 0.1

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function VTMakeMovies_sub()
{
  side=$1
  grp=$2
  what=$3
  MESHGRPDIR=$INITTEMP_DIR/group${grp}/meshwarp
  WORKGRPDIR=$INITTEMP_DIR/group${grp}/work
  LABELWARPDIR=$INITTEMP_DIR/group${grp}/labelwarp
  VTGEOSHOOTDIR=$INITTEMP_DIR/group${grp}/VTGeoShoot
  PCADIR=$VTGEOSHOOTDIR/shavg/pca/${side}
  GRP_IDS=$(cat $INITGRP_DIR/subjID_${side}_${grp}.txt)
  MOVIEDIR=$PCADIR/movie_${what}
  mkdir -p $MOVIEDIR

  # Template
  local TEMPMESH=$MESHGRPDIR/analysis/labels/thick_${side}_${grp}_MRG.vtk
  
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

##################################################
function VTLMMatching()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=VTLM${expid}
  # Iterate over side and subject
  for grp in ${groups[*]}; do
  for side in left right; do
    GRP_IDS=$(cat $INITGRP_DIR/subjID_${side}_${grp}.txt)
    for id in $GRP_IDS; do

      id=$(ls $ASHSRUNDIR | grep $id)
      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
           -q basic.q,all.q \
           -l h_vmem=6.1G,s_vmem=6G \
           -N "${PREFIX}_${id}_${side}_${grp}" \
           $0 VTLMMatching_sub \
           $id $side $grp
      sleep 0.1

    done
  done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function VTLMMatching_sub()
{
  id=$1
  fn=$(ls $ASHSRUNDIR | grep $id)
  side=$2
  grp=$3
  MESHGRPDIR=$INITTEMP_DIR/group${grp}/meshwarp
  GSDIR=$MESHGRPDIR/template/geoshoot
  WORKGRPDIR=$INITTEMP_DIR/group${grp}/work
  LABELWARPDIR=$INITTEMP_DIR/group${grp}/labelwarp
  VTGEOSHOOTDIR=$INITTEMP_DIR/group${grp}/VTGeoShoot
  PCADIR=$VTGEOSHOOTDIR/shavg/pca/${side}
  VTLMMATCHINGDIR=$INITTEMP_DIR/group${grp}/LMMatching
  ASHSSUBJDIR=$ASHSRUNDIR/$fn
  ITDIR=$WORKGRPDIR/$(printf iter_%s_%02d $side $((ITER-1)) )
  SUBJDIR=$VTLMMATCHINGDIR/$id
  mkdir -p $SUBJDIR

  # normalize posterior maps
  c3d $ASHSSUBJDIR/bootstrap/fusion/posterior_corr_usegray_${side}_*.nii.gz \
    -accum -add -endaccum \
    -o $SUBJDIR/posterior_corr_usegray_${side}_sum.nii.gz 

  # generate label probability maps
  for ((i=0;i<${#GSLABEL_IDS[*]};i++)); do

    LABELS=(${GSLABEL_MRG[i]})
    if [[ ${#LABELS[*]} == 1 ]]; then
      c3d $SUBJDIR/posterior_corr_usegray_${side}_sum.nii.gz \
        $ASHSSUBJDIR/bootstrap/fusion/posterior_corr_usegray_${side}_$(printf %03d ${LABELS[0]}}).nii.gz \
        -divide \
        -o $SUBJDIR/posterior_corr_usegray_${side}_${GSLABEL_IDS[i]}.nii.gz
    else
      c3d $(for k in ${GSLABEL_MRG[i]}; do echo "$ASHSSUBJDIR/bootstrap/fusion/posterior_corr_usegray_${side}_$(printf %03d $k).nii.gz "; done) \
        -accum -add -endaccum -popas SUM \
        $SUBJDIR/posterior_corr_usegray_${side}_sum.nii.gz \
        -push SUM -divide \
        -o $SUBJDIR/posterior_corr_usegray_${side}_${GSLABEL_IDS[i]}.nii.gz
    fi

    if [[ ${GSLABEL_IDS[i]} == 'BKG' ]]; then
      c3d $SUBJDIR/posterior_corr_usegray_${side}_${GSLABEL_IDS[i]}.nii.gz -replace nan 1 inf 1 -o $SUBJDIR/posterior_corr_usegray_${side}_${GSLABEL_IDS[i]}.nii.gz
    else
      c3d $SUBJDIR/posterior_corr_usegray_${side}_${GSLABEL_IDS[i]}.nii.gz -replace nan 0 inf 0 -o $SUBJDIR/posterior_corr_usegray_${side}_${GSLABEL_IDS[i]}.nii.gz
    fi

  done

  # generate MRG surface
  LABELDEF=(${GSMESHESDEF[4]})
  c3d $(for ((j=0;j<${#GSLABEL_IDS[*]};j++)); do echo "$SUBJDIR/posterior_corr_usegray_${side}_${GSLABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
      -accum -add -endaccum \
      -o $SUBJDIR/posterior_corr_usegray_${side}_MRG.nii.gz
  vtklevelset $SUBJDIR/posterior_corr_usegray_${side}_MRG.nii.gz \
    $SUBJDIR/posterior_corr_usegray_${side}_MRG.vtk 0.00001

  # warp the template mesh to subject space and rigidly warp back to template space? or just get the transformation is available.
  # we can get this in geodesic shooting step!
  c3d_affine_tool $VTGEOSHOOTDIR/$id/template_to_${id}_${side}_${grp}_rigid.mat -inv \
    -o $SUBJDIR/template_to_${id}_${side}_${grp}_rigid.mat

  # link the automatic seg and ground truth seg if available
  c3d $ASHSSUBJDIR/tse_native_chunk_${side}.nii.gz \
    $ASHSSUBJDIR/final/${fn}_${side}_lfseg_corr_usegray.nii.gz \
    -int 0 -reslice-identity \
    -o $SUBJDIR/autoseg_${side}.nii.gz
  if [[ -d $ASHSSUBJDIR/refseg/ ]]; then
    c3d $ASHSSUBJDIR/tse_native_chunk_${side}.nii.gz \
      $ASHSSUBJDIR/refseg/refseg_${side}.nii.gz \
      -int 0 -reslice-identity \
      -o $SUBJDIR/refseg_${side}.nii.gz
  fi

  # perform template landmark matching
  if [ ! -f $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching.vtk ]; then
  $LMMATCHINGDIR/lmtpmatch -d 3 \
    -m $PCADIR/template_allinfo_${side}_${grp}.vtk \
    -a $VTGEOSHOOTDIR/$id/template_to_${id}_${side}_${grp}_meshrigid.mat \
    -p $(for ((j=0;j<${#GSLABEL_IDS[*]};j++)); do echo "$SUBJDIR/posterior_corr_usegray_${side}_${GSLABEL_IDS[j]}.nii.gz"; done) \
    -o $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching.vtk \
    -im $VTGEOSHOOTDIR/$id/${id}_${side}_${grp}_momentum.vtk \
    -s 2.0 -l 5000 -t 50000 -n 40 -i 30 0 -npc 10
  fi

  #if [ ! -f $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching_meshwarp.nii.gz ]; then
  $LMTOWARPDIR/lmtowarp -d 3 -n 40 -s 2.0 \
    -r $ASHSRUNDIR/$id/mprage.nii.gz \
    -m $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching.vtk \
    -o $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching_meshwarp.nii.gz \
    -oinv $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching_invmeshwarp.nii.gz
  #fi

  greedy -d 3 \
    -rf $ASHSRUNDIR/$id/mprage.nii.gz \
    -rs $MESHGRPDIR/template/template_${side}_${grp}_MRG.vtk \
        $SUBJDIR/${id}_${side}_${grp}_MRG_LMMatching.vtk \
    -r $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching_meshwarp.nii.gz \
       $VTGEOSHOOTDIR/$id/template_to_${id}_${side}_${grp}_meshrigid.mat

  # inverse warp field
  #greedy -d 3 \
  #  -iw $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching_meshwarp.nii.gz \
  #      $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching_invmeshwarp.nii.gz \
  #  -invexp 4

  # warp template segmentation to subject space
  TRANS=""
  for ((i=0;i<${#GSLABEL_IDS[*]};i++)); do
    TRANS="$TRANS -rm $ITDIR/template_${side}_${grp}_${GSLABEL_IDS[i]}.nii.gz $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching_${GSLABEL_IDS[i]}.nii.gz"
  done

  greedy -d 3 \
    -rf $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray_dividedCS.nii.gz \
    $TRANS \
    -r $VTGEOSHOOTDIR/$id/${id}_to_template_${side}_${grp}_meshrigid.mat \
       $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching_invmeshwarp.nii.gz

  # Vote in original subject space to get the tempfit segmentation
  c3d $(for sub in ${GSLABEL_IDS[*]}; do echo $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching_${sub}.nii.gz; done) -vote -type ushort \
    -o $SUBJDIR/template_to_${id}_${side}_${grp}_LMMatching_seg.nii.gz
}

##############################################################################
function EvalLMMatching()
{
  VTLMEVALTMPDIR=$VTLMEVALDIR/tmp
  rm -rf $VTLMEVALTMPDIR
  mkdir -p $VTLMEVALTMPDIR

  PREFIX=EVVTLM${expid}
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
          $0 EvalLMMatching_sub $id $side $grp $VTLMEVALTMPDIR
        sleep 0.1

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # Combine scores
  header=""
  for sub in ${GSEVALLABELS[*]}; do
    header="${header},${sub}"
  done

  echo "ID,SIDE,GRP${header}" \
    > $VTLMEVALDIR/overlap_ashs_refseg.csv

  echo "ID,SIDE,GRP${header}" \
    > $VTLMEVALDIR/overlap_volreg_ashs.csv

  echo "ID,SIDE,GRP${header}" \
    > $VTLMEVALDIR/overlap_GSShoot_ashs.csv

  echo "ID,SIDE,GRP${header}" \
    > $VTLMEVALDIR/overlap_LMMatching_ashs.csv

  echo "ID,SIDE,GRP${header}" \
    > $VTLMEVALDIR/overlap_volreg_refseg.csv

  echo "ID,SIDE,GRP${header}" \
    > $VTLMEVALDIR/overlap_GSShoot_refseg.csv

  echo "ID,SIDE,GRP${header}" \
    > $VTLMEVALDIR/overlap_LMMatching_refseg.csv

  IDs=($(cat $SUBJ_TXT))
  for side in left right; do
    for ids in ${IDs[*]}; do

      id=$(ls $ASHSRUNDIR | grep $ids)

      # ashs vs refseg
      if [[ -f $VTLMEVALTMPDIR/${id}_${side}_ashs_overlap_refseg.csv ]]; then
        cat $VTLMEVALTMPDIR/${id}_${side}_ashs_overlap_refseg.csv \
          >> $VTLMEVALDIR/overlap_ashs_refseg.csv
      fi

      # volume registration
      cat $VTLMEVALTMPDIR/${id}_${side}_volreg_overlap_ashs.csv \
          >> $VTLMEVALDIR/overlap_volreg_ashs.csv
      if [[ -f $VTLMEVALTMPDIR/${id}_${side}_volreg_overlap_refseg.csv ]]; then
        cat $VTLMEVALTMPDIR/${id}_${side}_volreg_overlap_refseg.csv \
          >> $VTLMEVALDIR/overlap_volreg_refseg.csv
      fi

      # gs shooting
      cat $VTLMEVALTMPDIR/${id}_${side}_GSShoot_overlap_ashs.csv \
          >> $VTLMEVALDIR/overlap_GSShoot_ashs.csv
      if [[ -f $VTLMEVALTMPDIR/${id}_${side}_GSShoot_overlap_refseg.csv ]]; then
        cat $VTLMEVALTMPDIR/${id}_${side}_GSShoot_overlap_refseg.csv \
          >> $VTLMEVALDIR/overlap_GSShoot_refseg.csv
      fi

      # template lm matching
      cat $VTLMEVALTMPDIR/${id}_${side}_LMMatching_overlap_ashs.csv \
          >> $VTLMEVALDIR/overlap_LMMatching_ashs.csv
      if [[ -f $VTLMEVALTMPDIR/${id}_${side}_LMMatching_overlap_refseg.csv ]]; then
        cat $VTLMEVALTMPDIR/${id}_${side}_LMMatching_overlap_refseg.csv \
          >> $VTLMEVALDIR/overlap_LMMatching_refseg.csv
      fi

    done
  done
}


function EvalLMMatching_sub()
{
  id=$1
  side=$2
  grp=$3
  EVALTMPDIR=$4
  MESHGRPDIR=$INITTEMP_DIR/group${grp}/meshwarp
  GSDIR=$MESHGRPDIR/template/geoshoot
  WORKGRPDIR=$INITTEMP_DIR/group${grp}/work
  LABELWARPDIR=$INITTEMP_DIR/group${grp}/labelwarp
  VTGEOSHOOTDIR=$INITTEMP_DIR/group${grp}/VTGeoShoot
  PCADIR=$VTGEOSHOOTDIR/shavg/pca/${side}
  VTLMMATCHINGDIR=$INITTEMP_DIR/group${grp}/LMMatching
  #ALLOVL=""

  ###################################
  # ashs vs. refseg
  ASHSSEG=$TMPDIR/${id}_${side}_${grp}_seg_ashs.nii.gz
  c3d $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_ashs.nii.gz \
    -replace 9 8 -o $ASHSSEG

  REFSEG=$TMPDIR/${id}_${side}_${grp}_seg_refseg.nii.gz
  if [ -f $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_refseg.nii.gz ]; then
    c3d $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_refseg.nii.gz \
      -replace 9 8 -o $REFSEG

    do_pair_GS $ASHSSEG $REFSEG
    echo "${id},${side},${grp}${FULLOVL}" > \
      $EVALTMPDIR/${id}_${side}_ashs_overlap_refseg.csv
  fi

  ###################################
  # direct registration
  FITTED=$TMPDIR/template_to_${id}_${side}_${grp}_seg_tempfit.nii.gz
  c3d $INITTEMP_DIR/group${grp}/labelwarp/${id}_${side}_${grp}_seg_tempfit.nii.gz \
    -replace 9 8 -o $FITTED
  do_pair_GS $ASHSSEG $FITTED
  echo "${id},${side},${grp}${FULLOVL}" > \
    $EVALTMPDIR/${id}_${side}_volreg_overlap_ashs.csv

  if [ -f $REFSEG ]; then
    do_pair_GS $REFSEG $FITTED
    echo "${id},${side},${grp}${FULLOVL}" > \
      $EVALTMPDIR/${id}_${side}_volreg_overlap_refseg.csv
  fi

  ###################################
  # GS shoot
  FITTED=$VTGEOSHOOTDIR/$id/template_to_${id}_${side}_${grp}_GSShoot_seg.nii.gz
  do_pair_GS $ASHSSEG $FITTED
  echo "${id},${side},${grp}${FULLOVL}" > \
    $EVALTMPDIR/${id}_${side}_GSShoot_overlap_ashs.csv

  if [ -f $REFSEG ]; then
    do_pair_GS $REFSEG $FITTED
    echo "${id},${side},${grp}${FULLOVL}" > \
      $EVALTMPDIR/${id}_${side}_GSShoot_overlap_refseg.csv
  fi

  ###################################
  # LM Matching
  FITTED=$VTLMMATCHINGDIR/$id/template_to_${id}_${side}_${grp}_LMMatching_seg.nii.gz
  do_pair_GS $ASHSSEG $FITTED
  echo "${id},${side},${grp}${FULLOVL}" > \
    $EVALTMPDIR/${id}_${side}_LMMatching_overlap_ashs.csv

  if [ -f $REFSEG ]; then
    do_pair_GS $REFSEG $FITTED
    echo "${id},${side},${grp}${FULLOVL}" > \
      $EVALTMPDIR/${id}_${side}_LMMatching_overlap_refseg.csv
  fi
}

function do_pair_GS()
{
  # Get a pair of segmentations
  seg_a=$1
  seg_b=$2

  # out dice file
  #out_dice_file=$3

  # Iterate over all relevant labels
  FULLOVL=""
  for ((i=0; i<${#GSEVALLABELS[*]}; i++)); do

    # Do the analysis on full-size meshes
    REPRULE=$(for lab in ${GSRANGES[i]}; do echo $lab 99; done)

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

    # perform registration between templates
    for grp_src in ${groups[*]}; do
      for grp_tg in ${groups[*]}; do

        if [[ $grp_src -eq $grp_tg ]]; then
          continue
        fi
        
        qsub -cwd -o $DUMPDIR -j y \
          -N "${PREFIX}_${grp_src}_${grp_tg}_${side}" $0 \
          Insertion_sub \
          template${grp_src} template${grp_tg} $side

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

  if [[ $id_mov != $id_fix ]]; then

    if [[ -d $OUTDIR ]]; then

      echo "Seg file exists."
      rm -rf $OUTDIR

    fi

    mkdir -p $OUTDIR

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
  IDS=$(cat $SUBJ_WITHTEMP_TXT)

  PREFIX=TS${expid}
  for side in left right; do
    
    rm -rf ${SIM_DIR}/${side}
    mkdir -p ${SIM_DIR}/${side}
    
    for id_fix in $IDS; do

      qsub -cwd -o $DUMPDIR -j y \
           -l h_vmem=2.1G,s_vmem=2G \
           -N "${PREFIX}_${id_fix}_${side}" $0 \
           similarity_sub $id_fix $side $SUBJ_WITHTEMP_TXT
      sleep 0.1

    done
  done

  # wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  for side in left right; do

    # Concadenate the results
    #echo '' > ${SIM_DIR}/sim_PRC_seg_msq_${side}.txt
    rm -rf ${SIM_DIR}/sim_${SIM_TYPE}_withtemp_${side}.txt

    for id_fix in $IDS; do

        cat ${SIM_DIR}/${side}/${id_fix}.txt \
            >> ${SIM_DIR}/sim_${SIM_TYPE}_withtemp_${side}.txt

    done

    rm -rf ${SIM_DIR}/${side}

  done
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
  idx=$((grp_fix-1))
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
  # initial template to target template
  TRANSDIR=$FINALTEMPDIR/template${grp_fix}/tempwarps
##########################################################
  TRANSITERDIR=$TRANSDIR/${grp_mov}_to_${grp_fix}_${side}/${i}_template${grp_mov}/
  mkdir -p $TRANSITERDIR
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $DATADIR/template${grp_mov}_${side}_${sub}.nii.gz; done) \
    -vote -type ushort -as SEG \
    $DATADIR/template${grp_fix}_${side}_tse.nii.gz \
    -push SEG -interp 0 -reslice-identity \
    -o $TRANSITERDIR/template_${side}_${grp_mov}.nii.gz
  for ((ii=0;ii<${#EXTRAMESHES[*]};ii++)); do
    ln -sf $INITTEMP_DIR/group${grp_mov}/meshwarp/analysis/labels/thick_${side}_${grp_mov}_${EXTRAMESHES[ii]}.vtk $TRANSITERDIR/${i}_template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.vtk
  done
#########################################################
  warps=""   
  invwarps=""
  mkdir -p $TRANSDIR
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

########################################################
    # compose transform
    TRANSITERDIR=$TRANSDIR/${grp_mov}_to_${grp_fix}_${side}/${i}_${side}_${id_from}_to_${id_to}
    TRANSITERTMPDIR=$TRANSITERDIR/tmp
    mkdir -p $TRANSITERTMPDIR
    $ANTSPATH/ComposeMultiTransform 3 \
      $TRANSITERTMPDIR/template${grp_mov}_${side}_InitWarp.nii.gz \
      -R $DATADIR/template${grp_fix}_${side}_tse.nii.gz \
      $warps

    # apply init warps
    for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        $DATADIR/template${grp_mov}_${side}_${sub}.nii.gz \
        $TRANSITERTMPDIR/template${grp_mov}_${side}_inittotemp_reslice_${sub}.nii.gz \
        -R $DATADIR/template${grp_fix}_${side}_tse.nii.gz \
        $TRANSITERTMPDIR/template${grp_mov}_${side}_InitWarp.nii.gz

    done

    # Create the segmentation for the template
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $TRANSITERTMPDIR/template${grp_mov}_${side}_inittotemp_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $TRANSITERDIR/template${grp_mov}_${side}_inittotemp_reslice_seg.nii.gz

    # generate a surface mesh
    for ((ii=0;ii<${#EXTRAMESHES[*]};ii++)); do
      LABELDEF=(${EXTRAMESHESDEF[ii]})
      c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$TRANSITERTMPDIR/template${grp_mov}_${side}_inittotemp_reslice_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
        -accum -add -endaccum \
        -o $TRANSITERTMPDIR/template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.nii.gz

      vtklevelset \
        $TRANSITERTMPDIR/template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.nii.gz \
        $TRANSITERDIR/${i}_template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.vtk \
        $template_th

      # generate label for the mesh
      MESHES=""
      for ((jj=0;jj<${#LABEL_IDS[*]};jj++)); do
        if [[ ${LABELDEF[$jj]} == "1" ]]; then
          sub=${LABEL_IDS[jj]}
          mesh_image_sample \
            $TRANSITERDIR/${i}_template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.vtk \
            $TRANSITERTMPDIR/template${grp_mov}_${side}_inittotemp_reslice_${LABEL_IDS[jj]}.nii.gz \
            $TRANSITERTMPDIR/template_${side}_${grp_mov}_${EXTRAMESHES[ii]}_${sub}.vtk \
            PROB
          MESHES="$MESHES $TRANSITERTMPDIR/template_${side}_${grp_mov}_${EXTRAMESHES[ii]}_${sub}.vtk"
        fi
      done

      # merge prob arrays
      mesh_merge_arrays \
        -r $TRANSITERDIR/${i}_template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.vtk \
        $TRANSITERDIR/${i}_template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.vtk \
        PROB $MESHES

      # run matlab script to generate label
      $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
        addpath('/data/picsl/longxie/ASHS/thickness_newlabel/matlabcode');
        compute_mesh_label('$TRANSITERDIR/${i}_template_${side}_${grp_mov}_${EXTRAMESHES[$ii]}.vtk');
MATCODE

    done

    rm -rf $TRANSITERTMPDIR
##########################################################

  done

##########################################################
  TRANSITERDIR=$TRANSDIR/${grp_mov}_to_${grp_fix}_${side}/target_template${grp_fix}/
  mkdir -p $TRANSITERDIR
  ln -sf $DATADIR/template${grp_fix}_${side}_seg.nii.gz \
         $TRANSITERDIR/template_${side}_${grp_fix}.nii.gz
  for ((ii=0;ii<${#EXTRAMESHES[*]};ii++)); do
    ln -sf $INITTEMP_DIR/group${grp_fix}/meshwarp/analysis/labels/thick_${side}_${grp_fix}_${EXTRAMESHES[ii]}.vtk $TRANSITERDIR/${i}_template_${side}_${grp_fix}_${EXTRAMESHES[ii]}.vtk
  done

  # link the direct registration result from pairwise registration
  TRANSITERDIR=$TRANSDIR/${grp_mov}_to_${grp_fix}_${side}/direct/
  TRANSITERTMPDIR=$TRANSITERDIR/tmp
  mkdir -p $TRANSITERTMPDIR
  for ((ii=0;ii<${#EXTRAMESHES[*]};ii++)); do
    LABELDEF=(${EXTRAMESHESDEF[ii]})
    c3d $(for ((j=0;j<${#LABEL_IDS[*]};j++)); do echo "$PWDIR/template${grp_fix}_${side}/template${grp_mov}_to_template${grp_fix}/template${grp_mov}_to_template${grp_fix}_${side}_reslice_${LABEL_IDS[j]}.nii.gz -scale ${LABELDEF[j]}"; done) \
      -accum -add -endaccum \
      -o $TRANSITERTMPDIR/template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.nii.gz

    vtklevelset \
      $TRANSITERTMPDIR/template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.nii.gz \
      $TRANSITERDIR/direct_template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.vtk \
      $template_th

    # generate label for the mesh
    MESHES=""
    for ((jj=0;jj<${#LABEL_IDS[*]};jj++)); do
      if [[ ${LABELDEF[$jj]} == "1" ]]; then
        sub=${LABEL_IDS[jj]}
        mesh_image_sample \
          $TRANSITERDIR/direct_template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.vtk \
          $PWDIR/template${grp_fix}_${side}/template${grp_mov}_to_template${grp_fix}/template${grp_mov}_to_template${grp_fix}_${side}_reslice_${LABEL_IDS[jj]}.nii.gz \
          $TRANSITERTMPDIR/template_${side}_${grp_mov}_${EXTRAMESHES[ii]}_${sub}.vtk \
          PROB
        MESHES="$MESHES $TRANSITERTMPDIR/template_${side}_${grp_mov}_${EXTRAMESHES[ii]}_${sub}.vtk"
      fi
    done

    # merge prob arrays
    mesh_merge_arrays \
      -r $TRANSITERDIR/direct_template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.vtk \
      $TRANSITERDIR/direct_template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.vtk \
      PROB $MESHES

      # run matlab script to generate label
      $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
        addpath('/data/picsl/longxie/ASHS/thickness_newlabel/matlabcode');
        compute_mesh_label('$TRANSITERDIR/direct_template_${side}_${grp_mov}_${EXTRAMESHES[$ii]}.vtk');
MATCODE

   done
 
   rm -rf $TRANSITERTMPDIR 
##########################################################

  # compose transform
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
function PropagateTemplatesNew()
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
        -l h_vmem=8.1G,s_vmem=8G \
        -N "${PREFIX}_${side}_${grp_fix}_${grp_mov}" $0 \
        PropagateTemplatesNew_sub $side $grp_fix $grp_mov
      sleep 0.1

      done
    done
  done

  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function PropagateTemplatesNew_sub()
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
  idx=$((grp_fix-1))
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
  # initial template to target template
  TRANSDIR=$FINALTEMPDIR/template${grp_fix}/tempwarps
  TRANSITERDIR=$TRANSDIR/${grp_mov}_to_${grp_fix}_${side}/${i}_template${grp_mov}/
  mkdir -p $TRANSITERDIR
  for sub in ${LABEL_IDS[*]}; do
    ln -sf $DATADIR/template${grp_mov}_${side}_${sub}.nii.gz \
      $TRANSITERDIR/template${grp_mov}_${side}_${sub}.nii.gz
  done
  ln -sf $DATADIR/template${grp_mov}_${side}_seg.nii.gz \
    $TRANSITERDIR/template${grp_mov}_${side}_seg.nii.gz

  # reslice seg
  c3d  $DATADIR/template${grp_fix}_${side}_tse.nii.gz \
    $TRANSITERDIR/template${grp_mov}_${side}_seg.nii.gz \
    -interp 0 -reslice-identity \
    -o $TRANSITERDIR/${i}_template${grp_mov}_reslice_to_template${grp_fix}_${side}_seg.nii.gz
  for ((ii=0;ii<${#EXTRAMESHES[*]};ii++)); do
    ln -sf $INITTEMP_DIR/group${grp_mov}/meshwarp/template_${side}_${grp_mov}_${EXTRAMESHES[ii]}.vtk $TRANSITERDIR/${i}_template${grp_mov}_${side}_${EXTRAMESHES[ii]}.vtk
  done

  #c3d $(for sub in ${LABEL_IDS[*]}; do echo $DATADIR/template${grp_mov}_${side}_${sub}.nii.gz; done) \
  #  -vote -type ushort -as SEG \
  #  $DATADIR/template${grp_fix}_${side}_tse.nii.gz \
  #  -push SEG -interp 0 -reslice-identity \
  #  -o $TRANSITERDIR/template_${side}_${grp_mov}.nii.gz

  # go through the path
  warps=""
  invwarps=""
  mkdir -p $TRANSDIR
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

    # perform registration between steps
    PRETRANSITERDIR=$TRANSDIR/${grp_mov}_to_${grp_fix}_${side}/$((i-1))_${id_from}
    TRANSITERDIR=$TRANSDIR/${grp_mov}_to_${grp_fix}_${side}/${i}_${id_to}
    #TRANSITERTMPDIR=$TRANSITERDIR/tmp
    mkdir -p $TRANSITERDIR

    # link the original subject seg
    for sub in ${LABEL_IDS[*]}; do
    ln -sf $DATADIR/${id_to}_${side}_${sub}.nii.gz \
      $TRANSITERDIR/${id_to}_${side}_${sub}.nii.gz
    done
    ln -sf $DATADIR/${id_to}_${side}_seg.nii.gz \
      $TRANSITERDIR/${id_to}_${side}_seg.nii.gz

    # Use ml_affine for nice affine alignment
    #/data/picsl/pauly/wolk/ashs/ext/Linux/bin/ml_affine \
    #  $TRANSITERDIR/${id_to}_${side}_seg.nii.gz \
    #  $PRETRANSITERDIR/template${grp_mov}_${side}_reslice_seg.nii.gz \
    #  $TRANSITERDIR/${id_from}_to_${id_to}_${side}_mlaffine.txt

    /data/picsl/pauly/wolk/ashs/ext/Linux/bin/ml_affine \
      $TRANSITERDIR/${id_to}_${side}_seg.nii.gz \
      $PRETRANSITERDIR/${id_from}_${side}_seg.nii.gz \
      $TRANSITERDIR/${id_from}_to_${id_to}_${side}_mlaffine.txt

    # Convert that to ITK format
    c3d_affine_tool \
      $TRANSITERDIR/${id_from}_to_${id_to}_${side}_mlaffine.txt \
      -oitk \
      $TRANSITERDIR/${id_from}_to_${id_to}_${side}_mlaffine_itk.txt

    CMD=""
    for sub in ${REG_LABELS_1[*]}; do
      CMD="$CMD -m MSQ[$TRANSITERDIR/${id_to}_${side}_${sub}.nii.gz,$PRETRANSITERDIR/${id_from}_${side}_${sub}.nii.gz,$WGT_1]"
    done

    if [[ $ANTs_x_1 == "Y" ]]; then
      c3d $TRANSITERDIR/${id_to}_${side}_seg.nii.gz -dup \
          $PRETRANSITERDIR/${id_from}_${side}_seg.nii.gz \
          -int 0 -reslice-identity \
          -add -binarize -dilate 1 10x10x10vox \
          -o $TRANSITERDIR/${id_from}_to_${id_to}_${side}_mask.nii.gz
      ANTs_mask_1="-x $TRANSITERDIR/${id_from}_to_${id_to}_${side}_mask.nii.gz"
    else
      ANTs_mask_1=""
    fi

    # Perform ANTs registration
    ANTS 3 $CMD \
         -t $ANTs_t_1 \
         -r $ANTs_r_1 \
         -i $ANTs_i_1 \
         $ANTs_G_1 \
         $ANTs_SMOOTH_1 \
         $ANTs_mask_1 \
         -a $TRANSITERDIR/${id_from}_to_${id_to}_${side}_mlaffine_itk.txt \
         --continue-affine 0 \
         $ANTs_all_metrics_1 \
         -o $TRANSITERDIR/${id_from}_to_${id_to}_${side}_pairwise.nii.gz

    # compose warps
    warps="$TRANSITERDIR/${id_from}_to_${id_to}_${side}_pairwiseWarp.nii.gz $TRANSITERDIR/${id_from}_to_${id_to}_${side}_mlaffine_itk.txt $warps"
    invwarps="$invwarps -i $TRANSITERDIR/${id_from}_to_${id_to}_${side}_mlaffine_itk.txt $TRANSITERDIR/${id_from}_to_${id_to}_${side}_pairwiseInverseWarp.nii.gz"

    $ANTSPATH/ComposeMultiTransform 3 \
      $TRANSITERDIR/template${grp_mov}_${side}_InitWarp.nii.gz \
      -R $TRANSITERDIR/${id_to}_${side}_seg.nii.gz \
      $warps

    # apply init warps
    for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        $DATADIR/template${grp_mov}_${side}_${sub}.nii.gz \
        $TRANSITERDIR/template${grp_mov}_${side}_reslice_${sub}.nii.gz \
        -R $TRANSITERDIR/${id_to}_${side}_seg.nii.gz \
        $TRANSITERDIR/template${grp_mov}_${side}_InitWarp.nii.gz

    done

    # Create the segmentation for the template
    c3d $(for sub in ${LABEL_IDS[*]}; do echo $TRANSITERDIR/template${grp_mov}_${side}_reslice_${sub}.nii.gz; done) \
      -vote -type ushort \
      -o $TRANSITERDIR/template${grp_mov}_${side}_reslice_seg.nii.gz

    # resliced to output template
    c3d $DATADIR/template${grp_fix}_${side}_tse.nii.gz \
      $TRANSITERDIR/template${grp_mov}_${side}_reslice_seg.nii.gz \
      -interp 0 -reslice-identity \
      -o $TRANSITERDIR/${i}_template${grp_mov}_reslice_to_template${grp_fix}_${side}_seg.nii.gz

    # generate a surface mesh (TODO)
    $ANTSPATH/ComposeMultiTransform 3 \
      $TRANSITERDIR/template${grp_mov}_${side}_InitInverseWarp.nii \
      -R $DATADIR/template${grp_mov}_${side}_tse.nii.gz \
      $invwarps

    c3d -mcs $TRANSITERDIR/template${grp_mov}_${side}_InitInverseWarp.nii \
      -oo $TMPDIR/comp%02d.nii

    warpmesh -w ants -m ras \
      $TRANSDIR/${grp_mov}_to_${grp_fix}_${side}/1_template${grp_mov}/1_template${grp_mov}_${side}_MRG.vtk \
      $TRANSITERDIR/${i}_template${grp_mov}_${side}_MRG.vtk \
      $TMPDIR/comp00.nii $TMPDIR/comp01.nii $TMPDIR/comp02.nii

  done
    # get warps
    #warps="$PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_pairwiseWarp.nii.gz $PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_pairwiseAffine.txt $warps"

    #invwarps="$invwarps -i $PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_pairwiseAffine.txt $PWDIR/${id_to}_${side}/${id_from}_to_${id_to}/${id_from}_to_${id_to}_${side}_pairwiseInverseWarp.nii.gz"

  # compose transform
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
    for grp_temp in ${groups_temp[*]}; do

      # loop through each subject
      for ((i=0;i<${#IDS[*]};i++)); do

      id=${IDS[i]}
      grp_subj=${GRPS[i]}

      # submit jobs
      qsub -cwd -o $DUMPDIR -j y \
        -N "${PREFIX}_${side}_${grp_temp}_${id}_${grp_subj}" $0 \
        WarpSuperTemplate_sub $side $grp_temp $id $grp_subj
      sleep 0.1

      done
    done
  done

  # wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function AverageFinalTemplate()
{
  # generate the average template
  for side in left right; do
    for grp_temp in ${groups_temp[*]}; do

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
function shape_update_to_finaltemplate()
{
  side=$1
  grp_temp=$2
  FINALWORKDIR=$FINALTEMPDIR/template${grp_temp}/finalwork/

  # Borrowed from ANTS buildtemplateparallel

  # Average the warp fields
  local TEMPWARP=$FINALWORKDIR/template_${side}_warp.nii.gz
  $NEWANTSDIR/AverageImages 3 $TEMPWARP 0 \
    $FINALWORKDIR/*_${side}_totempWarp.nii.gz

  # Scale by -0.25 (gradient step)
  MultiplyImages 3 $TEMPWARP -0.25 $TEMPWARP

  # Create a new template affine
  local TEMPAFF=$FINALWORKDIR/template_${side}_${grp}_Affine.txt
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

function RegisterSuperTemplate()
{
  # choose the best template to be the super template
  PREFIX=RST${expid}
  for grp_temp in ${groups_temp[*]}; do
    FINALWORKDIR=$FINALTEMPDIR/template${grp_temp}/finalwork/
    for ((iter=0;iter<$FINALITER;iter++)); do

      for side in left right; do
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
             "${PREFIX}_${id}_${side}_${grp_temp}" $0 \
             RegisterSuperTemplate_sub $id $side $grp_temp $doants

        done
      done

      # Wait for completion
      qsub -cwd -o $DUMPDIR -j y \
        -hold_jid "${PREFIX}_*" -sync y -b y \
        sleep 1

      # If this is the last iteration, we don't want to recompute the template
      if [[ $iter -lt $((FINALITER-1)) ]]; then
        for side in left right; do

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
          if [[ $doants -eq 1 ]]; then
            shape_update_to_finaltemplate $side $grp_temp
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
    for grp_temp in ${groups_temp[*]}; do
      FINALWORKDIR=$FINALTEMPDIR/template${grp_temp}/finalwork/

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

      id=$(ls $ASHSRUNDIR | grep $id)

      for grp_temp in ${groups_temp[*]}; do

        mkdir -p $FINALTEMPDIR/template${grp_temp}/labelwarp
        # Submit job for this subject
        qsub -cwd -o $DUMPDIR -j y \
          -N "${PREFIX}_${id}_${side}_${grp_temp}" \
          $0 WarpFinalLabel_sub $id $side $grp_temp

      done
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
  grp_temp=$3

  FINALDIR=$FINALTEMPDIR/template${grp_temp}
  FINALWORKDIR=$FINALDIR/finalwork/

  # Apply the warp to each of the labels
  for sub in ${LABEL_IDS[*]}; do

    antsApplyTransforms -d 3 \
      -r $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray_dividedCS.nii.gz \
      -i $FINALWORKDIR/template_${side}_${sub}.nii.gz \
      -o $FINALDIR/labelwarp/${id}_${side}_${sub}_tempfit.nii.gz \
      -t [$ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt,1] \
      -t [$ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt,1] \
      -t $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii \
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
    c3d $ASHSRUNDIR/${id}/final/${id}_${side}_lfseg_corr_usegray_dividedCS.nii.gz \
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
function WarpFinalMesh()
{
  PREFIX=WFM${expid}
  for grp_temp in ${groups_temp[*]}; do
    FINALDIR=$FINALTEMPDIR/template${grp_temp}
    FINALWORKDIR=$FINALDIR/finalwork/

    mkdir -p $FINALDIR/meshwarp
    mkdir -p $FINALDIR/jacobian
    mkdir -p $FINALDIR/meshwarp/template

    # Iterate over side and subject
    for side in left right; do

      ITDIR=$FINALWORKDIR/$(printf iter_%s_%02d $side $((FINALITER-1)) )

      # Generate meshes for the individual subfields
      for sub in ${LABEL_FG[*]}; do

        vtklevelset \
          $ITDIR/template_${side}_${sub}.nii.gz \
          $FINALDIR/meshwarp/template_${side}_${sub}.vtk \
          $subfield_th

      done

      # Generate one mesh for all non-DG non-CS subfields
      for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

        vtklevelset \
          $ITDIR/template_${side}_${EXTRAMESHES[i]}.nii.gz \
          $FINALDIR/meshwarp/template_${side}_${EXTRAMESHES[i]}.vtk \
          $template_th

      done

      cp $FINALDIR/meshwarp/template*.vtk \
         $FINALDIR/meshwarp/template/

      IDS=$(cat $SUBJ_TXT)

      for id in $IDS; do
 
        id=$(ls $ASHSRUNDIR | grep $id)

        # Submit job for this subject
        qsub  -cwd -o $DUMPDIR -j y \
          -N "${PREFIX}_${id}_${side}_${grp_temp}" \
          $0 WarpFinalMesh_sub $id $side $grp_temp

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
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

  # Compose the transformation between the template and the subject
  $ANTSPATH/ComposeMultiTransform 3 \
     $TMPDIR/compose.nii \
     -R $FINALWORKDIR/template_${side}_tse.nii.gz \
     $FINALWORKDIR/${id}_${side}_totempWarp.nii.gz \
     $FINALWORKDIR/${id}_${side}_totempAffine.txt \
     $FINALDIR/work/${id}_${side}_totempWarp.nii.gz \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempWarp.nii \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
     $ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

  # Split the transformations into xyz
  c3d -mcs $TMPDIR/compose.nii -oo $TMPDIR/comp%02d.nii

  # Apply the warp to each of the meshes
  for sub in $ALLSF; do

    # Warp the subfield into subject space
    warpmesh -w ants -m ras \
      $FINALDIR/meshwarp/template_${side}_${sub}.vtk \
      $FINALDIR/meshwarp/${id}_${side}_${sub}_tempfit.vtk \
      $TMPDIR/comp00.nii $TMPDIR/comp01.nii $TMPDIR/comp02.nii

    # Extract the thickness of the subfield
    cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
      -T $FINALDIR/meshwarp/${id}_${side}_${sub}_thickmap.vtk \
      -p $thick_p -e $thick_e \
      $FINALDIR/meshwarp/${id}_${side}_${sub}_tempfit.vtk \
      $FINALDIR/meshwarp/skel_${id}_${side}_${sub}.vtk

  done

  # Compute the Jacobian map
  ANTSJacobian 3 $TMPDIR/compose.nii \
    $FINALDIR/jacobian/${id}_${side}_totemp 1
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
          -N "${PREFIX}_${id}_${side}_${grp_temp}" \
          $0 EvalFinalTemp_sub \
          $id $side $grp_temp $FINALEVALTMPDIR
        sleep 0.1

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # Combine scores
  for grp_temp in ${groups_temp[*]}; do

    FINALEVALDIR=$FINALTEMPDIR/template${grp_temp}/evaluation
    FINALEVALTMPDIR=$FINALEVALDIR/tmp
    for side in left right; do

      echo "ID ${EVALLABELS[*]} ${MESH_EVAL[*]} ${EVALLABELS[*]}" \
        > $FINALEVALDIR/overlap_${side}.txt

      echo "ID ${EVALLABELS[*]} ${MESH_EVAL[*]}" \
        > $FINALEVALDIR/overlap_refseg_${side}.txt

      for id in $IDS; do

        id=$(ls $ASHSRUNDIR | grep $id)
        cat $FINALEVALTMPDIR/${id}_${side}_overlap.txt \
            >> $FINALEVALDIR/overlap_${side}.txt

        if [[ -f $FINALEVALTMPDIR/${id}_${side}_refseg_overlap.txt ]]; then
          cat $FINALEVALTMPDIR/${id}_${side}_refseg_overlap.txt \
            >> $FINALEVALDIR/overlap_refseg_${side}.txt

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
  do_pair $FINALWORKDIR/template_${side}_seg.nii.gz \
    $FINALWORKDIR/${id}_${side}_totemp_reslice_seg.nii.gz
  ALLOVL="$ALLOVL $FULLOVL"

  ###################################
  # save to temporary file
  echo $id $ALLOVL > $EVALTMPDIR/${id}_${side}_overlap.txt

  ###################################
  ALLOVL=""
  # Compute dice overlap between warp seg and manual seg
  FITTED=$FINALDIR/labelwarp/${id}_${side}_seg_tempfit_ref.nii.gz
  REFSEG=$FINALDIR/labelwarp/${id}_${side}_seg_refseg.nii.gz
  if [[ -f $REFSEG ]]; then
  # compute dice
  c3d $FINALDIR/labelwarp/${id}_${side}_seg_tempfit.nii.gz -replace 8 9 -o $FITTED
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
  fi

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
          -N "${PREFIX}_${side}_${grp_temp}_${sub}" \
          $0 DispFinalStat_sub $side $grp_temp $sub

      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
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
    -r $FINALDIR/meshwarp/template_${side}_${sub}.vtk \
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
            -N "${PREFIX}_${exp}_${grp_temp}_${GRPNM[igrp]}" \
            $0 ThickFinalStat_sub $exp $grp_temp $igrp

        done
      done
    #done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
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
        $FINALDIR/meshwarp/template_${side}_${sub}.vtk \
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
        -p 1000 \
        -s P \
        -t 0.01 \
        -d 6 \
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
           -l h_vmem=2.1G,s_vmem=2G \
           -N "${PREFIX}_${grp_temp}_${side}" \
           $0 FinalMeanThickness_sub $grp_temp $side
      sleep 0.1
      #FinalMeanThickness_sub $grp_temp $side

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # generate header
  header="ID"
  for side in left right; do
    for ((i=0;i<${#MESH_LABEL[*]};i++)); do
      sub=${MESH_LABEL[$i]}
      if [[ $sub != "DG" ]]; then
        header="${header},${sub}_${side}"
      fi
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
          $MESHGRPDIR/template_${side}_${EXTRAMESHES[ii]}.vtk \
          $WORKGRPDIR/template_${side}_${sub}.nii.gz \
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

##################################################
function UTMeshSubSample()
{
  PREFIX=GS${expid}
  # Iterate over side and subject
  for grp_temp in ${groups_temp[*]}; do
    for side in left right; do
      
      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=6.1G,s_vmem=6G \
           -N "${PREFIX}_${side}_${grp_temp}" \
           $0 UTMeshSubSample_sub \
           $side $grp_temp
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function UTMeshSubSample_sub()
{
  side=$1
  grp_temp=$2
  FINALDIR=$FINALTEMPDIR/template${grp_temp}
  FINALWORKDIR=$FINALDIR/finalwork/
  FINALMESHDIR=$FINALDIR/meshwarp/
  FINALLABELWARPDIR=$FINALDIR/labelwarp
  ITDIR=$FINALWORKDIR/$(printf iter_%s_%02d $side $((FINALITER-1)) )

  # Generate one mesh for all non-DG non-CS subfields
  for ((i=0;i<${#EXTRAMESHES[*]};i++)); do

    vtklevelset \
      $ITDIR/template_${side}_${EXTRAMESHES[i]}.nii.gz \
      $FINALDIR/meshwarp/template_${side}_${EXTRAMESHES[i]}.stl \
      $template_th

  done

  cp $FINALDIR/meshwarp/template_${side}*.stl \
     $FINALDIR/meshwarp/template/

  # perform subsampling
  SAM=$FINALDIR/meshwarp/template/template_${side}_MRG_sampled.ply
  /data/picsl-build/pauly/vcg/gcc64rel/mesh_poisson_sample \
    $FINALDIR/meshwarp/template/template_${side}_MRG.stl \
    $SAM 2000

  # comvert into vtk mesh
  TEMPLATE=$FINALDIR/meshwarp/template/template_${side}_MRG_sampled.vtk
  NV=$(cat $SAM | grep 'element vertex' | awk '{print $3}')
  
  # Write the header of the VTK file
  echo "# vtk DataFile Version 3.0" > $TEMPLATE
  echo "vtk output" >> $TEMPLATE
  echo "ASCII" >> $TEMPLATE
  echo "DATASET POLYDATA" >> $TEMPLATE
  echo "POINTS ${NV} float" >> $TEMPLATE

  # Write the point coordinates
  grep -A $NV end_header $SAM | tail -n $NV >> $TEMPLATE

}

##################################################
function UTGeodesicShoot()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=GS${expid}
  # Iterate over side and subject
  for grp_temp in ${groups_temp[*]}; do
  for side in left right; do
    for id in $IDS; do
    #for id in DW104; do

      id=$(ls $ASHSRUNDIR | grep $id)
      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=6.1G,s_vmem=6G \
           -N "${PREFIX}_${id}_${side}_${grp_temp}" \
           $0 UTGeodesicShoot_sub \
           $id $side $grp_temp
      sleep 0.1

    done
  done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function UTGeodesicShoot_sub()
{
  id=$1
  side=$2
  grp_temp=$3
  FINALDIR=$FINALTEMPDIR/template${grp_temp}
  FINALWORKDIR=$FINALDIR/finalwork/
  FINALMESHDIR=$FINALDIR/meshwarp/
  FINALLABELWARPDIR=$FINALDIR/labelwarp
  UTGEOSHOOTDIR=$FINALDIR/UTGeoShoot
  SUBJDIR=$UTGEOSHOOTDIR/$id
  mkdir -p $SUBJDIR

  #if [ -f $SUBJDIR/${id}_${side}_momentum.vtk ]; then
  #  exit
  #fi

  # warp decimated template mesh back to subject space
  # Compose the transformation between the template and the subject
  $ANTSPATH/ComposeMultiTransform 3 \
     $TMPDIR/compose.nii \
     -R $FINALWORKDIR/template_${side}_tse.nii.gz \
     $FINALWORKDIR/${id}_${side}_totempWarp.nii.gz \
     $FINALWORKDIR/${id}_${side}_totempAffine.txt \
     $FINALDIR/work/${id}_${side}_totempWarp.nii.gz \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempWarp.nii \
     $ASHSRUNDIR/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
     $ASHSRUNDIR/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

  # Split the transformations into xyz
  #c3d -mcs $TMPDIR/compose.nii -oo $TMPDIR/comp%02d.nii

  # Apply the warp to each of the meshes
  for sub in MRG; do

    # Warp the subfield into subject space
    greedy -d 3 \
      -rf $ASHSRUNDIR/$id/mprage.nii.gz \
      -rs $FINALDIR/meshwarp/template/template_${side}_MRG_sampled.vtk \
          $SUBJDIR/${id}_${side}_${sub}_tempfit.vtk \
      -r $TMPDIR/compose.nii

    #warpmesh -w ants -m ras \
    #  $FINALDIR/meshwarp/template/template_${side}_MRG_sampled.vtk \
    #  $SUBJDIR/${id}_${side}_${sub}_tempfit.vtk \
    #  $TMPDIR/comp00.nii $TMPDIR/comp01.nii $TMPDIR/comp02.nii

    # Extract the thickness of the subfield
    #cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
    #  -T $SUBJDIR/${id}_${side}_${sub}_thickmap.vtk \
    #  -p $thick_p -e $thick_e \
    #  $SUBJDIR/${id}_${side}_${sub}_tempfit.vtk \
    #  $SUBJDIR/skel_${id}_${side}_${sub}.vtk

  done

  # rigid regsitration between template ans subject seg
  #greedy -d 3 \
  #  -i $FINALLABELWARPDIR/${id}_${side}_seg_tempfit.nii.gz \
  #     $FINALWORKDIR/template_${side}_seg.nii.gz \
  #  -o $SUBJDIR/template_to_${id}_${side}_rigid.mat \
  #  -a -ia-identity -dof 6 -n 100x50x20 -m NCC 4x4x4
  #warpmesh $SUBJDIR/${id}_${side}_MRG_tempfit.vtk \
  #  $SUBJDIR/${id}_${side}_MRG_tempfit_rigid.vtk \
  #  $SUBJDIR/template_to_${id}_${side}_rigid.mat 

  # use procrustes to align the two mesh
  vtkprocrustes $SUBJDIR/${id}_${side}_MRG_tempfit.vtk \
    $FINALDIR/meshwarp/template/template_${side}_MRG_sampled.vtk \
    $SUBJDIR/template_to_${id}_${side}_rigid.mat
  greedy -d 3 \
    -rf $ASHSRUNDIR/$id/mprage.nii.gz \
    -rs $SUBJDIR/${id}_${side}_MRG_tempfit.vtk \
        $SUBJDIR/${id}_${side}_MRG_tempfit_rigid.vtk \
    -r $SUBJDIR/template_to_${id}_${side}_rigid.mat
        
  #warpmesh $SUBJDIR/${id}_${side}_MRG_tempfit.vtk \
  #  $SUBJDIR/${id}_${side}_MRG_tempfit_rigid.vtk \
  #  $SUBJDIR/template_to_${id}_${side}_rigid.mat

  # perform geodesic shooting
  mkdir -p $SUBJDIR/movie_${side}
  lmshoot -d 3 \
    -m $FINALDIR/meshwarp/template/template_${side}_MRG_sampled.vtk \
       $SUBJDIR/${id}_${side}_MRG_tempfit_rigid.vtk \
    -o $SUBJDIR/${id}_${side}_momentum.vtk \
    -s 2.0 -l 5000 \
    -n 40 -i 200 0 \
    -f -r 1 \
    -O $SUBJDIR/movie_${side}/${id}_${side}_movie%04d.vtk
 
  # convert the shooting result into a warp
  #lmtowarp -d 3 -n 40 \
  #  -r $FINALWORKDIR/template_${side}_seg.nii.gz \
  #  -m $SUBJDIR/${id}_${side}_momentum.vtk \
  #  -o $SUBJDIR/${id}_${side}_momentum_warp.nii.gz \
  #  -s 2.0
}

##################################################
function UTGShootAvg()
{
  IDS=$(cat $SUBJ_TXT)

  PREFIX=GSA${expid}
  # Iterate over side and subject
  for grp_temp in ${groups_temp[*]}; do
    for side in left right; do

      # Submit job for this subject
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=20.1G,s_vmem=20G \
           -N "${PREFIX}_${grp_temp}_${side}" $0 \
           UTGShootAvg_sub $side $grp_temp
      sleep 0.1

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function UTGShootAvg_sub()
{
  side=$1
  grp_temp=$2
  FINALDIR=$FINALTEMPDIR/template${grp_temp}
  FINALWORKDIR=$FINALDIR/finalwork/
  FINALMESHDIR=$FINALDIR/meshwarp/
  FINALLABELWARPDIR=$FINALDIR/labelwarp
  UTGEOSHOOTDIR=$FINALDIR/UTGeoShoot
  mkdir -p $UTGEOSHOOTDIR/shavg

  # Average the momentum maps
  avgmesharr \
    $UTGEOSHOOTDIR/*/*_${side}_momentum.vtk \
    InitialMomentum \
    $FINALMESHDIR/template/template_${side}_MRG_sampled.vtk \
    $UTGEOSHOOTDIR/shavg/template_${side}_MRG_average_momenta.vtk
  
  # Perform the shooting and generate warp
  lmtowarp -d 3 -n 40 \
    -r $FINALWORKDIR/template_${side}_seg.nii.gz \
    -m $UTGEOSHOOTDIR/shavg/template_${side}_MRG_average_momenta.vtk \
    -o $UTGEOSHOOTDIR/shavg/template_${side}_MRG_average_momenta_warp.nii.gz \
    -s 2.0

  # Apply the warp to the template
  greedy -d 3 \
    -rf $FINALWORKDIR/template_${side}_seg.nii.gz \
    -rs $FINALMESHDIR/template_${side}_MRG.vtk \
        $UTGEOSHOOTDIR/shavg/template_${side}_MRG_Gshoot.vtk \
    -r $UTGEOSHOOTDIR/shavg/template_${side}_MRG_average_momenta_warp.nii.gz


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
