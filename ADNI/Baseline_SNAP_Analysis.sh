#//bin/bash
#$ -S /bin/bash
set -e
#set -e -x
#set -x -e

##############################################
# Setup environment
#source ~longxie/.bash_profile

# Software PATH
ANTSPATH=/data/picsl/longxie/pkg/bin/antsbin/bin
C3DPATH=/data/picsl/longxie/pkg/c3d_tool/bin
FSLROOT=/data/picsl/longxie/pkg/fsl
FSLPATH=$FSLROOT/bin
JAVAEXEBIN=/data/picsl/longxie/CardiacMR/JavaExecutable
JAVABIN=/data/picsl/longxie/pkg/jdk1.7.0_71/bin
#export ASHS_ROOT=/data/picsl/pauly/wolk/ashs
#export ASHS_ROOT=/data/picsl/longxie/pkg/ashs-fast
export ASHS_ROOT=/data/picsl/pauly/wolk/ashs-fast
export PATH=$PATH:$ASHS_ROOT/bin
MATLAB_BIN=/share/apps/matlab/R2016a/bin/matlab
export ALOHA_ROOT=/data/picsl/longxie/pkg/aloha
SRMATLABCODEDIR=/home/longxie/ASHS_PHC/SuperResolution/SRToolBox/matlabfunction
SHOOTWPRIORDIR=/home/longxie/LMTemplateMatching/PointSetGeodesicShooting/build_release
LMTOWARPDIR=/home/longxie/LMTemplateMatching/PointSetUtilities/build_release
export FREESURFER_HOME=$CFNAPPS/freesurfer/6.0.0
source $FREESURFER_HOME/SetUpFreeSurfer.sh
SRPATH=/home/longxie/pkg/PatchSuperResolution/build_release

TEMPLATEDIR=/data/jag/wolk_group/ADNI_longitudinal-Templates/Normal/

##############################################
# Directories
ROOT=/data/picsl/longxie/ADNI2018/
CODEDIR=$ROOT/scripts
ANALYSISDIR=$ROOT/analysis_input
#ALLINFO_TXT=$ANALYSISDIR/MRI3TListWithNIFTIPath_11262018.tsv
ALLINFO_TXT=$ANALYSISDIR/ADNIBL_Merge20190605_21_SNAP.csv
#ALLINFO_TXT=$ANALYSISDIR/MRILIST_T1T2_012_S_4128.tsv
ALLDATADIR=$ROOT/dataset
#LOCALDATADIR=$ROOT/dataset_local
DATADIR=$ROOT/dataset
OUTPUTDIR=$ROOT/output_BaselineSNAP_09112019
STATSDIR=$OUTPUTDIR/stats
DESIGNDIR=$STATSDIR/design
SHAPESTATSDIR=$OUTPUTDIR/shapestats
LOGDIR=$OUTPUTDIR/log
DUMPDIR=$OUTPUTDIR/dump
mkdir -p $LOGDIR

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

# multi template directories
GSTEMPDIR=/home/longxie/ASHS_T1/thickness/exp/exp601/
ATLASLIST=$GSTEMPDIR/GSTemplate/MST/paths/IDSide.txt
GTGROUPS=/home/longxie/ASHS_T1/thickness/group/group_xval_2group_all.txt
FLIPLRMAT=$CODEDIR/flip_LR.mat
FLIPLRITK=$CODEDIR/flip_LR_itk.txt

MSTUTTEMPDIR=/home/longxie/ASHS_T1/thickness/exp/exp702

#############################################
function main()
{
  reset_dir

  # Perform GLM analysis
  RegionalThickStatsUT






}

#############################################
function ReFormateDate()
{
  indate=$1

  if [[ $indate == "" || $indate == " " ]]; then

    outdate=$indate

  else

    DD=$(date -d "$indate" '+%d')
    MM=$(date -d "$indate" '+%m')
    YYYY=$(date -d "$indate" '+%Y')
    outdate="${YYYY}-${MM}-${DD}"

  fi

  echo $outdate
}

######################################################
function RegionalThickStatsUT()
{
  rm -rf $DESIGNDIR
  mkdir -p $DESIGNDIR

  # columns
  INFOCSV=$ALLINFO_TXT
  N=$(cat $ALLINFO_TXT | wc -l)
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  DX1Col=$(csvcol.sh $ALLINFO_TXT DX1)
  FourgroupsCol=$(csvcol.sh $ALLINFO_TXT Fourgroups)
  PTSEXCol=$(csvcol.sh $ALLINFO_TXT PTSEX)
  AgeatMRICol=$(csvcol.sh $ALLINFO_TXT AgeatMRI)
  MTL_LEFTCol=$(csvcol.sh $ALLINFO_TXT MTL_LEFT)
  MTL_RIGHTCol=$(csvcol.sh $ALLINFO_TXT MTL_RIGHT)

  # get information
  IDS=($(awk -F "," '{print $1}' $INFOCSV))
  ICVS=($(awk -F "," '{print $2}' $INFOCSV))
  DX=($(awk -F "," '{print $122}' $INFOCSV))
  AGE=($(awk -F "," '{print $107}' $INFOCSV))
  QCMTLBOTH=($(awk -F "," '{print $118}' $INFOCSV))
  QCMRI=($(awk -F "," '{print $112}' $INFOCSV))

  for ((i=2;i<${N};i++)); do

    # get information
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    dx1=$(echo $ROW | cut -f $DX1Col -d ",")
    fourgrp=$(echo $ROW | cut -f $FourgroupsCol -d ",")
    sex=$(echo $ROW | cut -f $PTSEXCol -d ",")
    age=$(echo $ROW | cut -f $AgeatMRICol -d ",")
    leftqc=$(echo $ROW | cut -f $MTL_LEFTCol -d ",")
    rightqc=$(echo $ROW | cut -f $MTL_RIGHTCol -d ",")
    echo "Processing $i $rid"

    for side in left right; do 

      # choose qc
      if [[ $side == "left" ]]; then
        qc=$leftqc
      else
        qc=$rightqc
      fi

      if [[ $qc != 0 || $dx1 == 3 || $fourgrp == 999 || $fourgrp == 1 ]]; then
        continue
      fi
   
      # preparation
      postfix="group_age_sexcr"
      if [[ $dx1 == 1 && $fourgrp == 0 ]]; then
        str="${id}|${PREFIX} 0 1 $age $sex"
        echo $str >> $DESIGNDIR/UTdesign_${side}_Control_3vs0_${postfix}.txt
        echo $str >> $DESIGNDIR/UTdesign_${side}_Control_2vs0_${postfix}.txt
      elif [[ $dx1 == 1 && $fourgrp == 2 ]]; then
        str="${id}|${PREFIX} 1 0 $age $sex"
        echo $str >> $DESIGNDIR/UTdesign_${side}_Control_2vs0_${postfix}.txt
        str="${id}|${PREFIX} 0 1 $age $sex"
        echo $str >> $DESIGNDIR/UTdesign_${side}_Control_3vs2_${postfix}.txt
      elif [[ $dx1 == 1 && $fourgrp == 3 ]]; then
        str="${id}|${PREFIX} 1 0 $age $sex"
        echo $str >> $DESIGNDIR/UTdesign_${side}_Control_3vs0_${postfix}.txt
        echo $str >> $DESIGNDIR/UTdesign_${side}_Control_3vs2_${postfix}.txt
      elif [[ $dx1 == 2 && $fourgrp == 0 ]]; then
        str="${id}|${PREFIX} 0 1 $age $sex"
        echo $str >> $DESIGNDIR/UTdesign_${side}_MCI_3vs0_${postfix}.txt
        echo $str >> $DESIGNDIR/UTdesign_${side}_MCI_2vs0_${postfix}.txt
      elif [[ $dx1 == 2 && $fourgrp == 2 ]]; then
        str="${id}|${PREFIX} 1 0 $age $sex"
        echo $str >> $DESIGNDIR/UTdesign_${side}_MCI_2vs0_${postfix}.txt
        str="${id}|${PREFIX} 0 1 $age $sex"
        echo $str >> $DESIGNDIR/UTdesign_${side}_MCI_3vs2_${postfix}.txt
      elif [[ $dx1 == 2 && $fourgrp == 3 ]]; then
        str="${id}|${PREFIX} 1 0 $age $sex"
        echo $str >> $DESIGNDIR/UTdesign_${side}_MCI_3vs0_${postfix}.txt
        echo $str >> $DESIGNDIR/UTdesign_${side}_MCI_3vs2_${postfix}.txt
      else
        continue
      fi

    done
  done

  # echo contrast
  for type in Control MCI; do
  echo "1 -1 0 0" > \
    $DESIGNDIR/UTcontrast_${type}_3vs0_${postfix}_3-0.txt
  echo "-1 1 0 0" > \
    $DESIGNDIR/UTcontrast_${type}_3vs0_${postfix}_0-3.txt
  echo "1 -1 0 0" > \
    $DESIGNDIR/UTcontrast_${type}_2vs0_${postfix}_2-0.txt
  echo "-1 1 0 0" > \
    $DESIGNDIR/UTcontrast_${type}_2vs0_${postfix}_0-2.txt
  echo "1 -1 0 0" > \
    $DESIGNDIR/UTcontrast_${type}_3vs2_${postfix}_3-2.txt
  echo "-1 1 0 0" > \
    $DESIGNDIR/UTcontrast_${type}_3vs2_${postfix}_2-3.txt
  done

  # submit jobs to run statistical analysis
  PREFIX=RTS${expid}
  if [[ 1  == 1 ]]; then
  for diag in Control MCI; do
    for side in left right; do
      rm -rf $STATSDIR/MTtemplate_${diag}_${side}
      for design in $(ls $DESIGNDIR | grep UTdesign | grep $side | grep $diag); do

        exp=$(echo $design | sed -e "s/^UTdesign_${side}_${diag}_//" | sed -e "s/\.txt//")

        qsub -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -l h_vmem=4.1G,s_vmem=4G \
          -N "${PREFIX}_${exp}_${side}_${diag}_MRG_UT" \
          $0 RegionalThickStats_sub $exp $side $diag UT
        echo "$exp $side $diag UT"
        sleep 0.1
      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
  fi
}

function RegionalThickStats_sub()
{
  exp=$1
  side=$2
  diag=$3
  type=$4

  # Create the work directory for this analysis
  WORK=$STATSDIR/${type}template_${diag}_${side}/design_${exp}_MRG
  mkdir -p $WORK

  # Get the list of subjects
  DESIGNTXT=$DESIGNDIR/${type}design_${side}_${diag}_${exp}.txt
  SUBJ=$(cat $DESIGNTXT | awk '{print $1}')

  # Generate the design matrix for meshglm
  cat $DESIGNTXT | awk '{$1=""; print}' > $WORK/design_${side}.txt

  # combine thickness measurements
  if [[ $type == "MT" ]]; then
  MESHES=$(for idtime in $SUBJ; do \
    id=$(echo $idtime | cut -d '|' -f 1)
    tp=$(echo $idtime | cut -d '|' -f 2)
    echo " $(ls $DATADIR/$id/$tp/ASHST1_MTLCORTEX_MSTTHK/GeoShoot/${side}/template_${diag}_to_${id}_${side}_GSShoot_MRG_thickmap.vtk) "; done)
  else
  MESHES=$(for idtime in $SUBJ; do \
    id=$(echo $idtime | cut -d '|' -f 1)
    tp=$(echo $idtime | cut -d '|' -f 2)
    mesh=$(ls $DATADIR/$id/$tp/ASHST1_MTLCORTEX_MSTTHK/MTUTTHK/${id}_${tp}_${side}_UT_template_fitted_mesh.vtk)
    if [[ $mesh == "" ]]; then
      echo "$id $tp UT thick map is missing"
      exit
    fi
    echo " $mesh "; done)
  fi

  #mesh=$(ls $DATADIR/$id/$tp/ASHST1_MTLCORTEX_MSTTHK/GeoShootUT/${side}/template_to_${id}_${side}_GSShoot_MRG_thickmap.vtk)

  if [[ $type == "MT" ]]; then
    TEMPLATE=$GSTEMPDIR/GSTemplate/gshoot/template_${diag}/template/iter_2/template_${diag}_gshoot_MRG.vtk
  else
    #TEMPLATE=$MSTUTTEMPDIR/GSTemplate/gshoot/template_1/template/iter_2/template_1_gshoot_MRG.vtk 
    TEMPLATE=/home/longxie/ASHS_T1/pipeline_package/MultiTempThkMTLWithHippo/template/GSUTemplate/gshoot/template_1/template/iter_2/template_1_gshoot_MRG.vtk
  fi
  mesh_merge_arrays -r \
    $TEMPLATE \
    $WORK/thick_${side}_${diag}_MRG.vtk Thickness $MESHES

    # Go through the list of contrasts
  for con in $(ls $DESIGNDIR | grep "${type}contrast_${diag}_${exp}_.*\.txt"); do

    # Get the suffix for this contrast
    suffix=$(echo $con | sed -e "s/^${type}contrast_${diag}_${exp}_//" -e "s/\.txt//")

    # Create the directory for this contrast
    CWORK=$WORK/contrast_${suffix}
    mkdir -p $CWORK

    FULLNM="design_${exp}_con_${suffix}"

    # Copy the contrast
    cp $DESIGNDIR/$con $CWORK/contrast.txt

    # Build the list of meshes to include
    MESHPARAM="-m $WORK/thick_${side}_${diag}_MRG.vtk $CWORK/thickstat_${FULLNM}_${side}_${diag}_MRG.vtk"
    meshglm $MESHPARAM \
      -g $WORK/design_${side}.txt $CWORK/contrast.txt \
      -a Thickness \
      -d 4 \
      -p 1000 \
      -s P \
      -t 0.05 \
      -e
#      -p 1000 \
#      -s T \
#      -t 2.4 \
#      -e

  done
}

######################################################
function reset_dir()
{
  rm -rf $DUMPDIR
  mkdir -p $DUMPDIR
}

######################################################
if [[ $# -lt 2 ]]; then

  main

else

  set -e -x
  cmd=$1
  shift
  $cmd $@

fi


