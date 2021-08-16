#//bin/bash
#$ -S /bin/bash
#set +e
#set -e -x
set -e

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
ROOT=/home/longxie/Janssen
CODEDIR=$ROOT/scripts
ANALYSISDIR=$ROOT/analysis_input
ALLINFO_TXT=$ANALYSISDIR/subj_list_11032019.txt
DATAROOT=/data/jux/wolk_group/janssen/CPSS_MRI_vols_UPenn/UPenn/
OUTPUTDIR=$ROOT/output_runT2_09282020
LOGDIR=$OUTPUTDIR/log
INFODIR=$OUTPUTDIR/info
QCDIR=$OUTPUTDIR/qc
DUMPDIR=$OUTPUTDIR/dump
mkdir -p $LOGDIR
mkdir -p $INFODIR

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

ASHST13TLABELIDS=(AHippo PHippo ERC BA35 BA36 PHC)
ASHST13TLABELNUMS=(1     2      10  11   12   13)
ASHST23TLABELIDS=(CA1 CA2 DG CA3 Tail SUB ERC BA35 BA36 PHC)
ASHST23TLABELNUMS=(1  2   3  4   5    8   9   10   11   12)

#############################################
function main()
{
  reset_dir

  ################################
  # summarize data
  ################################
  #CheckData

  ###############################
  # Run trim neck
  ###############################
  #TrimNeck

  ###############################
  # Run ASHS T2 (new atlas)
  #ASHST2SegSUB
  #CleanCortexASHS3TT2Seg

  ###############################
  # Clean up ASHS T2
  ###############################
  #CleanASHST2

  #summarize

  ###############################
  # Copy PNG file
  ###############################
  # Copy PNG to file
  QCASHST2Seg

}

######################################################
function CheckData()
{
  # Load id number
  Nbegin=2
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=20

  # initialize tmp dir
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # Submit job to copy data
  PREFIX=LD${expid}
  for ((i=$Nbegin;i<=${N};i++)); do

    qsub -cwd -o $DUMPDIR -j y \
         -q all.q \
         -l h_vmem=2.1G,s_vmem=2G \
         -N "${PREFIX}_${i}" \
         $0 CheckData_sub $i $OUTTMPDIR
    sleep 0.03

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # get all info
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  echo "$header,T1Orientation,T1Dimension,AnalysisType,T1,T2"  > $OUTPUTDIR/checkdata.csv
  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/checkdata.csv
}

function CheckData_sub()
{
  i=$1
  OUTTMPDIR=$2
  ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')

  # columns
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)

  # get information
  id=$(echo $ROW | cut -f $IDCol -d ",")
  SUBJDIR=$DATAROOT/$id
  T1NIFTI=$(ls $SUBJDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz)
  T2NIFTI=$(ls $SUBJDIR/${id}.Screening_MR.High_Res_T2.High-Resolution_2D_Coronal_T2.nii.gz)
 
  # if file does not exist, skip
  if [[ $T1NIFTI == ""  ]]; then
    ROW="$ROW,,,"
  else

    ORI=$(c3d $T1NIFTI -info \
              | cut -d ';' -f 5 | cut -d ' ' -f 5)
    if [[ $ORI == "Oblique," ]]; then
      ORI=$(c3d $T1NIFTI -info \
              | cut -d ';' -f 5 | cut -d ' ' -f 8)
    fi
    DIM=$(c3d $T1NIFTI -info \
      | cut -d ':' -f 2 | cut -d ';' -f 1 \
      | sed -e 's/, /_/g' | sed -e 's/ dim = //g')
  fi
  
  ROW="$ROW,$ORI,$DIM,T2Baseline,$T1NIFTI,$T2NIFTI"
  echo $ROW > $OUTTMPDIR/$(printf %04d $i)_${id}.csv
}

######################################################
function TrimNeck()
{
  # Load id number
  Nbegin=2
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=2

  # columns
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)


  # Submit job to copy data
  JPREFIX=TN
  for ((i=$Nbegin;i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # get information
    id=$(echo $ROW | cut -f $IDCol -d ",")
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    PREFIX="Screening"
    SUBJDIR=$DATAROOT/$id

    # perform trim neck
    TPDIR=$SUBJDIR
    echo "TrimNeck Processing:$i,$rid,$PREFIX"
    T1=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
    T1TRIM=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim.nii.gz

    if [[ -f $T1 && ! -f $T1TRIM ]]; then
      echo "        Adding $rid $PREFIX "
      echo "TrimNeck: Adding $rid $PREFIX " >> $LOGDIR/trimneck_log.txt
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${id}_${PREFIX}" \
           $0 TrimNeck_sub \
             $TPDIR \
             $T1 \
             $T1TRIM
      sleep 0.02
    fi

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}

function TrimNeck_sub()
{
  TPDIR=$1
  T1=$2
  T1TRIM=$3

  if [[ ! -f $T1TRIM ]]; then

    /data/picsl/pauly/bin/trim_neck_rf.sh \
      $T1 $T1TRIM

    echo "Using /data/picsl/pauly/bin/trim_neck_rf.sh" > \
      $TPDIR/neck_trimed_version.txt

  fi
}

######################################################
function ASHST2SegSUB()
{
  # Load id number
  Nbegin=2
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11
  ATLAS=/data/picsl/pauly/wolk/abc_prisma/exp03/final
  Nrun=10

  echo "" >> $LOGDIR/ASHST2_baseline_log.txt
  echo "$(date)" >> $LOGDIR/ASHST2_baseline_log.txt

  # Submit job to copy data
  JPREFIX=ASHST2
  for ((i=$Nbegin;i<=${N};i++)); do

    ROW=$(cat -A $OUTPUTDIR/checkdata.csv | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # get information
    id=$(echo $ROW | cut -f 1 -d,)
    Type=$(echo $ROW | cut -f 4 -d,)
    PREFIX="Screening"
    SUBJDIR=$DATAROOT/$id

    # perform seg
    TPDIR=$SUBJDIR
    echo "ASHS T2 Processing:$i,$rid,$PREFIX"
    if [[ $Type != "T2Baseline" && $Type != "T2Longitudinal" ]]; then
      continue
    fi
    T1=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
    T1TRIM=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim.nii.gz
    T1SR=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim_denoised_SR.nii.gz
    T2=$TPDIR/${id}.Screening_MR.High_Res_T2.High-Resolution_2D_Coronal_T2.nii.gz
    ASHST2DIR=$TPDIR/ASHST2
    
    if [[ -f $T1TRIM && -f $T2 && \
       ! -f $ASHST2DIR/final/${id}_left_heur_volumes.txt ]]; then

      echo "   Adding $id "
      echo "ASHS T2: Adding $id  " >> $LOGDIR/ASHST2_baseline_log.txt

      mkdir -p $ASHST2DIR
      $ASHS_ROOT/bin/ashs_main.sh \
        -a $ATLAS -d -T -I $id -g $T1TRIM -f $T2 \
        -Q -s 1-7 \
        -z $CODEDIR/ashs-fast-z.sh \
        -w $ASHST2DIR &
      sleep 0.1

      ALLJOBS=($(jobs -p))
      while [ ${#ALLJOBS[*]} -ge $Nrun ]; do
        sleep 60
        ALLJOBS=($(jobs -p))
      done
    fi

  done
}

######################################################
function CleanASHST2()
{
  # Load id number
  Nbegin=2
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11

  # Submit job to copy data
  JPREFIX=CLASHST2
  for ((i=$Nbegin;i<=${N};i++)); do

    ROW=$(cat -A $OUTPUTDIR/checkdata.csv | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # get information
    id=$(echo $ROW | cut -f $IDCol -d ",")
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    PREFIX="Screening"
    SUBJDIR=$DATAROOT/$id

    # perform seg
    TPDIR=$SUBJDIR
    echo "ASHS T2 Processing:$i,$rid,$PREFIX"
    if [[ $Type != "T2Baseline" && $Type != "T2Longitudinal" ]]; then
      continue
    fi
    T1=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
    T1TRIM=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim.nii.gz
    T1SR=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim_denoised_SR.nii.gz
    T2=$TPDIR/${id}.Screening_MR.High_Res_T2.High-Resolution_2D_Coronal_T2.nii.gz
    ASHST2DIR=$TPDIR/ASHST2

    if [[ -f $ASHST2DIR/final/${id}_left_lfseg_heur.nii.gz && -f $ASHST2DIR/tse.nii.gz ]]; then

      echo "   Cleaning $rid $PREFIX "
      echo "Cleaning ASHS T2: $rid $PREFIX " >> $LOGDIR/ASHST2_baseline_clean_log.txt

      rm -rf $ASHST1DIR/affine_t1_to_template \
        $ASHST1DIR/ants_t1_to_temp \
        $ASHST1DIR/bootstrap \
        $ASHST1DIR/dump \
        $ASHST1DIR/multiatlas \
        $ASHST1DIR/mprage* \
        $ASHST1DIR/tse.nii.gz \
        $ASHST1DIR/tse_to_chunktemp*.nii.gz

    fi

  done
}

######################################################
function CleanCortexASHS3TT2Seg()
{
  # Load id number
  Nbegin=2
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11

  IDCol=$(csvcol.sh $OUTPUTDIR/checkdata.csv ID)
  TYPECol=$(csvcol.sh $OUTPUTDIR/checkdata.csv AnalysisType)

  # Submit job to copy data
  JPREFIX=SUM
  for ((i=$Nbegin;i<=${N};i++)); do

    ROW=$(cat -A $OUTPUTDIR/checkdata.csv | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # get information
    id=$(echo $ROW | cut -f $IDCol -d ",")
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    PREFIX="Screening"
    SUBJDIR=$DATAROOT/$id

    # perform seg
    TPDIR=$SUBJDIR
    echo "ASHS T2 Clean Cortex Processing:$i,$id,$PREFIX"
    if [[ $Type != "T2Baseline" && $Type != "T2Longitudinal" ]]; then
      continue
    fi
    T1=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
    T1TRIM=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim.nii.gz
    T1SR=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim_denoised_SR.nii.gz
    T2=$TPDIR/${id}.Screening_MR.High_Res_T2.High-Resolution_2D_Coronal_T2.nii.gz
    ASHST2DIR=$TPDIR/ASHST2

    for side in left right; do
    if [[ -f $T2 && \
          -f $ASHST2DIR/final/${id}_${side}_lfseg_corr_nogray.nii.gz ]]; then
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${i}_${id}_${PREFIX}" \
           $0 CleanCortexASHS3TT2Seg_sub \
           $ASHST2DIR/final/${id}_${side}_lfseg_corr_nogray.nii.gz \
           $ASHST2DIR/final/${id}_${side}_lfseg_corr_nogray_CortexCleaned.nii.gz
      sleep 0.02
    fi
    done

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
}

function CleanCortexASHS3TT2Seg_sub()
{
  in=$1
  out=$2
  INFO=$(c3d $in -replace 16 9 -thresh 9 12 1 0 -o $TMPDIR/temp.nii.gz -info)
  NS=$(echo $INFO | sed -e "s/.*dim = .//g" -e "s/.;.*bb.*//g" | awk -F ',' '{print $3}')
  # -slice z 0 does not work -- memory allocation error, so skipping slice 0
  # slicecmd=$(for((i=0;i<$NS;i++)); do echo "-push X -slice z $i -voxel-sum "; done)
  slicecmd=$(for((i=0;i<$NS;i++)); do echo "-push X -slice z $i -voxel-sum "; done)
  c3d $TMPDIR/temp.nii.gz -popas X $slicecmd | grep Voxel | awk '{print $3}' > $TMPDIR/counts.txt
  NNZ=$(cat $TMPDIR/counts.txt | grep -v '^0$' | wc -l)
  MEDIAN=$(cat $TMPDIR/counts.txt | grep -v '^0$' | sort -n | tail -n $((NNZ/2)) | head -n 1)
  CUTOFF=$((MEDIAN / 6))
  RULE=$(cat $TMPDIR/counts.txt | awk "{print NR-1,int(\$1 < $CUTOFF)}")
  c3d $in $TMPDIR/temp.nii.gz -copy-transform -cmv -replace $RULE -popas M $in \
    -replace 16 9 -thresh 9 12 1 0 -push M -times -voxel-sum -scale -1 -shift 1 \
    $in -times -o $out
  NLEFT=$(cat $TMPDIR/counts.txt | awk "\$1 > $CUTOFF {k++} END {print k}")
  echo $NLEFT
#  for((;;)); do
#  sleep 999999999
#  done

}

######################################################
function summarize()
{
  # Load id number
  Nbegin=2
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11

  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  IDCol=$(csvcol.sh $OUTPUTDIR/checkdata.csv ID)
  TYPECol=$(csvcol.sh $OUTPUTDIR/checkdata.csv AnalysisType)

  # Submit job to copy data
  JPREFIX=SUM
  for ((i=$Nbegin;i<=${N};i++)); do

    ROW=$(cat -A $OUTPUTDIR/checkdata.csv | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # get information
    id=$(echo $ROW | cut -f $IDCol -d ",")
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    PREFIX="Screening"
    SUBJDIR=$DATAROOT/$id

    # perform seg
    TPDIR=$SUBJDIR
    echo "ASHS T2 Processing:$i,$rid,$PREFIX"
    if [[ $Type != "T2Baseline" && $Type != "T2Longitudinal" ]]; then
      continue
    fi
    T1=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
    T1TRIM=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim.nii.gz
    T1SR=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim_denoised_SR.nii.gz
    T2=$TPDIR/${id}.Screening_MR.High_Res_T2.High-Resolution_2D_Coronal_T2.nii.gz
    ASHST2DIR=$TPDIR/ASHST2

    if [[ -d $ASHST2DIR ]]; then
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${i}_${id}_${PREFIX}" \
           $0 summarize_sub $i $OUTTMPDIR
      sleep 0.02
    fi

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1

  # get all info
  header="ID,Type,ICV"
  header="$header,L_Modi_Rater_T2_3T,R_Modi_Rater_T2_3T"
  for type in NSlice VOL; do
  for side in L R; do
  for sub in ${ASHST23TLABELIDS[*]}; do
    header="$header,${side}_${sub}_${type}_ASHST2_3T"
  done
  done
  done
  side="M"
  type="VOL"
  for sub in ${ASHST23TLABELIDS[*]}; do
    header="$header,${side}_${sub}_${type}_ASHST2_3T"
  done

  echo $header > $OUTPUTDIR/Janssen_ASHST2_measurements.csv
  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/Janssen_ASHST2_measurements.csv
}

function summarize_sub()
{
  i=$1
  OUTTMPDIR=$2 

  ROW=$(cat -A $OUTPUTDIR/checkdata.csv | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

  IDCol=$(csvcol.sh $OUTPUTDIR/checkdata.csv ID)
  TYPECol=$(csvcol.sh $OUTPUTDIR/checkdata.csv AnalysisType)

  # get information
  id=$(echo $ROW | cut -f $IDCol -d ",")
  Type=$(echo $ROW | cut -f $TYPECol -d ",")
  PREFIX="Screening"
  SUBJDIR=$DATAROOT/$id
  T1=$SUBJDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
  T1TRIM=$SUBJDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim.nii.gz
  T1SR=$SUBJDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim_denoised_SR.nii.gz
  T2=$SUBJDIR/${id}.Screening_MR.High_Res_T2.High-Resolution_2D_Coronal_T2.nii.gz
  ASHST2DIR=$SUBJDIR/ASHST2

  ASHSICVDIR=$SUBJDIR/ASHSICV
  ICV=$(cat $ASHSICVDIR/final/${id}_left_corr_nogray_volumes.txt | cut -d " " -f 5)

  VOL=""
  NSLICE=""
  T2RATERS7T=""
  for side in left right; do
    SEG=$ASHST2DIR/final/${id}_${side}_lfseg_corr_nogray_CortexCleaned.nii.gz
    #SEG=$ASHST2DIR/final/${id}_${side}_lfseg_corr_nogray.nii.gz
    RATER="None"
    if [[ -f $SEG ]]; then
      STATS=$TMPDIR/rawvols_ASHST2_3T.txt
      c3d $SEG -dup -lstat | tee $STATS
      for ((ilab=0; ilab < ${#ASHST23TLABELIDS[*]}; ilab++)); do
        i=${ASHST23TLABELNUMS[ilab]}
        SUB=${ASHST23TLABELIDS[ilab]}
        NSLICE="$NSLICE,$(cat $STATS | awk -v id=$i '$1 == id {print $10}')"
        VOL="$VOL,$(cat $STATS | awk -v id=$i '$1 == id {print $7}')"
      done
    else
      NSLICE="$NSLICE,,,,,,,,,,"
      VOL="$VOL,,,,,,,,,,"
    fi
    T2RATERS7T="$T2RATERS7T,$RATER"
  done
  VOL="${VOL:1}"
  NSLICE="${NSLICE:1}"
  T2RATERS7T="${T2RATERS7T:1}"
  # compute mean
  MEA=$VOL
  NMEA=10
  for ((i=1;i<=$NMEA;i++)); do
    LMEA=$(echo $MEA | cut -d, -f $i)
    RMEA=$(echo $MEA | cut -d, -f $((i+NMEA)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=10;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    MEA="$MEA,$MMEA"
  done
  VOL=$MEA
  T2NSLICE=$NSLICE
  T2VOL=$VOL

  echo "$id,$Type,$ICV,$T2RATERS7T,$T2NSLICE,$T2VOL" > \
    $OUTTMPDIR/$(printf %04d $i)_${id}_${Type}.csv
}

######################################################
function QCASHST2Seg()
{
  # Load id number
  Nbegin=2
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11

  rm -rf $QCDIR
  mkdir -p $QCDIR/png
  echo "ID,TYPE,SIDE,SCANEXIST,ASHST2SEGEXIST,PNGEXIST" > $QCDIR/qc_subj_list.csv


  IDCol=$(csvcol.sh $OUTPUTDIR/checkdata.csv ID)
  TYPECol=$(csvcol.sh $OUTPUTDIR/checkdata.csv AnalysisType)

  # Submit job to copy data
  JPREFIX=SUM
  for ((i=$Nbegin;i<=${N};i++)); do

    ROW=$(cat -A $OUTPUTDIR/checkdata.csv | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # get information
    id=$(echo $ROW | cut -f $IDCol -d ",")
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    PREFIX="Screening"
    SUBJDIR=$DATAROOT/$id

    # perform seg
    TPDIR=$SUBJDIR
    echo "COPY ASHS T2 PNGS: $i,$id,$PRERFIX"
    if [[ $Type != "T2Baseline" && $Type != "T2Longitudinal" ]]; then
      continue
    fi
    T1=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
    T1TRIM=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim.nii.gz
    T1SR=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim_denoised_SR.nii.gz
    T2=$TPDIR/${id}.Screening_MR.High_Res_T2.High-Resolution_2D_Coronal_T2.nii.gz
    ASHST2DIR=$TPDIR/ASHST2

    for side in left right; do
     
      SCANEXIST="FALSE"
      ASHST2SEGEXIST="FALSE"
      if [[ -f $T1TRIM && -f $T2 ]]; then
        SCANEXIST="TRUE"
        if [[ -f $ASHST2DIR/final/${id}_${side}_lfseg_corr_usegray.nii.gz ]]; then
          ASHST2SEGEXIST="TRUE"
        fi
      fi

      PNGEXIST="FALSE"
      if [[ -f $ASHST2DIR/qa/qa_seg_bootstrap_corr_usegray_${side}_qa.png ]]; then
        PNGEXIST="TRUE"
        echo "    COPYING $i,$id,$Type,$side"
        cp $ASHST2DIR/qa/qa_seg_bootstrap_corr_usegray_${side}_qa.png \
           $QCDIR/png/${id}_${Type}_${side}_seg_bootstrap_corr_usegray.png
      else
        echo "    SKIPPING $i,$id,$Type,$side"
      fi

      echo "$id,$Type,$side,$SCANEXIST,$ASHST2SEGEXIST,$PNGEXIST" >> $QCDIR/qc_subj_list.csv

    done

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

