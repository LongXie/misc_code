#//bin/bash
#$ -S /bin/bash
#set +e
set -e -x
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
ROOT=/home/longxie/Janssen
CODEDIR=$ROOT/scripts
ANALYSISDIR=$ROOT/analysis_input
ALLINFO_TXT=$ANALYSISDIR/subj_list_11032019.txt
DATAROOT=/data/jux/wolk_group/janssen/CPSS_MRI_vols_UPenn/UPenn/
OUTPUTDIR=$ROOT/output_runT1_11032019
LOGDIR=$OUTPUTDIR/log
INFODIR=$OUTPUTDIR/info
DUMPDIR=$OUTPUTDIR/dump
mkdir -p $LOGDIR
mkdir -p $INFODIR

QCFILE=$ANALYSISDIR/Janssen_ASHST1_QC.csv

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

ASHST13TLABELIDS=(AHippo PHippo ERC BA35 BA36 PHC)
ASHST13TLABELNUMS=(1     2      10  11   12   13)



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
  # run SR
  ###############################
  #SR

  #ASHST1SegSUB
  #CleanASHST1

  summarize

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
         -q all.q,basic.q \
         -l h_vmem=2.1G,s_vmem=2G \
         -N "${PREFIX}_${i}" \
         $0 CheckData_sub $i $OUTTMPDIR
    sleep 0.03

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # get all info
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  echo "$header,T1Orientation,T1Dimension,AnalysisType"  > $OUTPUTDIR/checkdata.csv
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
  T1NIFTI=$SUBJDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
  
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
  
  ROW="$ROW,$ORI,$DIM,T1Baseline"
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
           -q all.q,basic.q \
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
       -q all.q,basic.q \
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
function SR()
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
  JPREFIX=SR
  for ((i=$Nbegin;i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # get information
    id=$(echo $ROW | cut -f $IDCol -d ",")
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    PREFIX="Screening"
    SUBJDIR=$DATAROOT/$id

    # perform SR
    TPDIR=$SUBJDIR
    echo "SR Processing:$i,$rid,$PREFIX"
    T1=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
    T1TRIM=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim.nii.gz
    T1SR=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim_denoised_SR.nii.gz

    if [[ -f $T1TRIM && ! -f $T1SR ]]; then
      echo "        Adding $rid $PREFIX "
      echo "SR: Adding $rid $PREFIX " >> $LOGDIR/SR_log.txt
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${id}_${PREFIX}" \
           $0 SR_sub \
             $TPDIR \
             $T1TRIM \
             $T1SR
      sleep 0.02
    fi

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}

function SR_sub()
{
  # perform SR
  TPDIR=$1
  T1TRIM=$2
  T1SR=$3

  if [[ ! -f $T1SR ]]; then

    $SRPATH/NLMDenoise \
      -i $T1TRIM \
      -o $TMPDIR/T1w_trim_denoised.nii.gz

    orient_code=$(c3d $T1TRIM -info | cut -d ';' -f 5 | cut -d ' ' -f 5)
    if [[ $orient_code == "Oblique," ]]; then
      orient_code=$(c3d $T1TRIM -info | cut -d ';' -f 5 | cut -d ' ' -f 8)
    fi

    c3d $TMPDIR/T1w_trim_denoised.nii.gz \
      -swapdim RPI \
      -o $TMPDIR/T1w_trim_denoised.nii.gz

    $SRPATH/NLMUpsample \
      -i $TMPDIR/T1w_trim_denoised.nii.gz \
      -o $TMPDIR/T1w_trim_denoised_SR.nii.gz \
      -lf 2 1 2

    c3d $TMPDIR/T1w_trim_denoised_SR.nii.gz\
      -swapdim $orient_code \
      -clip 0 inf \
      -type short \
      -o $T1SR

    echo "Using SRPATH=${SRPATH} to perform denoising and superresolution." > $TPDIR/super_resolution_version.txt

  fi
}

######################################################
function ASHST1SegSUB()
{
  # Load id number
  Nbegin=2
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11
  ATLAS=/home/longxie/ASHS_T1/ASHSexp/exp201/atlas/final
  Nrun=20

  # Submit job to copy data
  JPREFIX=ASHST1
  for ((i=$Nbegin;i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # get information
    id=$(echo $ROW | cut -f $IDCol -d ",")
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    PREFIX="Screening"
    SUBJDIR=$DATAROOT/$id

    # perform seg
    TPDIR=$SUBJDIR
    echo "ASHS T1 Processing:$i,$rid,$PREFIX"
    if [[ $Type != "T1Baseline" && $Type != "T1Longitudinal" ]]; then
      continue
    fi
    T1=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
    T1TRIM=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim.nii.gz
    T1SR=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim_denoised_SR.nii.gz
    ASHST1DIR=$TPDIR/ASHST1
    
    if [[ -f $T1TRIM && -f $SRT1 && \
       ! -f $ASHST1DIR/final/${id}_left_heur_volumes.txt ]]; then

      echo "   Adding $rid $PREFIX "
      echo "ASHS T1: Adding $rid $PREFIX " >> $LOGDIR/ASHST1_baseline_log.txt

      mkdir -p $OUTDIR
      $ASHS_ROOT/bin/ashs_main.sh \
        -a $ATLAS -d -T -I $id -g $T1TRIM -f $SRT1 \
        -Q -s 1-7 \
        -z $CODEDIR/ashs-fast-z.sh \
        -m $CODEDIR/identity.mat -M \
        -w $ASHST1DIR &
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
function CleanASHST1()
{
  # Load id number
  Nbegin=2
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11

  # Submit job to copy data
  JPREFIX=CLASHST1
  for ((i=$Nbegin;i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # get information
    id=$(echo $ROW | cut -f $IDCol -d ",")
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    PREFIX="Screening"
    SUBJDIR=$DATAROOT/$id

    # perform seg
    TPDIR=$SUBJDIR
    echo "ASHS T1 Processing:$i,$rid,$PREFIX"
    if [[ $Type != "T1Baseline" && $Type != "T1Longitudinal" ]]; then
      continue
    fi
    T1=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
    T1TRIM=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim.nii.gz
    T1SR=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim_denoised_SR.nii.gz
    ASHST1DIR=$TPDIR/ASHST1

    if [[ -f $ASHST1DIR/final/${id}_left_lfseg_heur.nii.gz && -f $ASHST1DIR/tse.nii.gz ]]; then

      echo "   Cleaning $rid $PREFIX "
      echo "Cleaning ASHS T1: $rid $PREFIX " >> $LOGDIR/ASHST1_baseline_clean_log.txt

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
    echo "ASHS T1 Processing:$i,$rid,$PREFIX"
    if [[ $Type != "T1Baseline" && $Type != "T1Longitudinal" ]]; then
      continue
    fi
    T1=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted.nii.gz
    T1TRIM=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim.nii.gz
    T1SR=$TPDIR/${id}.Screening_MR.High_Res_T2.3D_T1-weighted_trim_denoised_SR.nii.gz
    T2=$TPDIR/${id}.Screening_MR.High_Res_T2.High-Resolution_2D_Coronal_T2.nii.gz
    ASHST1DIR=$TPDIR/ASHST1
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
  header="$header,L_Modi_Rater_T1_3T,R_Modi_Rater_T1_3T"
  for type in VOL; do
  for side in L R M; do
  for sub in ${ASHST13TLABELIDS[*]}; do
    header="$header,${side}_${sub}_${type}_ASHST1_3T"
  done
  done
  done

  # thickness
  for side in L R M; do
  for type in MeanTHK MedianTHK; do
  for sub in Hippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_MSTTHKMT_ASHST1_3T"
  done
  done
  done

  # thickness QC
  for side in L R; do
  for sub in Hippo ERC BA35 BA36 PHC ALL; do
    header="$header,${side}_${sub}_fitquality_MSTTHKMT_ASHST1_3T"
  done
  done

  for side in L R; do
    header="$header,${side}_QA_Hippo_ASHST1_3T,${side}_QA_ERC_ASHST1_3T,${side}_QA_BA35_ASHST1_3T,${side}_QA_BA36_ASHST1_3T,${side}_QA_PHC_ASHST1_3T"
  done

  echo $header > $OUTPUTDIR/Janssen_ASHST1_measurements.csv
  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/Janssen_ASHST1_measurements.csv
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
  ASHST1DIR=$SUBJDIR/ASHST1
  ASHST2DIR=$SUBJDIR/ASHST2

  ASHSICVDIR=$SUBJDIR/ASHSICV
  ICV=$(cat $ASHSICVDIR/final/${id}_left_corr_nogray_volumes.txt | cut -d " " -f 5)

  VOL=""
  T1RATERS3T=""
  for side in left right; do
    SEG=$ASHST1DIR/final/${id}_${side}_lfseg_heur.nii.gz
    RATER="None"
    if [[ -f $SEG ]]; then
      STATS=$TMPDIR/rawvols_ASHST1_3T.txt
      c3d $SEG -dup -lstat | tee $STATS
      for ((ilab=0; ilab < ${#ASHST13TLABELIDS[*]}; ilab++)); do
        i=${ASHST13TLABELNUMS[ilab]}
        SUB=${ASHST13TLABELIDS[ilab]}
        VOL="$VOL,$(cat $STATS | awk -v id=$i '$1 == id {print $7}')"
      done
    else
      VOL="$VOL,,,,,,"
    fi
    T1RATERS3T="$T1RATERS3T,$RATER"
  done
  VOL="${VOL:1}"
  T1RATERS3T="${T1RATERS3T:1}"
  # compute mean
  MEA=$VOL
  NMEA=6
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
  T1VOL3T=$VOL

  # get MTL thickness
  THK=""
  MSTTHKASHST1DIR=$SUBJDIR/MSTThkASHST1
  for side in left right; do
    THKCSV=$MSTTHKASHST1DIR/${id}_${side}_thickness.csv
    THK="$THK,$(cat $THKCSV | grep MultiTemp \
      | cut -d , -f 5-14)"
  done
  THK="${THK:1}"
  # compute mean
  MEA=$THK
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
  THK=$MEA

  FITQUAL=""
  for side in left right; do
    THKCSV=$MSTTHKASHST1DIR/${id}_${side}_thickness.csv
    FITQUAL="$FITQUAL,$(cat $THKCSV | grep MultiTemp \
      | cut -d , -f 15-20)"
  done
  FITQUAL="${FITQUAL:1}"

  #####################################################################
  # ASHS T2 3T QC
  QCASHST13T=""
  set +e
  for side in left right; do
    QCROW=$(cat -A $QCFILE | grep $id | grep $side | sed -e "s/\^M\\$//g")
    QCASHST13T="$QCASHST13T,$(echo $QCROW | cut -d, -f 3)"
    QCASHST13T="$QCASHST13T,$(echo $QCROW | cut -d, -f 4)"
    QCASHST13T="$QCASHST13T,$(echo $QCROW | cut -d, -f 5)"
    QCASHST13T="$QCASHST13T,$(echo $QCROW | cut -d, -f 6)"
    QCASHST13T="$QCASHST13T,$(echo $QCROW | cut -d, -f 7)"
  done
  QCASHST13T="${QCASHST13T:1}"
  set -e


  echo "$id,$Type,$ICV,$T1RATERS3T,$T1VOL3T,$THK,$FITQUAL,$QCASHST13T" > \
    $OUTTMPDIR/$(printf %04d $i)_${id}_${Type}.csv
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

