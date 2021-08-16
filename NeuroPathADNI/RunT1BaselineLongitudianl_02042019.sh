#//bin/bash
#$ -S /bin/bash
set -e

set -e -x
#set -x -e

##############################################
# Setup environment
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
MATLAB_BIN=/share/apps/matlab/R2017a/bin/matlab
export ALOHA_ROOT=/data/picsl/longxie/pkg/aloha
SRMATLABCODEDIR=/home/longxie/ASHS_PHC/SuperResolution/SRToolBox/matlabfunction
SHOOTWPRIORDIR=/home/longxie/LMTemplateMatching/PointSetGeodesicShooting/build_release
LMTOWARPDIR=/home/longxie/LMTemplateMatching/PointSetUtilities/build_release
export FREESURFER_HOME=$CFNAPPS/freesurfer/6.0.0
source $FREESURFER_HOME/SetUpFreeSurfer.sh
SRPATH=/home/longxie/pkg/PatchSuperResolution/build_release

##############################################
# Directories
ROOT=/home/longxie/NeuroPathADNI/
CODEDIR=$ROOT/scripts
ANALYSISDIR=$ROOT/analysis_input
ALLINFO_TXT=$ROOT/output_organize/NeuroPathScanInfo_04112019.csv
ALLDATADIR=$ROOT/NeuroPath_dataset
OUTPUTDIR=$ROOT/output_runall_04112019
LOGDIR=$OUTPUTDIR/log
DUMPDIR=$OUTPUTDIR/dump
mkdir -p $LOGDIR

#RAWNIFTIDIR=$ROOT/rawdata/nifti
#DATADIR=$ROOT/dataset_local
#ASHSTEMPLATEDIR=/home/longxie/ASHS_T1/ASHSexp/exp008/headtailatlas/final
#BLINFODIR=$ROOT/info/T1Baseline
#SUBJ_TXT=$BLINFODIR/subjectlist_04072018.csv
#BLOUTPUTDIR=$ROOT/output_all
#BLTXT=$BLOUTPUTDIR/ADNI_T1_measurement_withFS.csv
#ALLINFO_TXT=$ANALYSISDIR/ADNI2018_baseline_subjectlist_04302018.csv
DEMOG_TXT=$ANALYSISDIR/NEUROPATH_04_12_18.csv

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi



#############################################
# parameters
#ylthresh="0.0"
#yhthresh="2.5"
#minfwtime="1.2"

#STATSDIR=$OUTPUTDIR/stats_${ylthresh}-${yhthresh}yrs_minfw${minfwtime}yrs
#DESIGNDIR=$STATSDIR/design

#GSTEMPDIR=/home/longxie/ASHS_T1/thickness/exp/exp601/

#MTSTATSDIR=$OUTPUTDIR/stats_${ylthresh}-${yhthresh}yrs_minfw${minfwtime}yrs_multitemp
#MTDESIGNDIR=$MTSTATSDIR/design

#UTTEMPDIR=/home/longxie/ASHS_T1/thickness/exp/exp510/FinalTemp/template3

#############################################
function main()
{
  reset_dir

  ################################
  # link and summarize data
  ################################
  #LinkData

  ###############################
  # Run trim neck
  ###############################
  #TrimNeck

  ###############################
  # run SR
  ###############################
  #SR


  ###############################
  # process baseline scans

  # perform segmentation on baseline scans
  #ASHSSegSUB
  #CleanASHS

  # perform ICV segmentagion
  #ASHSICVSegSUB
  #CleanASHSICV


  ###############################
  # perform longitudinal analysis
  ###############################
  # run aloha
  #AlohaMTL

  # perform QC
  #ComputePWNCC
  #SelectCasesToQC
  
  # rerun ALOHA for subjects whose segmentations are manually corrected
  #AlohaMTL_ManSegCorrected

  ###############################
  # summary all information analysis
  ###############################
  # summarize result
  Summarize

  # summarize result for Robin whole brain covariance analysis
  #SummarizeRobin








  ###############################
  # perform pointwise longitudinal analysis
  ###############################
  #MTLCortexRegionalThk

  #ComputeRegionalAtropyRate

  #RegionalStats_meshglm

  ###############################
  # perform pointwise longitudinal analysis using the unified template method
  ###############################
  #MTLCortexRegionalThkMultiTemp

  #ComputeRegionalAtropyRateMultiTemp

  #UTRegionalStats_meshglm

  ###############################
  # summary thickness analysis
  ###############################
  # summarize result
  #SummarizeMultiTemp

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
function LinkData()
{
  # Load id number
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=20

  # initialize tmp dir
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # Submit job to copy data
  PREFIX=LD${expid}
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  for ((i=2;i<=${N};i++)); do

    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=2.1G,s_vmem=2G \
         -N "${PREFIX}_${i}" \
         $0 LinkData_sub $i $OUTTMPDIR
    sleep 0.1

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # get all info
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  echo "$header,LinkExist,T1Orientation,T1Dimension,T1Resolution,LinkMatch,LinkT1NIFTI,AnalysisType,END"  > $OUTPUTDIR/checkdata.csv
  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/checkdata.csv
}

function LinkData_sub()
{
  i=$1
  OUTTMPDIR=$2
  ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  
  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  T1NIFTICol=$(csvcol.sh $ALLINFO_TXT FINALT1NIFTI)
  DATEDIFFTOBLCol=$(csvcol.sh $ALLINFO_TXT DATEDIFFTOBL)

  # get information
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX1="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
  blscandate=$(ReFormateDate $blscandate)
  T1NIFTI=$(echo $ROW | cut -f $T1NIFTICol -d ",")
  datediff=$(echo $ROW | cut -f $DATEDIFFTOBLCol -d ",")
  SUBJALLDIR=$ALLDATADIR/$id
  #SUBJDIR=$DATADIR/$id

  # if file does not exist, skip
  if [[ $T1NIFTI == "" || $T1NIFTI == "None"  ]]; then

    ROW="$ROW,,,,,,,,-1"

  else

    # link scan
    mkdir -p $SUBJALLDIR/${PREFIX1}
    T1SCAN=$SUBJALLDIR/${PREFIX1}/${PREFIX1}_${id}_T1w.nii.gz
    if [[ ! -f $T1SCAN ]]; then
      ln -sf $T1NIFTI $T1SCAN
      sleep 0.1
    fi
  
    # get scan info
    if [[ -f $T1SCAN ]]; then

      LINKNIFTI=$(readlink -sf $T1SCAN)
      if [[ $LINKNIFTI == $T1NIFTI ]]; then
        LINKMATCH=1
      else
        LINKMATCH=0
      fi

      EXIST=1
      ORI=$(c3d $T1SCAN -info \
              | cut -d ';' -f 5 | cut -d ' ' -f 5)
      if [[ $ORI == "Oblique," ]]; then
        ORI=$(c3d $T1SCAN -info \
                | cut -d ';' -f 5 | cut -d ' ' -f 8)
      fi
      DIM=$(c3d $T1SCAN -info \
        | cut -d ':' -f 2 | cut -d ';' -f 1 | sed -e 's/,/_/g')
      RES=$(c3d $T1SCAN -info \
        | cut -d ';' -f 3 | sed -e 's/,/_/g')

    else

      LINKNIFTI=""
      LINKMATCH=""
      EXIST=0
      ORI=""
      DIM=""
      RES=""

    fi

    if [[ $datediff == "0" ]]; then
      ROW="$ROW,$EXIST,$ORI,$DIM,$RES,$LINKMATCH,$LINKNIFTI,T1Baseline,-1"
    else
      ROW="$ROW,$EXIST,$ORI,$DIM,$RES,$LINKMATCH,$LINKNIFTI,T1Longitudinal,-1" 
    fi
  fi

  echo $ROW > $OUTTMPDIR/$(printf %04d $i)_${id}.csv
}

######################################################
function TrimNeck()
{
  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11

  # Submit job to copy data
  PREFIX1=TN${expid}
  for ((i=2;i<=${N};i++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # columns
    RIDCol=1
    IDCol=$(csvcol.sh $ALLINFO_TXT ID)
    ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
    TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    #if [[ $Type == "T1Baseline" ]]; then
    #  SUBJDIR=$ALLDATADIR/$id
    #elif [[ $Type == "T1Longitudinal" ]]; then
    #  SUBJDIR=$DATADIR/$id
    #else
    #  continue
    #fi

    # perform trim neck
    TPDIR=$SUBJALLDIR/${PREFIX}
    echo "TrimNeck Processing:$i,$rid,$PREFIX"
    #echo $TPDIR/${PREFIX}_${id}_T1w.nii.gz

    if [[ -f $TPDIR/${PREFIX}_${id}_T1w.nii.gz && \
       ! -f $TPDIR/${PREFIX}_${id}_T1w_trim.nii.gz ]]; then
      echo "        Adding $rid $PREFIX "
      echo "TrimNeck: Adding $rid $PREFIX " >> $LOGDIR/trimneck_log.txt
      qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX1}_${id}_${PREFIX}" \
           $0 TrimNeck_sub \
             $TPDIR \
             $TPDIR/${PREFIX}_${id}_T1w.nii.gz \
             $TPDIR/${PREFIX}_${id}_T1w_trim.nii.gz
      sleep 0.02
    fi

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX1}_*" -sync y -b y \
       sleep 1
}

######################################################
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
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11

  # Submit job to copy data
  PREFIX1=SR${expid}
  for ((i=2;i<=${N};i++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # columns
    RIDCol=1
    IDCol=$(csvcol.sh $ALLINFO_TXT ID)
    ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
    TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id


    #if [[ $Type == "T1Baseline" ]]; then
    #  SUBJDIR=$ALLDATADIR/$id
    #elif [[ $Type == "T1Longitudinal" ]]; then
    #  SUBJDIR=$DATADIR/$id
    #fi

    # perform SR
    TPDIR=$SUBJALLDIR/${PREFIX}
    #BLDIR=$ALLDATADIR/$id/${blscandate}
    echo "SR Processing: $i,$rid,$PREFIX"


    if [[ -f $TPDIR/${PREFIX}_${id}_T1w_trim.nii.gz && ! -f $TPDIR/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz  ]]; then

      echo "   Adding $rid $PREFIX "
      echo "SR: Adding $rid $PREFIX " >> $LOGDIR/SR_log.txt
      qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX1}_${id}_${PREFIX}" \
           $0 SR_sub $TPDIR \
              $TPDIR/${PREFIX}_${id}_T1w_trim.nii.gz \
              $TPDIR/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
      sleep 0.02
      #exit
    fi

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX1}_*" -sync y -b y \
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
      -o $T1SR

    echo "Using SRPATH=${SRPATH} to perform denoising and superresolution." > $TPDIR/super_resolution_version.txt

  fi
}

######################################################
function ASHSSegSUB()
{
  # Load id number
  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11
  ATLAS=/home/longxie/ASHS_T1/ASHSexp/exp201/atlas/final
  Nrun=20

  # Submit job to copy data
  PREFIX1=SR${expid}
  for ((i=2;i<=${N};i++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # columns
    RIDCol=1
    IDCol=$(csvcol.sh $ALLINFO_TXT ID)
    ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
    TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    echo "ASHS T1 Processing: $i,$rid,$PREFIX"
    if [[ $Type != "T1Baseline" ]]; then
      continue
    fi

    # perform SR
    TPDIR=$SUBJALLDIR/${PREFIX}
    OUTDIR=$TPDIR/ASHST1
    T1=$TPDIR/${PREFIX}_${id}_T1w_trim.nii.gz
    SRT1=$TPDIR/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz

    if [[ -f $T1 && -f $SRT1 && \
       ! -f $OUTDIR/final/${id}_left_heur_volumes.txt ]]; then

      echo "   Adding $rid $PREFIX "
      echo "ASHS T1: Adding $rid $PREFIX " >> $LOGDIR/ASHST1_baseline_log.txt  

      mkdir -p $OUTDIR
      $ASHS_ROOT/bin/ashs_main.sh \
        -a $ATLAS -d -T -I $id -g $T1 -f $SRT1 \
        -Q -s 1-7 \
        -z $CODEDIR/ashs-fast-z.sh \
        -m $CODEDIR/identity.mat -M \
        -w $OUTDIR &
      sleep 0.1

      ALLJOBS=($(jobs -p))
      while [ ${#ALLJOBS[*]} -ge $Nrun ]; do
        sleep 60
        ALLJOBS=($(jobs -p))
      done
    fi

  done
}

##################################################
function CleanASHS()
{
  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)

  # Submit job to copy data
  PREFIX=ASHSSeg${expid}
  for ((i=2;i<=${N};i++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # columns
    RIDCol=1
    IDCol=$(csvcol.sh $ALLINFO_TXT ID)
    ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
    TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    echo "Checking ASHS T1: $i,$rid,$PREFIX"
    if [[ $Type != "T1Baseline" ]]; then
      continue
    fi
    ASHST1DIR=$SUBJALLDIR/${PREFIX}/ASHST1

    #echo "$ASHST1DIR/final/${id}_left_lfseg_heur.nii.gz"
    #echo "$ASHST1DIR/tse.nii.gz"

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

##################################################
function ASHSICVSegSUB()
{
  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11
  ICVATLAS=/data/jet/vpiskin/ASHS/build_ICV_atlas_trimmed/ICV_atlas/final/
  Nrun=20

  # Submit job to copy data
  PREFIX=ASHSSeg${expid}
  for ((i=2;i<=${N};i++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # columns
    RIDCol=1
    IDCol=$(csvcol.sh $ALLINFO_TXT ID)
    ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
    TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    echo "ASHS ICV Processing: $i,$rid,$PREFIX"
    if [[ $Type != "T1Baseline" ]]; then
      continue
    fi

    # perform SR
    TPDIR=$SUBJALLDIR/${PREFIX}
    OUTDIR=$TPDIR/ASHSICV
    T1=$TPDIR/${PREFIX}_${id}_T1w_trim.nii.gz

    if [[ -f $T1 && ! -f $OUTDIR/final/${id}_left_corr_nogray_volumes.txt ]]; then

      echo "   Adding $rid $PREFIX "
      echo "ASHS ICV: Adding $rid $PREFIX " >> $LOGDIR/ASHSICV_baseline_log.txt

      mkdir -p $OUTDIR
      $ASHS_ROOT/bin/ashs_main.sh \
        -a $ICVATLAS -d -T -I $id -g $T1 -f $T1 \
        -Q -s 1-7 \
        -z $CODEDIR/ashs-fast-z.sh \
        -m $CODEDIR/identity.mat -M \
        -B \
        -w $OUTDIR &
      sleep 0.1

      ALLJOBS=($(jobs -p))
      while [ ${#ALLJOBS[*]} -ge $Nrun ]; do
        sleep 60
        ALLJOBS=($(jobs -p))
      done

    fi

  done
}

##################################################
function CleanASHSICV()
{
  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)

  # Submit job to copy data
  PREFIX=ASHSSeg${expid}
  for ((i=2;i<=${N};i++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # columns
    RIDCol=1
    IDCol=$(csvcol.sh $ALLINFO_TXT ID)
    ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
    TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    echo "Checking ASHS ICV: $i,$rid,$PREFIX"
    if [[ $Type != "T1Baseline" ]]; then
      continue
    fi
    ASHSICVDIR=$SUBJALLDIR/${PREFIX}/ASHSICV

    #echo "$ASHST1DIR/final/${id}_left_lfseg_heur.nii.gz"
    #echo "$ASHST1DIR/tse.nii.gz"

    if [[ -f $ASHSICVDIR/final/${id}_left_lfseg_heur.nii.gz && -f $ASHSICVDIR/tse.nii.gz ]]; then

      echo "   Cleaning $rid $PREFIX "
      echo "Cleaning ASHS ICV: $rid $PREFIX " >> $LOGDIR/ASHSICV_baseline_clean_log.txt

      rm -rf $ASHSICVDIR/affine_t1_to_template \
        $ASHSICVDIR/ants_t1_to_temp \
        $ASHSICVDIR/bootstrap \
        $ASHSICVDIR/dump \
        $ASHSICVDIR/multiatlas \
        $ASHSICVDIR/mprage* \
        $ASHSICVDIR/tse.nii.gz \
        $ASHSICVDIR/tse_to_chunktemp*.nii.gz \
        $ASHSICVDIR/flirt_t2_to_t1 \
        $ASHSICVDIR/tse_native_chunk*.nii.gz \
        $ASHSICVDIR/qa

    fi

  done
}


######################################################
function AlohaMTL()
{
  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)  
  N_begin=2

  # Submit job to copy data
  JPREFIX=ALOHA${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # columns
    RIDCol=1
    IDCol=$(csvcol.sh $ALLINFO_TXT ID)
    ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
    TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    #SUBJDIR=$DATADIR/$id

    # skip if it is not longitudinal
    echo "ALOHA Processing: $i,$rid,$PREFIX"
    if [[ $Type != "T1Longitudinal" ]]; then
      continue
    fi

    # perform SR
    FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz

    # baseline information
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    DD=$(date -d "$blscandate" '+%d')
    MM=$(date -d "$blscandate" '+%m')
    YYYY=$(date -d "$blscandate" '+%Y')
    BLPREFIX="${YYYY}-${MM}-${DD}"
    BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
    BLT1SRIMG=$BLSUBJDIR/${BLPREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
    BLASHST1DIR=$BLSUBJDIR/ASHST1
    BLASHST1LSEG=$BLASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
    BLASHST1RSEG=$BLASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

    # start running ALOHA
    SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
    if [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $BLASHST1LSEG && -f $BLASHST1RSEG && ! -d $SUBJALOHADIR ]]; then

        echo "   Adding $rid $PREFIX"
        echo "ALOHA: Adding $rid $PREFIX" >> $LOGDIR/ALOHA_log.txt

        qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q,basic.q \
           -l h_vmem=5.1G,s_vmem=5G \
           -N "${JPREFIX}_${i}_${rid}_${PREFIX}" \
           $0 AlohaMTL_sub \
               $BLT1SRIMG \
               $FUT1SRIMG \
               $BLASHST1LSEG \
               $BLASHST1RSEG \
               $SUBJALOHADIR
          sleep 0.1

      fi

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}

function AlohaMTL_sub()
{
  BLT1SRIMG=$1
  YT1SRIMG=$2
  BLASHST1LSEG=$3
  BLASHST1RSEG=$4
  SUBJALOHADIR=$5
  mkdir -p $SUBJALOHADIR

  if [[ ! -f $SUBJALOHADIR/results/volumes_left.txt || ! -f $SUBJALOHADIR/results/volumes_right.txt ]]; then

  c3d $BLASHST1LSEG \
     -replace 1 999 2 999 10 999 11 999 12 999 13 999 \
     -thresh 999 999 1 0 \
     -o $SUBJALOHADIR/bl_seg_left.nii.gz

  c3d $BLASHST1RSEG \
     -replace 1 999 2 999 10 999 11 999 12 999 13 999 \
     -thresh 999 999 1 0 \
     -o $SUBJALOHADIR/bl_seg_right.nii.gz

  $ALOHA_ROOT/scripts/aloha_main.sh \
    -b $BLT1SRIMG \
    -f $YT1SRIMG \
    -r $SUBJALOHADIR/bl_seg_left.nii.gz \
    -s $SUBJALOHADIR/bl_seg_right.nii.gz \
    -w $SUBJALOHADIR \
    -t 1-4

  fi

  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do 

    if [[ ! -f $SUBJALOHADIR/results_${sub}/volumes_left.txt || ! -f $SUBJALOHADIR/results_${sub}/volumes_right.txt ]]; then

      if [[ $sub == "ERC" ]]; then
        label=(10 10)
      elif [[ $sub == "BA35" ]]; then
        label=(11 11)
      elif [[ $sub == "BA36" ]]; then
        label=(12 12)
      elif [[ $sub == "PHC" ]]; then
        label=(13 13)
      elif [[ $sub == "AHippo" ]]; then
        label=(1 1)
      elif [[ $sub == "PHippo" ]]; then
        label=(2 2)
      elif [[ $sub == "Hippo" ]]; then
        label=(1 2)
      fi

      SUBJALOHASUBDIR=$TMPDIR/aloha_${sub}
      mkdir -p $SUBJALOHASUBDIR/results
      ln -sf $SUBJALOHADIR/deformable $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/init $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/global $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/dump $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/final $SUBJALOHASUBDIR
      
      c3d $BLASHST1LSEG \
        -thresh ${label[0]} ${label[1]} 1 0 \
        -o $SUBJALOHASUBDIR/bl_seg_left.nii.gz

      c3d $BLASHST1RSEG \
        -thresh ${label[0]} ${label[1]} 1 0 \
        -o $SUBJALOHASUBDIR/bl_seg_right.nii.gz

      $ALOHA_ROOT/scripts/aloha_main.sh \
        -b $BLT1SRIMG \
        -f $YT1SRIMG \
        -r $SUBJALOHASUBDIR/bl_seg_left.nii.gz \
        -s $SUBJALOHASUBDIR/bl_seg_right.nii.gz \
        -w $SUBJALOHASUBDIR \
        -t 4
      
      mv $SUBJALOHASUBDIR/results \
         $SUBJALOHADIR/results_${sub}

    fi
  
    # measure thickness
    for side in left right; do

      if [[ ! -f $SUBJALOHADIR/results_${sub}/mean_thickness_${side}.txt || ! -f $SUBJALOHADIR/results_${sub}/median_thickness_${side}.txt ]]; then

        # measure thickness of bl mesh
        cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
          -T $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_thickmap.vtk \
          -p 1.2 -e 6 \
          $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}.vtk \
          $SUBJALOHADIR/results_${sub}/skel_blmptrim_seg_${side}.vtk
       
        # measure thickness of fu mesh
        cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
          -T $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_warped_to_futrim_om_thickmap.vtk \
          -p 1.2 -e 6 \
          $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_warped_to_futrim_om.vtk \
          $SUBJALOHADIR/results_${sub}/skel_blmptrim_seg_${side}_warped_to_futrim_om.vtk

        # get thickness
        $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
        addpath('/home/longxie/ASHS_T1/Application/TAUPET/longitudinal/code');
        MeasureMeanMedianThickness('$SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_thickmap.vtk','$SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_warped_to_futrim_om_thickmap.vtk','$SUBJALOHADIR/results_${sub}/mean_thickness_${side}.txt','$SUBJALOHADIR/results_${sub}/median_thickness_${side}.txt');
MATCODE

     fi
   done
  done

  # clean up
  #cp $SUBJALOHADIR/deformable/mprage_global_long_*_omRAS_half.mat $SUBJALOHADIR/
  #cp $SUBJALOHADIR/deformable/mp_antsreg3d_*Warp*vec.nii.gz $SUBJALOHADIR/
  rm -rf $SUBJALOHADIR/final $SUBJALOHADIR/global \
         $SUBJALOHADIR/init \
         $SUBJALOHADIR/dump  $SUBJALOHADIR/bl_seg_*.nii.gz
  #mkdir -p $SUBJALOHADIR/deformable/
  #mv $SUBJALOHADIR/*.mat $SUBJALOHADIR/deformable/
  #mv $SUBJALOHADIR/*.nii.gz $SUBJALOHADIR/deformable/
}

######################################################
function AlohaMTL_ManSegCorrected()
{
  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=2
  N_begin=2

  # Submit job to copy data
  JPREFIX=ALOHARERUN
  for ((i=${N_begin};i<=${N};i++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # columns
    RIDCol=1
    IDCol=$(csvcol.sh $ALLINFO_TXT ID)
    ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
    TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    #SUBJDIR=$DATADIR/$id

    # perform SR
    FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz

    # baseline information
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    DD=$(date -d "$blscandate" '+%d')
    MM=$(date -d "$blscandate" '+%m')
    YYYY=$(date -d "$blscandate" '+%Y')
    BLPREFIX="${YYYY}-${MM}-${DD}"
    BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
    BLT1SRIMG=$BLSUBJDIR/${BLPREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
    BLASHST1DIR=$BLSUBJDIR/ASHST1
    BLASHST1LSEG=$BLASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
    BLASHST1RSEG=$BLASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

    # Modified segmentations by Laura. If exist, rerun aloha
    RERUN=0
    LWBLASHST1LSEG=/data/picsl/lwisse/ADNIpath/${id}_left_lfseg_heurLW.nii.gz
    if [[ -f $LWBLASHST1LSEG ]]; then
      cp $LWBLASHST1LSEG $BLASHST1DIR/final/
      BLASHST1LSEG=$BLASHST1DIR/final/${id}_left_lfseg_heurLW.nii.gz
      RERUN=1
    fi
    LWBLASHST1RSEG=/data/picsl/lwisse/ADNIpath/${id}_right_lfseg_heurLW.nii.gz
    if [[ -f $LWBLASHST1RSEG ]]; then
      cp $LWBLASHST1RSEG $BLASHST1DIR/final/
      BLASHST1RSEG=$BLASHST1DIR/final/${id}_right_lfseg_heurLW.nii.gz
      RERUN=1
    fi

    # skip if it is not longitudinal
    echo "ALOHA Rerun Processing: $i,$rid,$PREFIX"
    if [[ $Type != "T1Longitudinal" ]]; then
      continue
    fi


    # start running ALOHA
    ORIGSUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
    SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}_ManSegCorrected
    if [[ $RERUN == 1 && -f $FUT1SRIMG && -f $BLT1SRIMG && -f $BLASHST1LSEG && -f $BLASHST1RSEG && ! -d $SUBJALOHADIR ]]; then

        echo "   Adding $id $PREFIX"
        echo "ALOHA Rerun: Adding $id $PREFIX" >> $LOGDIR/ALOHA_rerun_log.txt

        qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q,basic.q \
           -l h_vmem=5.1G,s_vmem=5G \
           -N "${JPREFIX}_${i}_${rid}_${PREFIX}" \
           $0 AlohaMTL_sub \
               $BLT1SRIMG \
               $FUT1SRIMG \
               $BLASHST1LSEG \
               $BLASHST1RSEG \
               $SUBJALOHADIR \
               $ORIGSUBJALOHADIR
          sleep 0.1

      fi

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}

######################################################
function ComputePWNCC()
{
  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=20
  N_begin=2

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # header
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
  echo "$header,ALOHA_Success,NCOR_L,NCOR_R" > $OUTPUTDIR/QC_ALOHA.csv

  # Submit job to copy data
  JPREFIX=QCALOHA${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # columns
    RIDCol=1
    IDCol=$(csvcol.sh $ALLINFO_TXT ID)
    ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
    TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    SUBJDIR=$DATADIR/$id

    # skip if it is not longitudinal
    echo "ALOHA QC Processing: $i,$rid,$PREFIX"
    if [[ $Type != "T1Longitudinal" ]]; then
      echo "$ROW,-1,," > $OUTTMPDIR/$(printf %04d $i).csv
      continue
    fi

    # perform SR
    FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz

    # baseline information
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    DD=$(date -d "$blscandate" '+%d')
    MM=$(date -d "$blscandate" '+%m')
    YYYY=$(date -d "$blscandate" '+%Y')
    BLPREFIX="${YYYY}-${MM}-${DD}"
    BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
    BLT1SRIMG=$BLSUBJDIR/${BLPREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
    BLASHST1DIR=$BLSUBJDIR/ASHST1
    BLASHST1LSEG=$BLASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
    BLASHST1RSEG=$BLASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

    # start running ALOHA
    SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
    #echo $FUT1SRIMG
    #echo $BLT1SRIMG
    #echo $BLASHST1LSEG
    #echo $BLASHST1RSEG

    if [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $BLASHST1LSEG && -f $BLASHST1RSEG ]]; then

      echo "   Adding $rid $PREFIX"
      echo "ALOHA QC: Adding $rid $PREFIX" >> $LOGDIR/ALOHA_QC_log.txt

        qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${i}_${rid}_${PREFIX}" $0 \
               ComputePWNCC_sub \
               $i $SUBJALOHADIR $OUTTMPDIR
          sleep 0.1

    else

      echo "$ROW,-1,," > $OUTTMPDIR/$(printf %04d $i).csv

    fi
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1

  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/QC_ALOHA.csv

}


function ComputePWNCC_sub()
{
  i=$1
  SUBJALOHADIR=$2
  OUTTMPDIR=$3
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv 
  ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g') 
  echo "$ROW,0,," > $OUTTMPDIR/$(printf %04d $i).csv
 
  # get measure
  NCORALL=""
  for side in left right; do

    # get mask
    mesh2img -f \
      -vtk $SUBJALOHADIR/results/blmptrim_seg_${side}_tohw.vtk \
      -a 0.3 0.3 0.3 4 \
      $TMPDIR/blmptrim_seg_${side}_tohw.nii.gz

    c3d $SUBJALOHADIR/deformable/blmptrim_${side}_to_hw.nii.gz \
      $TMPDIR/blmptrim_seg_${side}_tohw.nii.gz \
      -int 0 -reslice-identity \
      -o $SUBJALOHADIR/results/blmptrim_seg_${side}_tohw.nii.gz\
      -dilate 1 20x20x20vox \
      -o $SUBJALOHADIR/results/blmptrim_seg_${side}_tohw_mask20vox.nii.gz

    NCOR=$(c3d $SUBJALOHADIR/deformable/blmptrim_${side}_to_hw.nii.gz \
      $SUBJALOHADIR/results/blmptrim_seg_${side}_tohw_mask20vox.nii.gz \
      -multiply \
      $SUBJALOHADIR/deformable/fumptrim_om_${side}to_hw.nii.gz \
      $SUBJALOHADIR/results/blmptrim_seg_${side}_tohw_mask20vox.nii.gz \
      -multiply \
      -ncor | cut -d = -f 2)

    echo $NCOR > $SUBJALOHADIR/qc_ncor_${side}.txt
    NCORALL="$NCORALL,$NCOR"

  done

  echo "$ROW,1$NCORALL" > $OUTTMPDIR/$(printf %04d $i).csv 
}

######################################################
function SelectCasesToQC()
{
  # get the selection first
  INCSV=$OUTPUTDIR/QC_ALOHA.csv
  OUTCSV=$OUTPUTDIR/QC_ALOHA_selected.csv
  if [[ ! -f $OUTCSV ]]; then
    $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('/home/longxie/ADNI2018/matlabcode');
    ALOHA_QC_Select('$INCSV', '$OUTCSV', 0.8, 0.1);
MATCODE
  fi

  # Load id number
  ALLINFO_TXT=$OUTCSV
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=2
  N_begin=2

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)
  NCORLCol=$(csvcol.sh $ALLINFO_TXT NCOR_L)
  NCORRCol=$(csvcol.sh $ALLINFO_TXT NCOR_R)
  QCSELLCol=$(csvcol.sh $ALLINFO_TXT QC_SEL_L)
  QCSELRCol=$(csvcol.sh $ALLINFO_TXT QC_SEL_R)

  for ((i=${N_begin};i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id

    # perform SR
    FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz

    # baseline information
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    DD=$(date -d "$blscandate" '+%d')
    MM=$(date -d "$blscandate" '+%m')
    YYYY=$(date -d "$blscandate" '+%Y')
    BLPREFIX="${YYYY}-${MM}-${DD}"
    BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
    BLT1SRIMG=$BLSUBJDIR/${BLPREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
    BLASHST1DIR=$BLSUBJDIR/ASHST1
    BLASHST1LSEG=$BLASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
    BLASHST1RSEG=$BLASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

    # ALOHA directory
    SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}

    for side in left right; do 

      echo "ALOHA QC Processing: $i, $id, $side, $scandate"
      if [[ $side == "left" ]]; then
        qcsel=$(echo $ROW | cut -f $QCSELLCol -d ",")
        ncor=$(echo $ROW | cut -f $NCORLCol -d ",")
      else
        qcsel=$(echo $ROW | cut -f $QCSELRCol -d ",")
        ncor=$(echo $ROW | cut -f $NCORRCol -d ",")
      fi

      # visualize when it is selected
      if [[ $qcsel -gt 0 ]]; then

        echo "       ncor = $ncor, selection type = $qcsel"
        echo ""

        # after rigid
        /home/srdas/bin/rview \
          $SUBJALOHADIR/deformable/blmptrim_${side}_to_hw.nii.gz \
          $SUBJALOHADIR/deformable/fumptrim_om_${side}to_hw.nii.gz \
          -res 6 &

        # after deformable
        /home/srdas/bin/rview \
          $SUBJALOHADIR/deformable/blmptrim_${side}_to_hw.nii.gz \
          $SUBJALOHADIR/deformable/fumptrim_om_to_hw_warped_3d_${side}.nii.gz \
          -res 6

      fi
    done
  done  
}

######################################################
function Summarize()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  QC_TXT=$OUTPUTDIR/QC_ALOHA_ADNIpath20190507.csv
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=10

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # Submit job to copy data
  if [[ 1 == 1 ]]; then
  PREFIX1=SM
  for ((i=2;i<=${N};i++)); do

    echo "$i submitted"
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q,gpu.q,himem.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX1}_${i}" $0 \
         Summarize_sub $i $OUTTMPDIR
    sleep 0.05

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 60

  fi

  # get all info
  # header
  #blheader=$(cat -A $DEMOG_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  #sed -e 's/\^M\$//g'
  header="$header,BaselineICV,ALOHA_Success,DateDiffFromBLScanDate,EXCL_L,EXCL_R,EXCL_MRI"
  longheader=""
  for type in VOL MeanTHK MedianTHK; do
  for side in L R; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  done
  done
  done
  #for type in MeanTHKFITTED MedianTHKFITTED; do
  #for side in L R; do
  #for sub in ERC BA35 BA36 PHC; do
  #  header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
  #  longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  #done
  #done
  #done
  header="$header$longheader,END"
  echo $header > $OUTPUTDIR/longitudinal.csv

  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/longitudinal.csv

}

function Summarize_sub()
{
  i=$1
  OUTTMPDIR=$2
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

  # get information
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  Type=$(echo $ROW | cut -f $TYPECol -d ",")
  SUBJALLDIR=$ALLDATADIR/$id
  FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz

  # baseline information
  blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
  blscandate=$(ReFormateDate $blscandate)
  DD=$(date -d "$blscandate" '+%d')
  MM=$(date -d "$blscandate" '+%m')
  YYYY=$(date -d "$blscandate" '+%Y')
  BLPREFIX="${YYYY}-${MM}-${DD}"
  BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
  BLT1SRIMG=$BLSUBJDIR/${BLPREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
  BLASHST1DIR=$BLSUBJDIR/ASHST1
  BLASHST1LSEG=$BLASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
  BLASHST1RSEG=$BLASHST1DIR/final/${id}_right_lfseg_heur.nii.gz
  BLASHSICVDIR=$BLSUBJDIR/ASHSICV

  SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}_ManSegCorrected
  if [[ ! -d $SUBJALOHADIR ]]; then
    SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
  fi

  # output status
  echo "Checking ASHS ICV: $i,$rid,$PREFIX"
  set +e

  # baseline row
  #BLROW=$(cat -A $DEMOG_TXT | grep $id | sed -e 's/\^M//g' | sed -e 's/\$//g')
  #if [[ $BLROW == "" ]]; then
  #  BLROW=$(cat -A $DEMOG_TXT | head -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | awk  -F "," '{for(i=1;i<=NF-1;i++)printf ","}')
  #fi

  # check if exist
  EXIST=0
  if [[ $Type != "T1Longitudinal" ]]; then
    EXIST="-1"
  elif [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $BLASHST1LSEG && -f $BLASHST1RSEG ]]; then
    EXIST=1
  else
    EXIST=0
  fi

  # compute date diff
  date_diff=$(( \
    ($(date -d $scandate +%s) - \
    $(date -d $blscandate +%s) )/(60*60*24) ))
  date_diff=$(echo ${date_diff#-})

  BLICV=$(cat $BLASHSICVDIR/final/${id}_left_corr_nogray_volumes.txt | awk -F' ' '{ print $5 }')

  # get QC information
  QC_TXT=$OUTPUTDIR/QC_ALOHA_ADNIpath20190507.csv
  QCROW=$(cat -A $QC_TXT | grep $id | grep $scandate | head -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
  EXCL_L_Col=$(csvcol.sh $QC_TXT EXCL_L)
  EXCL_R_Col=$(csvcol.sh $QC_TXT EXCL_R)
  EXCL_MRI_Col=$(csvcol.sh $QC_TXT EXCL_MRI)
  excl_l=$(echo $QCROW | cut -f $EXCL_L_Col -d ",")
  excl_r=$(echo $QCROW | cut -f $EXCL_R_Col -d ",")
  excl_mri=$(echo $QCROW | cut -f $EXCL_MRI_Col -d ",")
  if [[ $excl_l == "" ]]; then
    excl_l=0
  fi
  if [[ $excl_r == "" ]]; then
    excl_r=0
  fi
  if [[ $excl_mri == "" ]]; then
    excl_mri=0
  fi
  QCINFO="$excl_l,$excl_r,$excl_mri"


  # go through all the fields 
  RAWMEASURE=""
  LONGIMEASURE=""
  for type in volumes mean_thickness median_thickness; do
  for side in left right; do
    for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do

    if [[ $Type == "T1Longitudinal" ]]; then

      SUBJALOHASUBDIR=$SUBJALOHADIR/results_${sub}
      if [[ $EXIST == 1 ]]; then
        # get volume and thickness
        if [[ -f $SUBJALOHASUBDIR/${type}_${side}.txt ]]; then
          MEA=$(cat $SUBJALOHASUBDIR/${type}_${side}.txt)
          MEA="${MEA// /,}"
          BL=$(echo $MEA | cut -d , -f 1)
          FU=$(echo $MEA | cut -d , -f 2)
          CHANGE=$(echo "scale=7;($BL-$FU)*100*365/$BL/$date_diff" | bc)
          LONGIMEASURE="$LONGIMEASURE,$CHANGE"          

        else
          MEA=","
          EXIST=0
          LONGIMEASURE="$LONGIMEASURE,"
        fi
        RAWMEASURE="$RAWMEASURE,$MEA"
      else
        RAWMEASURE="$RAWMEASURE,,"
        LONGIMEASURE="$LONGIMEASURE,"
      fi

    else

      if [[ $type == "volumes" ]]; then

        SEG=$BLASHST1DIR/final/${id}_${side}_lfseg_heurLW.nii.gz
        if [[ -f $SEG ]]; then

          if [[ $sub == "ERC" ]]; then
            label=(10 10)
          elif [[ $sub == "BA35" ]]; then
            label=(11 11)
          elif [[ $sub == "BA36" ]]; then
            label=(12 12)
          elif [[ $sub == "PHC" ]]; then
            label=(13 13)
          elif [[ $sub == "AHippo" ]]; then
            label=(1 1)
          elif [[ $sub == "PHippo" ]]; then
            label=(2 2)
          elif [[ $sub == "Hippo" ]]; then
            label=(1 2)
          fi

          VOL=$(c3d $SEG -thresh ${label[0]} ${label[1]} 1 0 \
                -dup -lstat \
                | awk '{print $7}' | awk '{printf("%s ", $1)}' \
                | awk '{printf("%2.4f", $3)}')

        elif [[ -f $BLASHST1DIR/final/${id}_${side}_lfseg_heur.nii.gz ]]; then

          SEG=$BLASHST1DIR/final/${id}_${side}_lfseg_heur.nii.gz

          if [[ $sub == "AHippo" ]]; then
            txt="Anterior_hippocampus"
          elif [[ $sub == "PHippo" ]]; then
            txt="Posterior_hippocampus"
          elif [[ $sub == "ERC" ]]; then
            txt="ERC"
          elif [[ $sub == "BA35" ]]; then
            txt="Br35"
          elif [[ $sub == "BA36" ]]; then
            txt="Br36"
          elif [[ $sub == "PHC" ]]; then
            txt=" PHC"
          fi
      
          if [[ $sub == "Hippo" ]]; then
            VOL=$(c3d $BLASHST1DIR/final/${id}_${side}_lfseg_heur.nii.gz \
                -thresh 1 2 1 0 -dup -lstat \
                | awk '{print $7}' | awk '{printf("%s ", $1)}' \
                | awk '{printf("%2.4f", $3)}')
          else
            VOL=$(cat $BLASHST1DIR/final/${id}_${side}_heur_volumes.txt | grep "$txt" | cut -d ' ' -f 5)
          fi
        else
          VOL=""
        fi
        MEA="$VOL,$VOL"
        RAWMEASURE="$RAWMEASURE,$MEA"
        LONGIMEASURE="$LONGIMEASURE,"
    else
      RAWMEASURE="$RAWMEASURE,,"
      LONGIMEASURE="$LONGIMEASURE,"
    fi

  fi

  done
done
done

  # get thickness measurements from the fitted mesh
  if [[ 0 == 1 ]]; then
  SUBJALOHASUBDIR=$SUBJALOHADIR/results_MTLCortexThk
  for type in mean_thickness median_thickness; do
  for side in left right; do
    idx=3
    BLTXT=$SUBJALOHASUBDIR/${id}_${side}_${type}_bl.txt
    if [[ $type == "mean_thickness" ]]; then
      FUTXT=$SUBJALOHASUBDIR/${id}_${side}_${type}_fu.txt
    else
      FUTXT=$SUBJALOHASUBDIR/${id}_${side}_${type}_fu.vtk
    fi
    for sub in ERC BA35 BA36 PHC; do
      if [[ $EXIST == 1 ]]; then
        if [[ -f $BLTXT && -f $FUTXT ]]; then

          BL=$(cat $BLTXT | cut -d , -f $idx)
          FU=$(cat $FUTXT | cut -d , -f $idx)
          MEA="$BL,$FU"
          CHANGE=$(echo "scale=7;($BL-$FU)*100*365/$BL/$date_diff" | bc)
          LONGIMEASURE="$LONGIMEASURE,$CHANGE"

        else
          MEA=","
          #EXIST=0
          LONGIMEASURE="$LONGIMEASURE,"
        fi
        RAWMEASURE="$RAWMEASURE,$MEA"
      else
        RAWMEASURE="$RAWMEASURE,,"
        LONGIMEASURE="$LONGIMEASURE,"
      fi
      idx=$((idx+1))
    done
  done
  done
  fi


  echo "$ROW,$BLICV,$EXIST,$date_diff,$QCINFO$RAWMEASURE$LONGIMEASURE,-1" >> $OUTTMPDIR/$(printf %04d $i)_${id}.csv

  set -e

}

######################################################
function SummarizeRobin()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  ROBIN_TXT=/data/jux/rdeflores/scripts/thickness_ADNI/longit/list/list_longit_ID_Date.txt
  QC_TXT=$OUTPUTDIR/QC_ALOHA20190402.csv
  N=$(cat $ROBIN_TXT | wc -l)
  #N=10

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # Submit job to copy data
  if [[ 1 == 1 ]]; then
  PREFIX1=SMR
  for ((i=1;i<=${N};i++)); do

    echo "$i submitted"
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX1}_${i}" \
         $0 SummarizeRobin_sub $i $OUTTMPDIR
    sleep 0.05

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  fi

  # get all info
  # header
  blheader=$(cat -A $DEMOG_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header=$(cat -A $QC_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  #sed -e 's/\^M\$//g'
  header="$header,ALOHA_Success,DateDiffFromBLScanDate,ALOHAQC_EXCL_L,ALOHAQC_EXCL_R,ALOHAQC_EXCL_MRI"
  longheader=""
  for type in VOL MeanTHK MedianTHK; do
  for side in L R; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  done
  done
  done
  #for type in MeanTHKFITTED MedianTHKFITTED; do
  #for side in L R; do
  #for sub in ERC BA35 BA36 PHC; do
  #  header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
  #  longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  #done
  #done
  #done
  header="$header$longheader,$blheader,END"
  echo $header > $OUTPUTDIR/longitudinal_Robin.csv

  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/longitudinal_Robin.csv

}

function SummarizeRobin_sub()
{
  i=$1
  OUTTMPDIR=$2
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  ROBIN_TXT=/data/jux/rdeflores/scripts/thickness_ADNI/longit/list/list_longit_ID_Date.txt
  QC_TXT=$OUTPUTDIR/QC_ALOHA20190402.csv

  # get row from Robin spreadsheet
  RBROW=$(cat -A $ROBIN_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g') 
  id=$(echo $RBROW | cut -f 1 -d ",")
  scandate=$(echo $RBROW | cut -f 3 -d ",")
  blscandate=$(echo $RBROW | cut -f 2 -d ",")

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $QC_TXT ID)
  Vis2Col=$(csvcol.sh $QC_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $QC_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $QC_TXT BLSANDATE)
  TYPECol=$(csvcol.sh $QC_TXT AnalysisType)
  EXCL_L_Col=$(csvcol.sh $QC_TXT EXCL_L)
  EXCL_R_Col=$(csvcol.sh $QC_TXT EXCL_R)
  EXCL_MRI_Col=$(csvcol.sh $QC_TXT EXCL_MRI)

  # search allinfo for this line
  ROW=$(cat -A $QC_TXT | grep $id | grep $scandate | grep $blscandate | while read line; do
    tmpid=$(echo $line | cut -f $IDCol -d ",")
    tmpscandate=$(echo $line | cut -f $ScanDateCol -d ",")
    tmpblscandate=$(echo $line | cut -f $BLScanDateCol -d ",")
    if [[ $id == $tmpid && $scandate == $tmpscandate && $blscandate == $tmpblscandate ]]; then
      echo $line
      break
    fi
    done)

  if [[ $ROW == "" ]]; then
    #echo "$id,$scandate,$blscandate,-1" >> $OUTTMPDIR/$(printf %04d $i)_${id}.csv
    exit
  fi
  ROW=$(echo $ROW | sed -e 's/\^M//g' | sed -e 's/\$//g')

  # get information
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  vis2=$(echo $ROW | cut -f $Vis2Col -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  Type=$(echo $ROW | cut -f $TYPECol -d ",")
  SUBJALLDIR=$ALLDATADIR/$id
  FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz

  # baseline information
  blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
  blscandate=$(ReFormateDate $blscandate)
  DD=$(date -d "$blscandate" '+%d')
  MM=$(date -d "$blscandate" '+%m')
  YYYY=$(date -d "$blscandate" '+%Y')
  BLPREFIX="${YYYY}-${MM}-${DD}"
  BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
  BLT1SRIMG=$BLSUBJDIR/${BLPREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
  BLASHST1DIR=$BLSUBJDIR/ASHST1
  BLASHST1LSEG=$BLASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
  BLASHST1RSEG=$BLASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

  SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}

  # output status
  echo "Checking ASHS ICV: $i,$rid,$PREFIX"
  set +e

  # baseline row
  BLROW=$(cat -A $DEMOG_TXT | grep $id | sed -e 's/\^M//g' | sed -e 's/\$//g')
  if [[ $BLROW == "" ]]; then
    BLROW=$(cat -A $DEMOG_TXT | head -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | awk  -F "," '{for(i=1;i<=NF-1;i++)printf ","}')
  fi

  # check if exist
  EXIST=0
  if [[ $Type != "T1Longitudinal" ]]; then
    EXIST="-1"
  elif [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $BLASHST1LSEG && -f $BLASHST1RSEG ]]; then
    EXIST=1
  else
    EXIST=0
  fi

  # compute date diff
  date_diff=$(( \
    ($(date -d $scandate +%s) - \
    $(date -d $blscandate +%s) )/(60*60*24) ))
  date_diff=$(echo ${date_diff#-})

  # get QC information
  excl_l=$(echo $ROW | cut -f $EXCL_L_Col -d ",")
  excl_r=$(echo $ROW | cut -f $EXCL_R_Col -d ",")
  excl_mri=$(echo $ROW | cut -f $EXCL_MRI_Col -d ",")
  if [[ $excl_l == "" ]]; then
    excl_l=0
  fi
  if [[ $excl_r == "" ]]; then
    excl_r=0
  fi
  if [[ $excl_mri == "" ]]; then
    excl_mri=0
  fi
  QCINFO="$excl_l,$excl_r,$excl_mri"


  # go through all the fields
  RAWMEASURE=""
  LONGIMEASURE=""
  for type in volumes mean_thickness median_thickness; do
  for side in left right; do
    for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    SUBJALOHASUBDIR=$SUBJALOHADIR/results_${sub}
      if [[ $EXIST == 1 ]]; then
        # get volume and thickness
        if [[ -f $SUBJALOHASUBDIR/${type}_${side}.txt ]]; then
          MEA=$(cat $SUBJALOHASUBDIR/${type}_${side}.txt)
          MEA="${MEA// /,}"
          BL=$(echo $MEA | cut -d , -f 1)
          FU=$(echo $MEA | cut -d , -f 2)
          CHANGE=$(echo "scale=7;($BL-$FU)*100*365/$BL/$date_diff" | bc)
          LONGIMEASURE="$LONGIMEASURE,$CHANGE"

        else
          MEA=","
          EXIST=0
          LONGIMEASURE="$LONGIMEASURE,"
        fi
        RAWMEASURE="$RAWMEASURE,$MEA"
      else
        RAWMEASURE="$RAWMEASURE,,"
        LONGIMEASURE="$LONGIMEASURE,"
      fi
    done
  done
  done

  echo "$ROW,$EXIST,$date_diff,$QCINFO$RAWMEASURE$LONGIMEASURE,$BLROW,-1" >> $OUTTMPDIR/$(printf %04d $i)_${id}.csv

  set -e

}

######################################################
function MTLCortexRegionalThk()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.tsv
  N=$(cat $ALLINFO_TXT | wc -l)
  N_begin=2
  #N=1064
  rm -f $OUTPUTDIR/log.txt


  # Submit job to copy data
  JPREFIX=RTHK${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')

    # columns
    RIDCol=1
    IDCol=$(csvcol_tab.sh $ALLINFO_TXT ID)
    Vis2Col=$(csvcol_tab.sh $ALLINFO_TXT VISCODE2)
    ScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT BLSANDATE)
    TYPECol=$(csvcol_tab.sh $ALLINFO_TXT AnalysisType)

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    vis2=$(echo $ROW | cut -f $Vis2Col -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    if [[ $Type == "T1Longitudinal" ]]; then
      SUBJDIR=$DATADIR/$id
    else
      continue
    fi
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    DD=$(date -d "$blscandate" '+%d')
    MM=$(date -d "$blscandate" '+%m')
    YYYY=$(date -d "$blscandate" '+%Y')
    BLPREFIX="${YYYY}-${MM}-${DD}"

    # perform thickness analysis
    echo "$i,$rid, $PREFIX"
    SUBJALOHADIR=$SUBJDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}

    if [[ -d $SUBJALOHADIR ]]; then

      qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${i}_${rid}_${PREFIX}" \
           $0 MTLCortexRegionalThk_sub $i \

    fi
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}

function MTLCortexRegionalThk_sub()
{
  i=$1
  ALLINFO_TXT=$OUTPUTDIR/checkdata.tsv
  ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')

  # columns
  RIDCol=1
  IDCol=$(csvcol_tab.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol_tab.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT BLSANDATE)
  TYPECol=$(csvcol_tab.sh $ALLINFO_TXT AnalysisType)

  # get information
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  vis2=$(echo $ROW | cut -f $Vis2Col -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  Type=$(echo $ROW | cut -f $TYPECol -d ",")
  SUBJDIR=$DATADIR/$id

  # baseline info
  blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
  blscandate=$(ReFormateDate $blscandate)
  DD=$(date -d "$blscandate" '+%d')
  MM=$(date -d "$blscandate" '+%m')
  YYYY=$(date -d "$blscandate" '+%Y')
  BLPREFIX="${YYYY}-${MM}-${DD}"
  SUBJALLDIR=$ALLDATADIR/$id
  BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
  BLT1SRIMG=$BLSUBJDIR/${BLPREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
  BLASHST1DIR=$BLSUBJDIR/ASHST1
  BLMSTTHKDIR=$BLSUBJDIR/ASHST1_MTLCORTEX_MSTTHK

  # outdir
  SUBJALOHADIR=$SUBJDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
  OUTDIR=$SUBJALOHADIR/results_MTLCortexThk
  mkdir -p $OUTDIR

  # warp baseline mesh to followup image space
  for side in left right; do

    # baseline info
    BLMESHORIG=$(ls $BLMSTTHKDIR/GeoShoot/$side/template_?_to_${id}_${side}_GSShoot_MRG_thickmap_withlabel.vtk)
    group=$(basename $BLMESHORIG | cut -d _ -f 2)
    BLMESH=$OUTDIR/template_${group}_to_${id}_${side}_GSShoot_MRG_thickmap_withlabel_bl.vtk
    ln -sf $BLMESHORIG $BLMESH
    ln -sf $BLMSTTHKDIR/GeoShoot/$side/${id}_${side}_mean_thickness.vtk \
      $OUTDIR/${id}_${side}_mean_thickness_bl.txt
    ln -sf $BLMSTTHKDIR/GeoShoot/$side/${id}_${side}_median_thickness.vtk \
      $OUTDIR/${id}_${side}_median_thickness_bl.txt
    
    # follow up
    FUMESH=$OUTDIR/template_${group}_to_${id}_${side}_GSShoot_MRG_thickmap_withlabel_fu.vtk
    FUSKEL=$OUTDIR/template_${group}_to_${id}_${side}_GSShoot_MRG_thickmap_withlabel_fu_skel.vtk
    WDDEF=$SUBJALOHADIR/deformable
    warpmesh $BLMESH \
      $TMPDIR/blmptrim_seg_${side}_tohw.vtk \
      $WDDEF/mprage_global_long_${side}_omRAS_half.mat
    warpmesh $TMPDIR/blmptrim_seg_${side}_tohw.vtk \
      $TMPDIR/blmptrim_seg_${side}_warped.vtk \
      $WDDEF/mp_antsreg3d_${side}Warp?vec.nii.gz
    warpmesh $TMPDIR/blmptrim_seg_${side}_warped.vtk \
      $FUMESH \
      $WDDEF/mprage_global_long_${side}_omRAS_half.mat
    cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
      -T $FUMESH \
      -p 1.2 -e 6 \
      $FUMESH \
      $FUSKEL

    # get summary thickness
    $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('/home/longxie/ADNI2018/matlabcode');
    compute_meshlabel_meanthickness('$FUMESH','$id','$side','$OUTDIR/${id}_${side}_mean_thickness_fu.txt','$OUTDIR/${id}_${side}_median_thickness_fu.vtk');
MATCODE


  done
}

######################################################
function ComputeRegionalAtropyRate()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/longitudinal.csv
  #N=$(cat $ALLINFO_TXT | wc -l)
  #N_begin=2
  #N=1064
  rm -f $OUTPUTDIR/log.txt

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)
  SUCCESS=$(csvcol.sh $ALLINFO_TXT ALOHA_Success)
  DateDiffCol=$(csvcol.sh $ALLINFO_TXT DateDiffFromBLScanDate)

  # IDs
  IDs=($(cat $ALLINFO_TXT | awk -F, '{print $2}' | uniq))
  #echo ${IDs[*]}
  N=${#IDs[*]}
  #N=290
  N_begin=0

  # Submit job to copy data
  JPREFIX=RTHK${expid}
  for ((i=${N_begin};i<${N};i++)); do

    id=${IDs[$i]}

    # get cases
    set +e
    fcases=($(cat $ALLINFO_TXT | grep $id | grep T1Longitudinal | sed -e 's/\ //g'))
    set -e
    # if no case is found, skip
    if [[ ${#fcases[*]} == 0 ]]; then
      continue
    fi

    # check how many qualified cases
    cases=""
    idx=0
    lthresh=$(printf "%.0f" $(echo "scale=2;${ylthresh}*365" | bc))
    hthresh=$(printf "%.0f" $(echo "scale=2;${yhthresh}*365" | bc))
    for ((j=0;j<${#fcases[*]};j++)); do
      fcase=${fcases[j]}
      date_diff=$(echo $fcase | cut -f $DateDiffCol -d ",")
      success=$(echo $fcase | cut -f $SUCCESS -d ",")
      if [[ $date_diff -le $hthresh && $date_diff -ge $lthresh && $success -eq 1 ]]; then
        cases[${idx}]=$fcase
        idx=$((idx+1))
      fi
    done

    # if no qualified cases, skip
    if [[ $cases == "" ]]; then
      continue
    fi

    # submit job
    for side in left right; do
      qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${id}_${side}" \
           $0 ComputeRegionalAtropyRate_sub $id $side
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}

function ComputeRegionalAtropyRate_sub()
{
  id=$1
  side=$2
  ALLINFO_TXT=$OUTPUTDIR/longitudinal.csv
  
  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)
  SUCCESS=$(csvcol.sh $ALLINFO_TXT ALOHA_Success)
  DateDiffCol=$(csvcol.sh $ALLINFO_TXT DateDiffFromBLScanDate)

  # get cases
  fcases=($(cat $ALLINFO_TXT | grep $id | grep T1Longitudinal | sed -e 's/\ //g'))
  # if no case is found, skip
  if [[ ${#fcases[*]} == 0 ]]; then
    exit
  fi

  # check how many qualified cases
  cases=""
  idx=0
  lthresh=$(printf "%.0f" $(echo "scale=2;${ylthresh}*365" | bc))
  hthresh=$(printf "%.0f" $(echo "scale=2;${yhthresh}*365" | bc))
  for ((j=0;j<${#fcases[*]};j++)); do
    fcase=${fcases[j]}
    date_diff=$(echo $fcase | cut -f $DateDiffCol -d ",")
    success=$(echo $fcase | cut -f $SUCCESS -d ",")
    if [[ $date_diff -le $hthresh && $date_diff -ge $lthresh && $success -eq 1 ]]; then
      cases[${idx}]=$fcase
      idx=$((idx+1))
    fi
  done

  # if no qualified cases, skip
  if [[ $cases == "" ]]; then
    exit
  fi

  # get baseline information
  # get information
  case=${cases[0]}
  rid=$(echo $case | cut -f $RIDCol -d ",")
  blscandate=$(echo $case | cut -f $BLScanDateCol -d ",")
  blscandate=$(ReFormateDate $blscandate)
  DD=$(date -d "$blscandate" '+%d')
  MM=$(date -d "$blscandate" '+%m')
  YYYY=$(date -d "$blscandate" '+%Y')
  BLPREFIX="${YYYY}-${MM}-${DD}"  
  SUBJALLDIR=$ALLDATADIR/$id
  SUBJDIR=$DATADIR/$id
  BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
  BLASHST1DIR=$BLSUBJDIR/ASHST1
  BLMSTTHKDIR=$BLSUBJDIR/ASHST1_MTLCORTEX_MSTTHK
  REGIONALTHKDIR=$SUBJDIR/AlohaMTL_ASHST1_${BLPREFIX}_MTLCortexRegionalThk
  mkdir -p $REGIONALTHKDIR

  # baseline mesh
  BLMESHORIG=$(ls $BLMSTTHKDIR/GeoShoot/$side/template_?_to_${id}_${side}_GSShoot_MRG_thickmap_withlabel.vtk)
  group=$(basename $BLMESHORIG | cut -d _ -f 2)
  MESHES="$BLMESHORIG"

  # get all the meshes
  INFO=${#cases[*]}
  for ((i=0;i<${#cases[*]};i++)); do

    case=${cases[i]}
    date_diff=$(echo $case | cut -f $DateDiffCol -d ",")
    scandate=$(echo $case | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    SUBJALOHADIR=$SUBJDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
    FUMESH=$SUBJALOHADIR/results_MTLCortexThk/template_${group}_to_${id}_${side}_GSShoot_MRG_thickmap_withlabel_fu.vtk
    MESHES="$MESHES $FUMESH"
    INFO="$INFO,$date_diff"

  done   
 
  # support information
  echo $INFO > $REGIONALTHKDIR/info_${side}.csv
 
  # merge the meshes
  mesh_merge_arrays \
    -r $BLMESHORIG \
    $REGIONALTHKDIR/template_${group}_to_${id}_${side}_GSShoot_MRG_thickmap_withlabel_alltp_${ylthresh}-${yhthresh}yrs.vtk \
    Thickness $MESHES

  # compute annualized change
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('/home/longxie/ADNI2018/matlabcode');
    compute_regional_annualized_atrophy_rate('$REGIONALTHKDIR/template_${group}_to_${id}_${side}_GSShoot_MRG_thickmap_withlabel_alltp_${ylthresh}-${yhthresh}yrs.vtk','$REGIONALTHKDIR/info_${side}.csv');
MATCODE
}

######################################################
function RegionalStats_meshglm()
{

  if [[ 1 == 1 ]]; then
  rm -rf $DESIGNDIR
  mkdir -p $DESIGNDIR
  ALLINFO_TXT=$OUTPUTDIR/longitudinal_multiteimpoint_2,5year.csv

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSANDATE)
  AgeCol=$(csvcol.sh $ALLINFO_TXT AgeatMRIbl)
  DXCol=$(csvcol.sh $ALLINFO_TXT DX2AMYDATEDIFFQUALFIT)
  QCMRICol=$(csvcol.sh $ALLINFO_TXT QC_MRI)
  MaxDateDiffCol=$(csvcol.sh $ALLINFO_TXT MaxDateDiff)
  IncludeCol=$(csvcol.sh $ALLINFO_TXT INCLUDE)
  PTAUBLCol=$(csvcol.sh $ALLINFO_TXT PTAU_bl)

  # IDs
  IDs=($(cat $ALLINFO_TXT | awk -F, '{print $2}' | uniq))
  N=${#IDs[*]}  
  #echo ${IDs[*]}
  # generate design matrix
  for ((i=0;i<${N};i++)); do

    id=${IDs[$i]}
    echo $id
    # get cases
    set +e
    BLINFO=$(cat $ALLINFO_TXT | grep $id)
    set -e

    # if no case is found, skip
    if [[ $BLINFO == "" ]]; then
      continue
    fi

    # get baseline information
    rid=$(echo $BLINFO | cut -f $RIDCol -d ",")
    blscandate=$(echo $BLINFO | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    DD=$(date -d "$blscandate" '+%d')
    MM=$(date -d "$blscandate" '+%m')
    YYYY=$(date -d "$blscandate" '+%Y')
    BLPREFIX="${YYYY}-${MM}-${DD}"
    SUBJALLDIR=$ALLDATADIR/$id
    SUBJDIR=$DATADIR/$id
    BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
    BLASHST1DIR=$BLSUBJDIR/ASHST1
    BLMSTTHKDIR=$BLSUBJDIR/ASHST1_MTLCORTEX_MSTTHK
    REGIONALTHKDIR=$SUBJDIR/AlohaMTL_ASHST1_${BLPREFIX}_MTLCortexRegionalThk

    for side in left right; do    


      # baseline mesh
      set +e
      MESH=$(ls $REGIONALTHKDIR/template_?_to_${id}_${side}_GSShoot_MRG_thickmap_withlabel_alltp_${ylthresh}-${yhthresh}yrs.vtk)
      set -e

      # if mesh does not exist, skip
      if [[ $MESH == "" ]]; then 
        continue
      fi

      # get related information
      grp=$(basename $MESH | cut -d _ -f 2)
      age=$(echo $BLINFO | cut -f $AgeCol -d ",")
      dx=$(echo $BLINFO | cut -f $DXCol -d ",")
      qcmri=$(echo $BLINFO | cut -f $QCMRICol -d ",")
      maxdatediff=$(echo $BLINFO | cut -f $MaxDateDiffCol -d ",")
      include=$(echo $BLINFO | cut -f $IncludeCol -d ",")
      ptaubl=$(echo $BLINFO | cut -f $PTAUBLCol -d ",")
         
      if [[ $qcmri == 0 && $(echo "$maxdatediff > $minfwtime" | bc) == 1  && $include == 1 ]]; then
        echo "$id,$side,$maxdatediff"
        if [[ $dx == 1 ]]; then
          str="$id 1 0 $age"
          echo $str >> $DESIGNDIR/design_${side}_${grp}_PCLINICALADvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_PCLINICALADvsNC_group-agecr.txt
          echo $str >> $DESIGNDIR/design_${side}_${grp}_ABPOSEMCIvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_ABPOSEMCIvsNC_group-agecr.txt
          echo $str >> $DESIGNDIR/design_${side}_${grp}_ABPOSLMCIvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_ABPOSLMCIvsNC_group-agecr.txt
          echo $str >> $DESIGNDIR/design_${side}_${grp}_ABPOSADvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_ABPOSADvsNC_group-agecr.txt
 
          # use ptau to diconamize
          echo $str >> $DESIGNDIR/design_${side}_${grp}_PTAU24POSPCLINICALADvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_PTAU24POSPCLINICALADvsNC_group-agecr.txt
          echo $str >> $DESIGNDIR/design_${side}_${grp}_PTAUMEDIANPOSPCLINICALADvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_PTAUMEDIANPOSPCLINICALADvsNC_group-agecr.txt


        elif [[ $dx == 2 ]]; then
          str="$id 0 1 $age"
          echo $str >> $DESIGNDIR/design_${side}_${grp}_PCLINICALADvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_PCLINICALADvsNC_group-agecr.txt

          # use ptau to diconamize
          if [[ $(echo "$ptaubl >= 24" | bc) == 1 ]]; then
          echo $str >> $DESIGNDIR/design_${side}_${grp}_PTAU24POSPCLINICALADvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_PTAU24POSPCLINICALADvsNC_group-agecr.txt
          fi

          if [[ $(echo "$ptaubl >= 19.97" | bc) == 1 ]]; then
          echo $str >> $DESIGNDIR/design_${side}_${grp}_PTAUMEDIANPOSPCLINICALADvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_PTAUMEDIANPOSPCLINICALADvsNC_group-agecr.txt
          fi


        elif [[ $dx == 4 ]]; then
          str="$id 0 1 $age"
          echo $str >> $DESIGNDIR/design_${side}_${grp}_ABPOSEMCIvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_ABPOSEMCIvsNC_group-agecr.txt
        elif [[ $dx == 6 ]]; then
          str="$id 0 1 $age"
          echo $str >> $DESIGNDIR/design_${side}_${grp}_ABPOSLMCIvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_ABPOSLMCIvsNC_group-agecr.txt
        elif [[ $dx == 8 ]]; then
          str="$id 0 1 $age"
          echo $str >> $DESIGNDIR/design_${side}_${grp}_ABPOSADvsNC_group-agecr.txt
          echo $MESH >> $DESIGNDIR/meshes_${side}_${grp}_ABPOSADvsNC_group-agecr.txt
        fi
      fi
    done
  done

  # echo contrast
  echo "1 -1 0" > \
    $DESIGNDIR/contrast_PCLINICALADvsNC_group-agecr_nc-pat.txt
  echo "-1 1 0" > \
    $DESIGNDIR/contrast_PCLINICALADvsNC_group-agecr_pat-nc.txt
  echo "1 -1 0" > \
    $DESIGNDIR/contrast_PTAU24POSPCLINICALADvsNC_group-agecr_nc-pat.txt
  echo "-1 1 0" > \
    $DESIGNDIR/contrast_PTAU24POSPCLINICALADvsNC_group-agecr_pat-nc.txt
  echo "1 -1 0" > \
    $DESIGNDIR/contrast_PTAUMEDIANPOSPCLINICALADvsNC_group-agecr_nc-pat.txt
  echo "-1 1 0" > \
    $DESIGNDIR/contrast_PTAUMEDIANPOSPCLINICALADvsNC_group-agecr_pat-nc.txt
  echo "1 -1 0" > \
    $DESIGNDIR/contrast_ABPOSEMCIvsNC_group-agecr_nc-pat.txt
  echo "-1 1 0" > \
    $DESIGNDIR/contrast_ABPOSEMCIvsNC_group-agecr_pat-nc.txt
  echo "1 -1 0" > \
    $DESIGNDIR/contrast_ABPOSLMCIvsNC_group-agecr_nc-pat.txt
  echo "-1 1 0" > \
    $DESIGNDIR/contrast_ABPOSLMCIvsNC_group-agecr_pat-nc.txt
  echo "1 -1 0" > \
    $DESIGNDIR/contrast_ABPOSADvsNC_group-agecr_nc-pat.txt
  echo "-1 1 0" > \
    $DESIGNDIR/contrast_ABPOSADvsNC_group-agecr_pat-nc.txt

  fi

  # submit jobs to run statistical analysis
  PREFIX=RTS${expid}
  if [[ 1  == 1 ]]; then
  for grp in 1 2; do
    for side in left right; do
      rm -rf $STATSDIR/template${grp}_${side}
      for design in $(ls $DESIGNDIR | grep design | grep $side | grep $grp); do

        exp=$(echo $design | sed -e "s/^design_${side}_${grp}_//" | sed -e "s/\.txt//")

        qsub -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -l h_vmem=4.1G,s_vmem=4G \
          -N "${PREFIX}_${exp}_${side}_${grp}_Atrophy" \
          $0 RegionalStats_meshglm_sub $exp $side $grp
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

function RegionalStats_meshglm_sub()
{
  exp=$1
  side=$2
  grp=$3

  # Create the work directory for this analysis
  WORK=$STATSDIR/template${grp}_${side}/design_${exp}_MRG
  mkdir -p $WORK

  # Get the list of subjects
  DESIGNTXT=$DESIGNDIR/design_${side}_${grp}_${exp}.txt
  SUBJ=$(cat $DESIGNTXT | awk '{print $1}')

  # Generate the design matrix for meshglm
  cat $DESIGNTXT | awk '{$1=""; print}' > $WORK/design_${side}.txt

  # combine thickness measurements
  MESHES=$(cat $DESIGNDIR/meshes_${side}_${grp}_${exp}.txt)
  #MESHES=$(for id in $SUBJ; do \
  #  echo " $(ls $DATADIR/$id/*/ASHST1_MTLCORTEX_MSTTHK/GeoShoot/${side}/template_${grp}_to_${id}_${side}_GSShoot_MRG_thickmap.vtk) "; done)

  mesh_merge_arrays -r \
    $GSTEMPDIR/GSTemplate/gshoot/template_${grp}/template/iter_2/template_${grp}_gshoot_MRG.vtk \
    $WORK/thick_${side}_${grp}_MRG.vtk ChangeAnnualizedPercentageAllTP $MESHES

  # Go through the list of contrasts
  for con in $(ls $DESIGNDIR | grep "contrast_${exp}_.*\.txt"); do

    # Get the suffix for this contrast
    suffix=$(echo $con | sed -e "s/^contrast_${exp}_//" -e "s/\.txt//")

    # Create the directory for this contrast
    CWORK=$WORK/contrast_${suffix}
    mkdir -p $CWORK

    FULLNM="design_${exp}_con_${suffix}"

    # Copy the contrast
    cp $DESIGNDIR/$con $CWORK/contrast.txt

    # Build the list of meshes to include
    MESHPARAM="-m $WORK/thick_${side}_${grp}_MRG.vtk $CWORK/thickstat_${FULLNM}_${side}_${grp}_MRG.vtk"
    meshglm $MESHPARAM \
      -g $WORK/design_${side}.txt $CWORK/contrast.txt \
      -a ChangeAnnualizedPercentageAllTP \
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
function MTLCortexRegionalThkMultiTemp()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.tsv
  N=$(cat $ALLINFO_TXT | wc -l)
  N_begin=2
  #N=1064
  rm -f $OUTPUTDIR/log.txt


  # Submit job to copy data
  JPREFIX=RTHK${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')

    # columns
    RIDCol=1
    IDCol=$(csvcol_tab.sh $ALLINFO_TXT ID)
    Vis2Col=$(csvcol_tab.sh $ALLINFO_TXT VISCODE2)
    ScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT BLSANDATE)
    TYPECol=$(csvcol_tab.sh $ALLINFO_TXT AnalysisType)

    # get information
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    vis2=$(echo $ROW | cut -f $Vis2Col -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    if [[ $Type == "T1Longitudinal" ]]; then
      SUBJDIR=$DATADIR/$id
    else
      continue
    fi
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    DD=$(date -d "$blscandate" '+%d')
    MM=$(date -d "$blscandate" '+%m')
    YYYY=$(date -d "$blscandate" '+%Y')
    BLPREFIX="${YYYY}-${MM}-${DD}"
  
    # perform thickness analysis
    echo "$i,$rid, $PREFIX"
    SUBJALOHADIR=$SUBJDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}

    if [[ -d $SUBJALOHADIR ]]; then

      qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${i}_${rid}_${PREFIX}" \
           $0 MTLCortexRegionalThkMultiTemp_sub $i \

    fi
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}

function MTLCortexRegionalThkMultiTemp_sub()
{
  i=$1
  ALLINFO_TXT=$OUTPUTDIR/checkdata.tsv
  ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')

  # columns
  RIDCol=1
  IDCol=$(csvcol_tab.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol_tab.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT BLSANDATE)
  TYPECol=$(csvcol_tab.sh $ALLINFO_TXT AnalysisType)

  # get information
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  vis2=$(echo $ROW | cut -f $Vis2Col -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  Type=$(echo $ROW | cut -f $TYPECol -d ",")
  SUBJDIR=$DATADIR/$id

  # baseline info
  blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
  blscandate=$(ReFormateDate $blscandate)
  DD=$(date -d "$blscandate" '+%d')
  MM=$(date -d "$blscandate" '+%m')
  YYYY=$(date -d "$blscandate" '+%Y')
  BLPREFIX="${YYYY}-${MM}-${DD}"
  SUBJALLDIR=$ALLDATADIR/$id
  BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
  BLT1SRIMG=$BLSUBJDIR/${BLPREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
  BLASHST1DIR=$BLSUBJDIR/ASHST1
  BLMTTHKDIR=$BLSUBJDIR/ASHST1_MTLCORTEX_MultiTempThk

  # outdir
  SUBJALOHADIR=$SUBJDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
  OUTDIR=$SUBJALOHADIR/results_MTLCortexThkMultiTemp
  mkdir -p $OUTDIR

  # warp baseline mesh to followup image space
  for type in Mult Uni; do
  for side in left right; do

    # baseline info
    BLMESHORIG=$(ls $BLMTTHKDIR/RegTo${type}Temp/${id}_${side}*MRG_thickmap_withlabel.vtk)
    BLMESH=$OUTDIR/${id}_${side}_MRG_thickmap_withlabel_${type}Temp_bl.vtk
    ln -sf $BLMESHORIG $BLMESH
    ln -sf $BLMTTHKDIR/RegTo${type}Temp/${id}_${side}_mean_thickness.txt $OUTDIR/${id}_${side}_mean_thickness_${type}Temp_bl.txt
    ln -sf $BLMTTHKDIR/RegTo${type}Temp/${id}_${side}_median_thickness.txt $OUTDIR/${id}_${side}_median_thickness_${type}Temp_bl.txt

    # follow up
    FUMESH=$OUTDIR/${id}_${side}_MRG_thickmap_withlabel_${type}Temp_fu.vtk
    FUSKEL=$OUTDIR/${id}_${side}_MRG_thickmap_withlabel_${type}Temp_fu_skel.vtk
    WDDEF=$SUBJALOHADIR/deformable
    warpmesh $BLMESH \
      $TMPDIR/blmptrim_seg_${side}_tohw.vtk \
      $WDDEF/mprage_global_long_${side}_omRAS_half.mat
    warpmesh $TMPDIR/blmptrim_seg_${side}_tohw.vtk \
      $TMPDIR/blmptrim_seg_${side}_warped.vtk \
      $WDDEF/mp_antsreg3d_${side}Warp?vec.nii.gz
    warpmesh $TMPDIR/blmptrim_seg_${side}_warped.vtk \
      $FUMESH \
      $WDDEF/mprage_global_long_${side}_omRAS_half.mat
    cmrep_vskel -Q /data/picsl/pauly/bin/qvoronoi \
      -T $FUMESH \
      -p 1.2 -e 6 \
      $FUMESH \
      $FUSKEL

    # get summary thickness
    $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('/home/longxie/ADNI2018/matlabcode');
    compute_meshlabel_meanthickness('$FUMESH','$id','$side','$OUTDIR/${id}_${side}_mean_thickness_${type}Temp_fu.txt','$OUTDIR/${id}_${side}_median_thickness_${type}Temp_fu.txt');
MATCODE

  done
  done
}
######################################################
function ComputeRegionalAtropyRateMultiTemp()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/longitudinal.csv
  #N=$(cat $ALLINFO_TXT | wc -l)
  #N_begin=2
  #N=1064
  rm -f $OUTPUTDIR/log.txt

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)
  SUCCESS=$(csvcol.sh $ALLINFO_TXT ALOHA_Success)
  DateDiffCol=$(csvcol.sh $ALLINFO_TXT DateDiffFromBLScanDate)

  # IDs
  IDs=($(cat $ALLINFO_TXT | awk -F, '{print $2}' | uniq))
  #echo ${IDs[*]}
  N=${#IDs[*]}
  #N=290
  N_begin=0

  # Submit job to copy data
  JPREFIX=RETHK${expid}
  for ((i=${N_begin};i<${N};i++)); do

    id=${IDs[$i]}

    # get cases
    set +e
    fcases=($(cat $ALLINFO_TXT | grep $id | grep T1Longitudinal | sed -e 's/\ //g'))
    set -e
    # if no case is found, skip
    if [[ ${#fcases[*]} == 0 ]]; then
      continue
    fi

    # check how many qualified cases
    cases=""
    idx=0
    lthresh=$(printf "%.0f" $(echo "scale=2;${ylthresh}*365" | bc))
    hthresh=$(printf "%.0f" $(echo "scale=2;${yhthresh}*365" | bc))
    for ((j=0;j<${#fcases[*]};j++)); do
      fcase=${fcases[j]}
      date_diff=$(echo $fcase | cut -f $DateDiffCol -d ",")
      success=$(echo $fcase | cut -f $SUCCESS -d ",")
      if [[ $date_diff -le $hthresh && $date_diff -ge $lthresh && $success -eq 1 ]]; then
        cases[${idx}]=$fcase
        idx=$((idx+1))
      fi
    done

    # if no qualified cases, skip
    if [[ $cases == "" ]]; then
      continue
    fi

    # submit job
    for side in left right; do
      qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${id}_${side}" \
           $0 ComputeRegionalAtropyRateMultiTemp_sub $id $side
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}

function ComputeRegionalAtropyRateMultiTemp_sub()
{
  id=$1
  side=$2
  ALLINFO_TXT=$OUTPUTDIR/longitudinal.csv

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)
  SUCCESS=$(csvcol.sh $ALLINFO_TXT ALOHA_Success)
  DateDiffCol=$(csvcol.sh $ALLINFO_TXT DateDiffFromBLScanDate)

  # get cases
  fcases=($(cat $ALLINFO_TXT | grep $id | grep T1Longitudinal | sed -e 's/\ //g'))
  # if no case is found, skip
  if [[ ${#fcases[*]} == 0 ]]; then
    exit
  fi

  # check how many qualified cases
  cases=""
  idx=0
  lthresh=$(printf "%.0f" $(echo "scale=2;${ylthresh}*365" | bc))
  hthresh=$(printf "%.0f" $(echo "scale=2;${yhthresh}*365" | bc))
  for ((j=0;j<${#fcases[*]};j++)); do
    fcase=${fcases[j]}
    date_diff=$(echo $fcase | cut -f $DateDiffCol -d ",")
    success=$(echo $fcase | cut -f $SUCCESS -d ",")
    if [[ $date_diff -le $hthresh && $date_diff -ge $lthresh && $success -eq 1 ]]; then
      cases[${idx}]=$fcase
      idx=$((idx+1))
    fi
  done

  # if no qualified cases, skip
  if [[ $cases == "" ]]; then
    exit
  fi

  # get baseline information
  # get information
  case=${cases[0]}
  rid=$(echo $case | cut -f $RIDCol -d ",")
  blscandate=$(echo $case | cut -f $BLScanDateCol -d ",")
  blscandate=$(ReFormateDate $blscandate)
  DD=$(date -d "$blscandate" '+%d')
  MM=$(date -d "$blscandate" '+%m')
  YYYY=$(date -d "$blscandate" '+%Y')
  BLPREFIX="${YYYY}-${MM}-${DD}"
  SUBJALLDIR=$ALLDATADIR/$id
  SUBJDIR=$DATADIR/$id
  BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
  BLASHST1DIR=$BLSUBJDIR/ASHST1
  BLMTTHKDIR=$BLSUBJDIR/ASHST1_MTLCORTEX_MultiTempThk
  REGIONALTHKDIR=$SUBJDIR/AlohaMTL_ASHST1_${BLPREFIX}_MTLCortexRegionalThk_MultiTemp
  mkdir -p $REGIONALTHKDIR


  # baseline mesh
  type="Uni"
  BLMESHORIG=$(ls $BLMTTHKDIR/RegTo${type}Temp/${id}_${side}*MRG_thickmap_withlabel.vtk)
  MESHES="$BLMESHORIG"

  # get all the meshes
  INFO=${#cases[*]}
  for ((i=0;i<${#cases[*]};i++)); do

    case=${cases[i]}
    date_diff=$(echo $case | cut -f $DateDiffCol -d ",")
    scandate=$(echo $case | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    SUBJALOHADIR=$SUBJDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
    FUMESH=$SUBJALOHADIR/results_MTLCortexThkMultiTemp/${id}_${side}_MRG_thickmap_withlabel_${type}Temp_fu.vtk
    MESHES="$MESHES $FUMESH"
    INFO="$INFO,$date_diff"

  done

  # support information
  echo $INFO > $REGIONALTHKDIR/info_${side}.csv

  # merge the meshes
  mesh_merge_arrays \
    -r $BLMESHORIG \
    $REGIONALTHKDIR/${id}_${side}_MRG_thickmap_withlabel_${type}Temp_alltp_${ylthresh}-${yhthresh}yrs.vtk \
    Thickness $MESHES

  # compute annualized change
  $MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
    addpath('/home/longxie/ADNI2018/matlabcode');
    compute_regional_annualized_atrophy_rate('$REGIONALTHKDIR/${id}_${side}_MRG_thickmap_withlabel_${type}Temp_alltp_${ylthresh}-${yhthresh}yrs.vtk','$REGIONALTHKDIR/info_${side}.csv');
MATCODE
}

######################################################
function UTRegionalStats_meshglm()
{

  if [[ 1 == 1 ]]; then
  rm -rf $MTDESIGNDIR
  mkdir -p $MTDESIGNDIR
  ALLINFO_TXT=$OUTPUTDIR/longitudinal_multiteimpoint_2,5year.csv

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSANDATE)
  AgeCol=$(csvcol.sh $ALLINFO_TXT AgeatMRIbl)
  DXCol=$(csvcol.sh $ALLINFO_TXT DX2AMYDATEDIFFQUALFIT)
  QCMRICol=$(csvcol.sh $ALLINFO_TXT QC_MRI)
  MaxDateDiffCol=$(csvcol.sh $ALLINFO_TXT MaxDateDiff)
  IncludeCol=$(csvcol.sh $ALLINFO_TXT INCLUDE)
  PTAUBLCol=$(csvcol.sh $ALLINFO_TXT PTAU_bl)

  # IDs
  IDs=($(cat $ALLINFO_TXT | awk -F, '{print $2}' | uniq))
  N=${#IDs[*]}
  #echo ${IDs[*]}
  # generate design matrix
  for ((i=0;i<${N};i++)); do

    id=${IDs[$i]}
    echo $id
    # get cases
    set +e
    BLINFO=$(cat $ALLINFO_TXT | grep $id)
    set -e

    # if no case is found, skip
    if [[ $BLINFO == "" ]]; then
      continue
    fi

    # get baseline information
    rid=$(echo $BLINFO | cut -f $RIDCol -d ",")
    blscandate=$(echo $BLINFO | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    DD=$(date -d "$blscandate" '+%d')
    MM=$(date -d "$blscandate" '+%m')
    YYYY=$(date -d "$blscandate" '+%Y')
    BLPREFIX="${YYYY}-${MM}-${DD}"
    SUBJALLDIR=$ALLDATADIR/$id
    SUBJDIR=$DATADIR/$id
    BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
    SUBJALLDIR=$ALLDATADIR/$id
    SUBJDIR=$DATADIR/$id
    BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
    BLASHST1DIR=$BLSUBJDIR/ASHST1
    BLMSTTHKDIR=$BLSUBJDIR/ASHST1_MTLCORTEX_MultiTempThk
    REGIONALTHKDIR=$SUBJDIR/AlohaMTL_ASHST1_${BLPREFIX}_MTLCortexRegionalThk_MultiTemp

    for side in left right; do


      # baseline mesh
      set +e
      type="Uni"
      MESH=$(ls $REGIONALTHKDIR/${id}_${side}_MRG_thickmap_withlabel_${type}Temp_alltp_${ylthresh}-${yhthresh}yrs.vtk)
      set -e

      # if mesh does not exist, skip
      if [[ $MESH == "" ]]; then
        continue
      fi

      # get related information
      #grp=$(basename $MESH | cut -d _ -f 2)
      age=$(echo $BLINFO | cut -f $AgeCol -d ",")
      dx=$(echo $BLINFO | cut -f $DXCol -d ",")
      qcmri=$(echo $BLINFO | cut -f $QCMRICol -d ",")
      maxdatediff=$(echo $BLINFO | cut -f $MaxDateDiffCol -d ",")
      include=$(echo $BLINFO | cut -f $IncludeCol -d ",")
      ptaubl=$(echo $BLINFO | cut -f $PTAUBLCol -d ",")

      if [[ $qcmri == 0 && $(echo "$maxdatediff > $minfwtime" | bc) == 1  && $include == 1 ]]; then
        echo "$id,$side,$maxdatediff"
        if [[ $dx == 1 ]]; then
          str="$id 1 0 $age"
          echo $str >> $MTDESIGNDIR/design_${side}_PCLINICALADvsNC_group-agecr.txt
          echo $MESH >> $MTDESIGNDIR/meshes_${side}_PCLINICALADvsNC_group-agecr.txt
          echo $str >> $MTDESIGNDIR/design_${side}_ABPOSEMCIvsNC_group-agecr.txt
          echo $MESH >> $MTDESIGNDIR/meshes_${side}_ABPOSEMCIvsNC_group-agecr.txt
          echo $MESH >> $MTDESIGNDIR/meshes_${side}_ABPOSEMCIvsNC_group-agecr.txt
          echo $str >> $MTDESIGNDIR/design_${side}_ABPOSLMCIvsNC_group-agecr.txt
          echo $MESH >> $MTDESIGNDIR/meshes_${side}_ABPOSLMCIvsNC_group-agecr.txt
          echo $str >> $MTDESIGNDIR/design_${side}_ABPOSADvsNC_group-agecr.txt
          echo $MESH >> $MTDESIGNDIR/meshes_${side}_ABPOSADvsNC_group-agecr.txt

          # use ptau to diconamize
          echo $str >> $MTDESIGNDIR/design_${side}_PTAU24POSPCLINICALADvsNC_group-agecr.txt
          echo $MESH >> $MTDESIGNDIR/meshes_${side}_PTAU24POSPCLINICALADvsNC_group-agecr.txt


        elif [[ $dx == 2 ]]; then
          str="$id 0 1 $age"
          echo $str >> $MTDESIGNDIR/design_${side}_PCLINICALADvsNC_group-agecr.txt
          echo $MESH >> $MTDESIGNDIR/meshes_${side}_PCLINICALADvsNC_group-agecr.txt

          # use ptau to diconamize
          if [[ $(echo "$ptaubl >= 24" | bc) == 1 ]]; then
          echo $str >> $MTDESIGNDIR/design_${side}_PTAU24POSPCLINICALADvsNC_group-agecr.txt
          echo $MESH >> $MTDESIGNDIR/meshes_${side}_PTAU24POSPCLINICALADvsNC_group-agecr.txt
          fi

        elif [[ $dx == 4 ]]; then
          str="$id 0 1 $age"
          echo $str >> $MTDESIGNDIR/design_${side}_ABPOSEMCIvsNC_group-agecr.txt
          echo $MESH >> $MTDESIGNDIR/meshes_${side}_MCIvsNC_group-agecr.txt
        elif [[ $dx == 6 ]]; then
          str="$id 0 1 $age"
          echo $str >> $MTDESIGNDIR/design_${side}_ABPOSLMCIvsNC_group-agecr.txt
          echo $MESH >> $MTDESIGNDIR/meshes_${side}_ABPOSLMCIvsNC_group-agecr.txt
        elif [[ $dx == 8 ]]; then
          str="$id 0 1 $age"
          echo $str >> $MTDESIGNDIR/design_${side}_ABPOSADvsNC_group-agecr.txt
          echo $MESH >> $MTDESIGNDIR/meshes_${side}_ABPOSADvsNC_group-agecr.txt
        fi
      fi
    done
  done

  # echo contrast
  echo "1 -1 0" > \
    $MTDESIGNDIR/contrast_PCLINICALADvsNC_group-agecr_nc-pat.txt
  echo "-1 1 0" > \
    $MTDESIGNDIR/contrast_PCLINICALADvsNC_group-agecr_pat-nc.txt
  echo "1 -1 0" > \
    $MTDESIGNDIR/contrast_PTAU24POSPCLINICALADvsNC_group-agecr_nc-pat.txt
  echo "-1 1 0" > \
    $MTDESIGNDIR/contrast_PTAU24POSPCLINICALADvsNC_group-agecr_pat-nc.txt
  echo "1 -1 0" > \
    $MTDESIGNDIR/contrast_ABPOSEMCIvsNC_group-agecr_nc-pat.txt
  echo "-1 1 0" > \
    $MTDESIGNDIR/contrast_ABPOSEMCIvsNC_group-agecr_pat-nc.txt
  echo "1 -1 0" > \
    $MTDESIGNDIR/contrast_ABPOSLMCIvsNC_group-agecr_nc-pat.txt
  echo "-1 1 0" > \
    $MTDESIGNDIR/contrast_ABPOSLMCIvsNC_group-agecr_pat-nc.txt
  echo "1 -1 0" > \
    $MTDESIGNDIR/contrast_ABPOSADvsNC_group-agecr_nc-pat.txt
  echo "-1 1 0" > \
    $MTDESIGNDIR/contrast_ABPOSADvsNC_group-agecr_pat-nc.txt

  fi

  # submit jobs to run statistical analysis
  PREFIX=RTS${expid}
  if [[ 1  == 1 ]]; then
    for side in left right; do
      rm -rf $MTSTATSDIR/template_${side}
      for design in $(ls $MTDESIGNDIR | grep design | grep $side); do

        exp=$(echo $design | sed -e "s/^design_${side}_//" | sed -e "s/\.txt//")

        qsub -cwd -o $DUMPDIR -j y \
          -q all.q,basic.q \
          -l h_vmem=4.1G,s_vmem=4G \
          -N "${PREFIX}_${exp}_${side}_Atrophy" \
          $0 UTRegionalStats_meshglm_sub $exp $side
        sleep 0.1
      done
    done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1
  fi

}

function UTRegionalStats_meshglm_sub()
{
  exp=$1
  side=$2

  # Create the work directory for this analysis
  WORK=$MTSTATSDIR/UTtemplate_${side}/design_${exp}_MRG
  mkdir -p $WORK

  # Get the list of subjects
  DESIGNTXT=$MTDESIGNDIR/design_${side}_${exp}.txt
  SUBJ=$(cat $DESIGNTXT | awk '{print $1}')

  # Generate the design matrix for meshglm
  cat $DESIGNTXT | awk '{$1=""; print}' > $WORK/design_${side}.txt

  # combine thickness measurements
  MESHES=$(cat $MTDESIGNDIR/meshes_${side}_${exp}.txt)
  #MESHES=$(for id in $SUBJ; do \
  #  echo " $(ls $DATADIR/$id/*/ASHST1_MTLCORTEX_MSTTHK/GeoShoot/${side}/template_${grp}_to_${id}_${side}_GSShoot_MRG_thickmap.vtk) "; done)

  mesh_merge_arrays -r \
    $UTTEMPDIR/meshwarp/template_${side}_MRG.vtk \
    $WORK/thick_${side}_MRG.vtk ChangeAnnualizedPercentageAllTP $MESHES

  # Go through the list of contrasts
  for con in $(ls $MTDESIGNDIR | grep "contrast_${exp}_.*\.txt"); do

    # Get the suffix for this contrast
    suffix=$(echo $con | sed -e "s/^contrast_${exp}_//" -e "s/\.txt//")

    # Create the directory for this contrast
    CWORK=$WORK/contrast_${suffix}
    mkdir -p $CWORK

    FULLNM="design_${exp}_con_${suffix}"

    # Copy the contrast
    cp $MTDESIGNDIR/$con $CWORK/contrast.txt

    # Build the list of meshes to include
    MESHPARAM="-m $WORK/thick_${side}_MRG.vtk $CWORK/thickstat_${FULLNM}_${side}_MRG.vtk"
    meshglm $MESHPARAM \
      -g $WORK/design_${side}.txt $CWORK/contrast.txt \
      -a ChangeAnnualizedPercentageAllTP \
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
function SummarizeMultiTemp()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.tsv
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=10

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # Submit job to copy data
  if [[ 1 == 1 ]]; then
  PREFIX1=SM
  for ((i=2;i<=${N};i++)); do
  #for ((i=2000;i<=2010;i++)); do

    echo "$i submitted"
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX1}_${i}" \
         $0 SummarizeMultiTemp_sub $i $OUTTMPDIR
    sleep 0.01

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  fi

  # get all info
  # header
  blheader=$(cat -A $DEMOG_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  #sed -e 's/\^M\$//g'
  header="$header,ALOHA_Success,DateDiffFromBLScanDate"
  longheader=""
  for type in VOL MeanTHK MedianTHK; do
  for side in L R; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  done
  done
  done
  for temptype in MT UT; do
  for type in MeanTHKFITTED MedianTHKFITTED; do
  for side in L R; do
  for sub in ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_${temptype}_BL,${side}_${sub}_${type}_${temptype}_FU"
    longheader="$longheader,${side}_${sub}_${type}_${temptype}_ChangeAnnualized"
  done
  done
  done
  done
  header="$header$longheader,$blheader,END"
  echo $header > $OUTPUTDIR/longitudinal_multitemp.csv

  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/longitudinal_multitemp.csv
}

function SummarizeMultiTemp_sub()
{
  i=$1
  OUTTMPDIR=$2
  ALLINFO_TXT=$OUTPUTDIR/checkdata.tsv
  ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')

  # columns
  RIDCol=1
  IDCol=$(csvcol_tab.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol_tab.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT BLSANDATE)
  TYPECol=$(csvcol_tab.sh $ALLINFO_TXT AnalysisType)

  set +e

  # get information
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  vis2=$(echo $ROW | cut -f $Vis2Col -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  Type=$(echo $ROW | cut -f $TYPECol -d ",")
  SUBJALLDIR=$ALLDATADIR/$id
  SUBJDIR=$DATADIR/$id
  FUT1SRIMG=$SUBJDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz

  # baseline information
  blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
  blscandate=$(ReFormateDate $blscandate)
  DD=$(date -d "$blscandate" '+%d')
  MM=$(date -d "$blscandate" '+%m')
  YYYY=$(date -d "$blscandate" '+%Y')
  BLPREFIX="${YYYY}-${MM}-${DD}"
  BLSUBJDIR=$SUBJALLDIR/$BLPREFIX
  BLT1SRIMG=$BLSUBJDIR/${BLPREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
  BLASHST1DIR=$BLSUBJDIR/ASHST1
  BLASHST1LSEG=$BLASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
  BLASHST1RSEG=$BLASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

  SUBJALOHADIR=$SUBJDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}

  # baseline row
  BLROW=$(cat -A $DEMOG_TXT | grep $id | sed -e 's/\^M//g' | sed -e 's/\$//g')
  if [[ $BLROW == "" ]]; then
    BLROW=$(cat -A $DEMOG_TXT | head -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | awk  -F "," '{for(i=1;i<=NF-1;i++)printf ","}')
  fi

  # check if exist
  EXIST=0
  if [[ $Type != "T1Longitudinal" ]]; then
    EXIST="-1"
  elif [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $BLASHST1LSEG && -f $BLASHST1RSEG ]]; then
    EXIST=1
  else
    EXIST=0
  fi

  # compute date diff
  date_diff=$(( \
    ($(date -d $scandate +%s) - \
    $(date -d $blscandate +%s) )/(60*60*24) ))
  date_diff=$(echo ${date_diff#-})

  # go through all the fields
  RAWMEASURE=""
  LONGIMEASURE=""
  for type in volumes mean_thickness median_thickness; do
  for side in left right; do
    for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    SUBJALOHASUBDIR=$SUBJALOHADIR/results_${sub}
      if [[ $EXIST == 1 ]]; then
        # get volume and thickness
        if [[ -f $SUBJALOHASUBDIR/${type}_${side}.txt ]]; then
          MEA=$(cat $SUBJALOHASUBDIR/${type}_${side}.txt)
          MEA="${MEA// /,}"
          BL=$(echo $MEA | cut -d , -f 1)
          FU=$(echo $MEA | cut -d , -f 2)
          CHANGE=$(echo "scale=7;($BL-$FU)*100*365/$BL/$date_diff" | bc)
          LONGIMEASURE="$LONGIMEASURE,$CHANGE"

        else
          MEA=","
          EXIST=0
          LONGIMEASURE="$LONGIMEASURE,"
        fi
        RAWMEASURE="$RAWMEASURE,$MEA"
      else
        RAWMEASURE="$RAWMEASURE,,"
        LONGIMEASURE="$LONGIMEASURE,"
      fi
    done
  done
  done

  # get thickness measurements from the fitted mesh
  SUBJALOHASUBDIR=$SUBJALOHADIR/results_MTLCortexThkMultiTemp
  for temptype in Mult Uni; do
  for type in mean_thickness median_thickness; do
  for side in left right; do
    idx=3
    BLTXT=$SUBJALOHASUBDIR/${id}_${side}_${type}_${temptype}Temp_bl.txt
    FUTXT=$SUBJALOHASUBDIR/${id}_${side}_${type}_${temptype}Temp_fu.txt
    for sub in ERC BA35 BA36 PHC; do
      if [[ $EXIST == 1 ]]; then
        if [[ -f $BLTXT && -f $FUTXT ]]; then

          BL=$(cat $BLTXT | cut -d , -f $idx)
          FU=$(cat $FUTXT | cut -d , -f $idx)
          MEA="$BL,$FU"
          CHANGE=$(echo "scale=7;($BL-$FU)*100*365/$BL/$date_diff" | bc)
          LONGIMEASURE="$LONGIMEASURE,$CHANGE"

        else
          MEA=","
          EXIST=0
          LONGIMEASURE="$LONGIMEASURE,"
        fi
        RAWMEASURE="$RAWMEASURE,$MEA"
      else
        RAWMEASURE="$RAWMEASURE,,"
        LONGIMEASURE="$LONGIMEASURE,"
      fi
      idx=$((idx+1))
    done
  done
  done
  done


  echo "$ROW,$EXIST,$date_diff$RAWMEASURE$LONGIMEASURE,$BLROW,-1" >> $OUTTMPDIR/$(printf %04d $i)_${id}.csv

  set -e

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
