#//bin/bash
#$ -S /bin/bash
#set +e
set -e -x
#set -x -e

##############################################
# Setup environment
#source ~longxie/.bash_profile

# Software PATH
LXROOT=/home/lxie
ANTSPATH=$LXROOT/pkg/bin/antsbin/bin
C3DPATH=$LXROOT/pkg/c3d_tool/bin
FSLROOT=$LXROOT/pkg/fsl
FSLPATH=$FSLROOT/bin
PAULROOT=/project/hippogang_2/pauly/
#export ASHS_ROOT=/data/picsl/pauly/wolk/ashs
#export ASHS_ROOT=/data/picsl/longxie/pkg/ashs-fast
PAULROOT=/project/hippogang_2/pauly/
#export ASHS_ROOT=$PAULROOT/wolk/ashs-fast
export ASHS_ROOT=/project/hippogang_2/longxie/pkg/ashs/ashs-fast
export PATH=$PATH:$ASHS_ROOT/bin
#MATLAB_BIN=/share/apps/matlab/R2016a/bin/matlab
export ALOHA_ROOT=/project/hippogang_2/longxie/pkg/aloha
SRMATLABCODEDIR=$LXROOT/ASHS_PHC/SuperResolution/SRToolBox/matlabfunction
SHOOTWPRIORDIR=$LXROOT/LMTemplateMatching/PointSetGeodesicShooting/build_release
LMTOWARPDIR=$LXROOT/LMTemplateMatching/PointSetUtilities/build_release
#export FREESURFER_HOME=$CFNAPPS/freesurfer/6.0.0
#source $FREESURFER_HOME/SetUpFreeSurfer.sh
SRPATH=$LXROOT/pkg/PatchSuperResolution/build_release

#TEMPLATEDIR=/data/jag/wolk_group/ADNI_longitudinal-Templates/Normal/

##############################################
# Directories
ROOT=$LXROOT/ADNI2018/
CODEDIR=$ROOT/scripts
ANALYSISDIR=$ROOT/analysis_input
#ALLINFO_TXT=$ANALYSISDIR/MRI3TListWithNIFTIPath_11262018.tsv
#ALLINFO_TXT=/home/longxie/ADNI2018/info/T1Longitudinal_allTP/MRI3TListWithNIFTIPath_06132019.tsv
#ALLINFO_TXT=$ANALYSISDIR/MRILIST_T1T2_012_S_4128.tsv
LXROOT=/home/lxie/
ALLINFO_TXT=$LXROOT/ADNI2018/RefreshT1T2NIFTI_08022021/MRI3TListWithNIFTIPath_08022021.tsv
DATAROOT=/project/wolk_2/ADNI2018/
ALLDATADIR=$DATAROOT/dataset
OUTPUTDIR=$ROOT/output_runall_08022021
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
#DEMOG_TXT=$ANALYSISDIR/ADNI_T1_AllInfo_MTL_20180815.csv
DEMOG_TXT=$LXROOT/ADNI2018/output_all/ADNI_T1_measurement.csv

# TMPDIR
#if [[ ! $TMPDIR ]]; then

TMPDIR=$(mktemp -d /tmp/foo.XXXXXXXXXXX)
mkdir -p $TMPDIR
echo "Temporary directory: $TMPDIR"
#fi


#############################################
# parameters
ylthresh="0.0"
yhthresh="2.5"
smooth="6vox"
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
  # process scans of each timepoint indepedently
  ###############################

  # perform segmentation on baseline scans
  #ASHSSegSUB
  #CleanASHS

  # perform MTL thickness analysis
  #MSTMTUTThickness

  ###############################
  # generate ICV of the baseline timepoint
  ###############################

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

  ###############################
  # perform hotspot analysis
  ###############################
  # warp hotspot to baseline T1 SR space
  #WarpHotspotToBL

  # perform ALOHA one more time to estimate change
  #AlohaHotspot


  # perform analysis using variant-templates
  #WarpVTHotspotToBL

  # perform ALOHA one more time to estimate change for VT
  #AlohaVTHotspot

  ###############################
  # summary all information for various analyses
  ###############################
  # summarize result
  Summarize

  # summarize result with hotspot
  #SummarizeHotspot

  # summarize result for Robin whole brain covariance analysis
  #SummarizeRobin

  # summarize information for Mengjin's longitudinal project
  #SummarizeMengjin

  # summarize all results for comparisons (thickness and volume)
  #SummarizeComparisons


  # summarize volume and thickness measurements for Dave and Sandy's mismatch work (05132020)
  #SummarizeMismatchCohort


  # summarize baseline ASHS T1 info and QC info (ADNI_T1_AllInfo_MTL_20181025) for Mengjin's longitudinal paper
  #SummarizeBLMengjin


  ###############################
  # additional thickness measurements
  ###############################
  #MSTMTUTThicknessForMisMatchCohort


  ###############################
  # perform cleanup to save some space
  ###############################
  #CleanUpSpace







  ###############################
  # perform voxelwise longitudinal analysis
  # using the ANTs longitudinal pipeline 
  # ran by Robin (also Robin's spreadsheet, SummarizeRobin)
  ###############################
  #VoxelwiseLongiRobin
  #CheckVoxelwiseResult












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
  #N=3020

  # initialize tmp dir
  OUTTMPDIR=$OUTPUTDIR/tmp

  # Submit job to copy data
  PREFIX=LD${expid}
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  if [[ 0 == 1 ]]; then
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  set +e
  for ((i=2;i<=${N};i++)); do

    #file=$(ls $OUTTMPDIR/$(printf %04d $i)_*.csv)
    
    #if [[ $file == "" ]]; then

    pybatch.sh \
      -m "2G" \
      -n 1 \
      -N "${PREFIX}_${i}" \
      -o $DUMPDIR \
      $0 LinkData_sub $i $OUTTMPDIR
   
    #fi

    sleep 0.1

  done
  set -e

  # Wait for completion
  ALLJOBS=$(bjobs | wc -l)
  while [ $ALLJOBS -gt 0 ]; do
    sleep 60
    ALLJOBS=$(bjobs | wc -l)
  done
  #pybatch.sh \
  #  -n 1 \
  #  -w "${PREFIX}_*"
  fi

  # get all info
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  echo "$header,LinkExist,T1Orientation,T1Dimension,T1Resolution,LinkMatch,LinkT1NIFTI,AnalysisType"  > $OUTPUTDIR/checkdata.csv
  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/checkdata.csv
}

function LinkData_sub()
{
  i=$1
  OUTTMPDIR=$2
  ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  
  # columns
  RIDCol=1
  IDCol=$(csvcol_tab.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol_tab.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT BLSCANDATE)
  T1NIFTICol=$(csvcol_tab.sh $ALLINFO_TXT FINALT1NIFTI)
  #DATEDIFFTOBLCol=$(csvcol_tab.sh $ALLINFO_TXT DATEDIFFTOBL)

  # get information
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  vis2=$(echo $ROW | cut -f $Vis2Col -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX1="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
  blscandate=$(ReFormateDate $blscandate)
  T1NIFTI=$(echo $ROW | cut -f $T1NIFTICol -d ",")
  
  # correct for T1NIFTI because of the move
  T1NIFTI=$(echo $T1NIFTI | sed 's/\/data\/tesla-data\/PUBLIC\/ADNI2018/\/project\/wolk_1\/PUBLIC/g')

  datediff=$(( \
          ($(date -d $scandate +%s) - \
           $(date -d $blscandate +%s) )/(60*60*24) ))
  #datediff=$(echo $ROW | cut -f $DATEDIFFTOBLCol -d ",")
  SUBJALLDIR=$ALLDATADIR/$id
  #SUBJDIR=$DATADIR/$id

  # if file does not exist, skip
  if [[ $T1NIFTI == ""  ]]; then

    ROW="$ROW,,,,,,"

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
        | cut -d ':' -f 2 | cut -d ';' -f 3 | sed -e 's/,/_/g')

    else

      LINKNIFTI=""
      LINKMATCH=""
      EXIST=0
      ORI=""
      DIM=""
      RES=""

    fi

    if [[ $datediff == "0" ]]; then
      ROW="$ROW,$EXIST,$ORI,$DIM,$RES,$LINKMATCH,$LINKNIFTI,T1Baseline"
    else
      ROW="$ROW,$EXIST,$ORI,$DIM,$RES,$LINKMATCH,$LINKNIFTI,T1Longitudinal" 
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

  # columns
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)


  # Submit job to copy data
  PREFIX1=TN${expid}
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id

    # perform trim neck
    TPDIR=$SUBJALLDIR/${PREFIX}
    echo "TrimNeck Processing:$i,$rid,$PREFIX"

    if [[ -f $TPDIR/${PREFIX}_${id}_T1w.nii.gz && \
       ! -f $TPDIR/${PREFIX}_${id}_T1w_trim.nii.gz ]]; then
      echo "        Adding $rid $PREFIX "
      echo "TrimNeck: Adding $rid $PREFIX " >> $LOGDIR/trimneck_log.txt
      pybatch.sh \
        -m "4G" \
        -n 1 \
        -N "${PREFIX1}_${id}_${PREFIX}" \
        -o $DUMPDIR \
        $0 TrimNeck_sub \
           $TPDIR \
           $TPDIR/${PREFIX}_${id}_T1w.nii.gz \
           $TPDIR/${PREFIX}_${id}_T1w_trim.nii.gz
      sleep 0.1
    fi

  done

  # Wait for completion
  ALLJOBS=$(bjobs | wc -l)
  while [ $ALLJOBS -gt 0 ]; do
    sleep 60
    ALLJOBS=$(bjobs | wc -l)
  done
  #pybatch.sh \
  #  -n 1 \
  #  -w "${PREFIX1}_*"
}

######################################################
function TrimNeck_sub()
{
  TPDIR=$1
  T1=$2
  T1TRIM=$3

  if [[ ! -f $T1TRIM ]]; then

    /home/lxie/pkg/trim_neck_rf.sh \
      $T1 $T1TRIM

    echo "Using /home/lxie/pkg/trim_neck_rf.sh" > \
      $TPDIR/neck_trimed_version.txt

  fi
}

######################################################
function SR()
{
  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11

  # columns
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

  # Submit job to copy data
  PREFIX1=SR${expid}
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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


    if [[ -f $TPDIR/${PREFIX}_${id}_T1w_trim.nii.gz && ! -f $TPDIR/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz && ! -f $SUBJALLDIR/${PREFIX}/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz ]]; then

      echo "   Adding $rid $PREFIX "
      echo "SR: Adding $rid $PREFIX " >> $LOGDIR/SR_log.txt
      pybatch.sh \
        -m "4G" \
        -n 1 \
        -N "${PREFIX1}_${id}_${PREFIX}" \
        -o $DUMPDIR \
        $0 SR_sub $TPDIR \
           $TPDIR/${PREFIX}_${id}_T1w_trim.nii.gz \
           $TPDIR/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
      sleep 0.1

    #elif [[ -f $TPDIR/${PREFIX}_${id}_T1w_trim.nii.gz && ! -f $TPDIR/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz && -f $SUBJDIR/${PREFIX}/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz ]]; then
    #  echo "   Copying $rid $PREFIX from dataset_local"
    #  echo "SR: Copying $rid $PREFIX from dataset_local" >> $LOGDIR/SR_log.txt
    #  cp $SUBJDIR/${PREFIX}/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz $TPDIR/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
    fi

  done

  # Wait for completion
  ALLJOBS=$(bjobs | wc -l)
  while [ $ALLJOBS -gt 0 ]; do
    sleep 60
    ALLJOBS=$(bjobs | wc -l)
  done
  #pybatch.sh \
  #  -n 1 \
  #  -w "${PREFIX}_*"
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
  ATLAS=/home/lxie/ASHS_atlases/PMC_3TT1_atlas
  Nrun=5

  # Submit job to copy data
  PREFIX1=SR${expid}
  for ((i=2;i<=${N};i++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # columns
    RIDCol=1
    IDCol=$(csvcol.sh $ALLINFO_TXT ID)
    Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
    ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
    TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

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
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    echo "ASHS T1 Processing: $i,$rid,$PREFIX"
    if [[ $Type != "T1Baseline" && $Type != "T1Longitudinal" ]]; then
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
        -l -s 1-7 \
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
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  
  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

  # Submit job to copy data
  PREFIX=ASHSSeg${expid}
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    echo "Checking ASHS T1: $i,$rid,$PREFIX"
    if [[ $Type != "T1Baseline" && $Type != "T1Longitudinal" ]]; then
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

######################################################
function MSTMTUTThickness()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=10
  N_begin=2

  echo "" >> $LOGDIR/MSTThkASHST1_log.txt
  echo $(date) >> $LOGDIR/MSTThkASHST1_log.txt

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

  # Submit job to copy data
  JPREFIX=MSTTHK${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    #SUBJDIR=$DATADIR/$id

    # skip if it is longitudinal
    echo "MST thickness ASHS T1 Processing: $i,$rid,$PREFIX $Type"
    if [[ $Type != "T1Baseline" ]]; then
      continue
    fi

    SUBJASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
    SUBJMSTTHKDIR=$SUBJALLDIR/$PREFIX/ASHST1_MTLCORTEX_MSTTHK
    #SUBJMSTDATADIR=$SUBJMSTTHKDIR/data
    #SUBJMSTREGATLASDIR=$SUBJMSTTHKDIR/RegToAtlases
    #SUBJMSTVTREGDIR=$SUBJMSTTHKDIR/RegToInitTemp
    #SUBJMSTUTREGDIR=$SUBJMSTTHKDIR/RegToUT

    # temporary only process ADNI GO and 2
    #if [[ $rid -lt 2000 || $rid -gt 5999 ]]; then
    #  continue
    #fi


    # this is temporary, prioratise the ones needed
    #if [[ ! -d $SUBJMSTTHKDIR/GeoShoot ]]; then
    #  continue
    #fi
    #SNAPROW=$(cat $ANALYSISDIR/ADNIBL_Merge20190605_21_SNAP.csv | grep $id | grep $scandate)
    #DX=$(echo $SNAPROW | cut -f 5 -d ,)
    #GRP=$(echo $SNAPROW | cut -f 6 -d ,)
    #if [[ $DX -gt 2 || $GRP == 1 ]]; then
    #  continue
    #fi

    # report status
    echo "   Adding $rid $PREFIX $Type"
    echo "MST thickness ASHS T1: Adding $rid $PREFIX $Type" >> $LOGDIR/MSTThkASHST1_log.txt

    # submit jobs
    for side in left right; do
      AUTOSEG=$SUBJASHST1DIR/final/${id}_${side}_lfseg_heur.nii.gz
      if [[ ! -f $SUBJMSTTHKDIR/work_${side}/RegToUT/UTemp/${id}_${PREFIX}_${side}_totempWarp.nii.gz ]]; then
      pybatch.sh \
        -m "8G" \
        -n 1 \
        -N "${JPREFIX}_${id}_${PREFIX}_${side}" \
        -o $DUMPDIR \
        $0 MSTMTUTThickness_sub $id $PREFIX $side $AUTOSEG $SUBJMSTTHKDIR
       else
         echo "skip!"
       fi
    done

  done

  # Wait for completion
  ALLJOBS=$(bjobs | wc -l)
  while [ $ALLJOBS -ge 100 ]; do
    sleep 60
    ALLJOBS=$(bjobs | wc -l)
  done
  pybatch.sh \
    -n 1 \
    -w "${PREFIX}_*"
}

function MSTMTUTThickness_sub()
{
  id=$1
  PREFIX=$2
  side=$3
  AUTOSEG=$4
  OUTDIR=$5
  OUTCSV=$OUTDIR/${id}_${PREFIX}_${side}_thickness.csv
  
  if [[ ! -f $OUTCSV ]]; then
    mkdir -p $OUTDIR
    /home/lxie/ASHS_T1/pipeline_package/MultiTempThkMTLWithHippo/runpipeline_MultiAndUnifiTemplate_single.sh \
      ${id}_${PREFIX} $AUTOSEG $side $OUTDIR
  fi
}

##################################################
function ASHSICVSegSUB()
{
  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=11
  ICVATLAS=/home/lxie/ASHS_atlases/ICVatlas_3TT1
  Nrun=20

  # Submit job to copy data
  PREFIX=ASHSICVSeg${expid}
  for ((i=2;i<=${N};i++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

    # columns
    RIDCol=1
    IDCol=$(csvcol.sh $ALLINFO_TXT ID)
    Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
    ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
    BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
    TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

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
        -l -s 1-7 \
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
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)


  # Submit job to copy data
  PREFIX=ASHSSeg${expid}
  for ((i=2;i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)  
  #N=10
  N_begin=2

  # columns   RIDCol=1
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

  # Submit job to copy data
  JPREFIX=ALOHA${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    #SUBJDIR=$DATADIR/$id

    # skip if it is not longitudinal
    echo "ALOHA Processing: $i,$rid,$PREFIX $Type"
    if [[ $Type != "T1Longitudinal" ]]; then
      continue
    fi

    # perform SR
    FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
    FUASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
    FUASHST1LSEG=$FUASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
    FUASHST1RSEG=$FUASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

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

    # start running ALOHA  forward
    SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
    if [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $BLASHST1LSEG && -f $BLASHST1RSEG && ! -d $SUBJALOHADIR ]]; then

      #if [[ -d $SUBJDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX} ]]; then

      #  echo "   Copying $rid $PREFIX from dataset_local"
      #  echo "ALOHA: Copying $rid $PREFIX from dataset_local" >> $LOGDIR/ALOHA_log.txt

      #  cp -r $SUBJDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX} \
      #     $SUBJALLDIR/${PREFIX}/

      #else

      echo "   Adding $rid $PREFIX $Type"
      echo "ALOHA: Adding $rid $PREFIX $Type" >> $LOGDIR/ALOHA_log.txt

      pybatch.sh \
        -m "5G" \
        -n 1 \
        -N "${JPREFIX}_${i}_${rid}_${PREFIX}" \
        -o $DUMPDIR \
        $0 AlohaMTL_sub \
           $BLT1SRIMG \
           $FUT1SRIMG \
           $BLASHST1LSEG \
           $BLASHST1RSEG \
           $SUBJALOHADIR
        sleep 0.1

      #fi
    fi

    # start running ALOHA  forward
    SUBJALOHADIR=$SUBJALLDIR/${BLPREFIX}/AlohaMTL_ASHST1_${PREFIX}
    if [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $FUASHST1LSEG && -f $FUASHST1RSEG && ! -d $SUBJALOHADIR ]]; then

      echo "   Adding $rid $PREFIX $Type reverse"
      echo "ALOHA: Adding $rid $PREFIX $Type reverse" >> $LOGDIR/ALOHA_log.txt

      pybatch.sh \
        -m "5G" \
        -n 1 \
        -N "${JPREFIX}_${i}_${rid}_${PREFIX}_rev" \
        -o $DUMPDIR \
        $0 AlohaMTL_sub \
           $FUT1SRIMG \
           $BLT1SRIMG \
           $FUASHST1LSEG \
           $FUASHST1RSEG \
           $SUBJALOHADIR
      sleep 0.1


    fi

  done

  # Wait for completion
  ALLJOBS=$(bjobs | wc -l)
  while [ $ALLJOBS -gt 0 ]; do
    sleep 60
    ALLJOBS=$(bjobs | wc -l)
  done
  pybatch.sh \
    -n 1 \
    -w "${JPREFIX}_*"
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
<<COMMENT
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
COMMENT
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
function ComputePWNCC()
{
  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=3020
  N_begin=2

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # header
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
  echo "$header,ALOHA_Success,ALOHA_run_type,NCOR_L,NCOR_R" > $OUTPUTDIR/QC_ALOHA.csv

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)


  # Submit job to copy data
  JPREFIX=QCALOHA${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    #SUBJDIR=$DATADIR/$id

    # skip if it is not longitudinal
    echo "ALOHA QC Processing: $i,$rid,$PREFIX"
    if [[ $Type != "T1Longitudinal" ]]; then
      echo "$ROW,-1,forward,," > $OUTTMPDIR/$(printf %04d $i)_0.csv
      echo "$ROW,-1,reverse,," > $OUTTMPDIR/$(printf %04d $i)_1.csv
      continue
    else
      echo "$ROW,0,forward,," > $OUTTMPDIR/$(printf %04d $i)_0.csv
      echo "$ROW,0,reverse,," > $OUTTMPDIR/$(printf %04d $i)_1.csv
    fi

    # perform SR
    FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
    FUASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
    FUASHST1LSEG=$FUASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
    FUASHST1RSEG=$FUASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

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

    # compute measure: forward
    SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
    if [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $BLASHST1LSEG && -f $BLASHST1RSEG ]]; then

      echo "   Adding $rid $PREFIX forward"
      echo "ALOHA QC: Adding $rid $PREFIX forward" >> $LOGDIR/ALOHA_QC_log.txt

        qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${i}_${rid}_${PREFIX}_forward" $0 \
               ComputePWNCC_sub \
               $i forward $SUBJALOHADIR $OUTTMPDIR
          sleep 0.1

    fi

    
    # compute measure: backward
    SUBJALOHADIR=$SUBJALLDIR/${BLPREFIX}/AlohaMTL_ASHST1_${PREFIX}
    if [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $FUASHST1LSEG && -f $FUASHST1RSEG ]]; then

      echo "   Adding $rid $PREFIX reverse"
      echo "ALOHA QC: Adding $rid $PREFIX reverse" >> $LOGDIR/ALOHA_QC_log.txt

        qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${i}_${rid}_${PREFIX}_reverse" $0 \
               ComputePWNCC_sub \
               $i reverse $SUBJALOHADIR $OUTTMPDIR
          sleep 0.1

    fi
  done


  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1

  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/QC_ALOHA.csv
}


function ComputePWNCC_sub()
{
  i=$1
  type=$2
  SUBJALOHADIR=$3
  OUTTMPDIR=$4
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  if [[ $type == "forward" ]]; then
    postfix=0
  else
    postfix=1
  fi
  ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
      -dilate 1 20x20x20vox \
      -o $TMPDIR/blmptrim_seg_${side}_tohw_mask20vox.nii.gz

    NCOR=$(c3d $SUBJALOHADIR/deformable/blmptrim_${side}_to_hw.nii.gz \
      $TMPDIR/blmptrim_seg_${side}_tohw_mask20vox.nii.gz \
      -multiply \
      $SUBJALOHADIR/deformable/fumptrim_om_${side}to_hw.nii.gz \
      $TMPDIR/blmptrim_seg_${side}_tohw_mask20vox.nii.gz \
      -multiply \
      -ncor | cut -d = -f 2)

    echo $NCOR > $SUBJALOHADIR/qc_ncor_${side}.txt
    NCORALL="$NCORALL,$NCOR"

  done
  echo "$ROW,1,$type$NCORALL" > $OUTTMPDIR/$(printf %04d $i)_${postfix}.csv
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
    ALOHA_QC_Select('$INCSV', '$OUTCSV', 0.7, 0.05);
MATCODE
  fi
  ALLINFO_TXT=$OUTCSV

  # Load id number
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=3020
  N_begin=2

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)
  NCORLCol=$(csvcol.sh $ALLINFO_TXT NCOR_L)
  NCORRCol=$(csvcol.sh $ALLINFO_TXT NCOR_R)
  QCSELLCol=$(csvcol.sh $ALLINFO_TXT QC_SEL_L)
  QCSELRCol=$(csvcol.sh $ALLINFO_TXT QC_SEL_R)
  RUNTYPECol=$(csvcol.sh $ALLINFO_TXT ALOHA_run_type)

  # Submit job to copy data
  JPREFIX=QCALOHA${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id

    # perform SR
    FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
    FUASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
    FUASHST1LSEG=$FUASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
    FUASHST1RSEG=$FUASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

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

    # compute measure: forward
    runtype=$(echo $ROW | cut -f $RUNTYPECol -d ",")
    if [[ $runtype == "forward" ]]; then
      SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
    else
      SUBJALOHADIR=$SUBJALLDIR/${BLPREFIX}/AlohaMTL_ASHST1_${PREFIX}
    fi

    for side in left right; do

      echo "ALOHA QC Processing: $i, $id, $side, $scandate, $runtype"
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
          -res 6

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
function WarpHotspotToBL()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  #N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  N=1001
  N_begin=1000

  # log
  echo "" >> $LOGDIR/warphotspot_log.txt
  echo $(date) >> $LOGDIR/warphotspot_log.txt

  # columns   RIDCol=1
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

  # Submit job to copy data
  JPREFIX=WHSP${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    #SUBJDIR=$DATADIR/$id

    # skip if it is longitudinal
    echo "Warp Hotspot Processing: $i,$rid,$PREFIX $Type"
    if [[ $Type != "T1Baseline" ]]; then
      continue
    fi

    # temporary only process ADNI GO and 2
    if [[ $rid -lt 2000 || $rid -gt 5999 ]]; then
      continue
    fi

    SUBJASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
    SUBJMSTTHKDIR=$SUBJALLDIR/$PREFIX/ASHST1_MTLCORTEX_MSTTHK
    SUBJMSTDATADIR=$SUBJMSTTHKDIR/data
    SUBJMSTREGATLASDIR=$SUBJMSTTHKDIR/RegToAtlases
    SUBJMSTVTREGDIR=$SUBJMSTTHKDIR/RegToInitTemp

    # report status
    echo "   Adding $rid $PREFIX $Type"
    echo "Warp Hotspot: Adding $rid $PREFIX $Type" >> $LOGDIR/warphotspot_log.txt

    # submit jobs
    for side in left right; do

      if [[ -d $SUBJMSTTHKDIR/work_${side}/RegToUT/UTemp ]]; then

        qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q \
           -l h_vmem=8.1G,s_vmem=8G \
           -N "${JPREFIX}_${id}_${PREFIX}_${side}" $0 \
           WarpHotspotToBL_sub $id $PREFIX $side $SUBJALLDIR

        # report status
        echo "   Adding $rid $PREFIX $Type"
        echo "Warp Hotspot: Adding $rid $PREFIX $Type" >> $LOGDIR/warphotspot_log.txt


      else
        # report status
        echo "   Skipping $rid $PREFIX $Type"
        echo "Warp Hotspot: Skipping $rid $PREFIX $Type" >> $LOGDIR/warphotspot_log.txt
      fi   
        

    done

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}

function WarpHotspotToBL_sub()
{
  id=$1
  PREFIX=$2
  side=$3
  SUBJALLDIR=$4
  T1SR=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
  SUBJASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
  SUBJMSTTHKDIR=$SUBJALLDIR/$PREFIX/ASHST1_MTLCORTEX_MSTTHK
  SUBJMSTTHKWORKDIR=$SUBJMSTTHKDIR/work_${side}
  SUBJDATADIR=$SUBJMSTTHKWORKDIR/data
  SUBJREGUTTEMPDIR=$SUBJMSTTHKWORKDIR/RegToUT 
  SUBJHOTSPOTDIR=$SUBJMSTTHKDIR/hotspot
  mkdir -p $SUBJHOTSPOTDIR

  # Reference space (root node in cm-rep space)
  REFSPACE=$SUBJHOTSPOTDIR/refspace_${side}.nii.gz
  if [ ! -f $REFSPACE ]; then
  PAD=60
  c3d $SUBJDATADIR/${id}_${PREFIX}_${side}_seg_orig.nii.gz \
    -pad ${PAD}x${PAD}x${PAD}vox ${PAD}x${PAD}x${PAD}vox 0 \
    -o $TMPDIR/${id}_${PREFIX}_${side}_seg_orig.nii.gz

  MSTUTTEMPDIR=/home/longxie/ASHS_T1/pipeline_package/MultiTempThkMTLWithHippo/template/GSUTemplate/
  /data/picsl/longxie/pkg/bin/antsbin/bin/AverageImages 3 \
    $REFSPACE 0 \
    $MSTUTTEMPDIR/InitTemp/template_1/work/iter_04/template_1_seg.nii.gz \
    $MSTUTTEMPDIR/gshoot/template_1/refspace_1.nii.gz \
    $TMPDIR/${id}_${PREFIX}_${side}_seg_orig.nii.gz

  c3d $REFSPACE \
    -thresh 0.0001 inf 1 0 \
    -trim 20vox \
    -resample-mm 0.4x0.4x0.4mm \
    -o $REFSPACE
  fi

  # warp hotspot mesh to subject space
  WARPCHAIN=$(cat $SUBJREGUTTEMPDIR/UTemp/chain_unwarp_to_final_${side}.txt)
  HOTSPOTDIR=/home/longxie/ASHS_T1/exvivo_atlas_hotspot/t1_template
  for label in 2 3; do
    hotspot=$HOTSPOTDIR/ev_to_ivtemplate_hotspot_resliced_label${label}.vtk
    TGHSP=$SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_hotspot_label${label}.vtk
    TGHNII=$SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_hotspot_label${label}.nii.gz

    if [[ ! -f $TGHNII ]]; then
    greedy -d 3 \
      -rf $REFSPACE \
      -rs $hotspot $TGHSP \
      -r $WARPCHAIN \
         $SUBJDATADIR/${id}_${PREFIX}_${side}_to_MSTInitTemp_mlaffine.txt
      
    # convert mesh to nifti 
    mesh2img -f -vtk $TGHSP \
      -a 0.3 0.3 0.3 4 \
      $TGHNII
    c3d $T1SR \
      $TGHNII \
      -int 0 -reslice-identity \
      -o $TGHNII
    fi
  done

  c3d $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_hotspot_label2.nii.gz \
    -scale 2 \
    $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_hotspot_label3.nii.gz \
    -scale 3 \
    -add \
    -replace 5 2 \
    -o $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_hotspot.nii.gz

<<COMMENT
  # compose the deformation from native space to unified template space
  if [[ ! -f $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_native_to_UTTemp_warp.nii.gz ]]; then

  WARPCHAIN=($(cat $SUBJREGUTTEMPDIR/UTemp/chain_unwarp_to_final_${side}.txt))
  CWARPCHAIN=""
  for ((j=0;j<${#WARPCHAIN[*]};j++)); do
    str=${WARPCHAIN[j]}
    comp1=$(echo $str | cut -d , -f 1)
    comp2=$(echo $str | cut -d , -f 2)
    if [[ $comp2 == "-1" ]]; then
      prefix="-i"
    else
      prefix=""
    fi

    filename=`basename $comp1`
    fileext=${filename##*.}
    if [[ $fileext == "mat" || $fileext == "txt" ]]; then
      c3d_affine_tool $comp1 \
        -oitk $TMPDIR/${filename}_affine.txt
      outstr=$TMPDIR/${filename}_affine.txt
    else
      outstr=$comp1
    fi

    CWARPCHAIN="$CWARPCHAIN $prefix $outstr"
  done

  c3d_affine_tool $SUBJDATADIR/${id}_${PREFIX}_${side}_to_MSTInitTemp_mlaffine.txt \
    -oitk $TMPDIR/${id}_${PREFIX}_${side}_to_MSTInitTemp_mlaffine.txt

  $ANTSPATH/ComposeMultiTransform 3 \
    $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_native_to_UTTemp_warp.nii.gz \
    -R $REFSPACE \
    $CWARPCHAIN \
    $TMPDIR/${id}_${PREFIX}_${side}_to_MSTInitTemp_mlaffine.txt
  
  fi

  # compute inverse warp
  if [[ ! -f $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_native_to_UTTemp_invwarp.nii.gz ]]; then
  greedy -d 3 -iw $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_native_to_UTTemp_warp.nii.gz \
    $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_native_to_UTTemp_invwarp.nii.gz
  fi
 
  # warp template hotspt ROI to subject space
  HOTSPOTDIR=/data/jux/sravikumar/atlasPHG2019/map_hotspots/register_invivotemp/t1_template
  greedy -d 3 \
    -rf $T1SR \
    -ri LABEL 0.3vox \
    -rm $HOTSPOTDIR/ev_to_ivtemplate_hotspot_resliced.nii.gz \
        $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_ev_to_ivtemplate_hotspot_native.nii.gz \
    -r $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_native_to_UTTemp_invwarp.nii.gz
COMMENT


}

######################################################
function WarpVTHotspotToBL()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=1110
  N_begin=1100

  # log
  echo "" >> $LOGDIR/warpVThotspot_log.txt
  echo $(date) >> $LOGDIR/warpVThotspot_log.txt

  # columns   RIDCol=1
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

  # Submit job to copy data
  JPREFIX=WHSP${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    #SUBJDIR=$DATADIR/$id

    # skip if it is longitudinal
    echo "Warp Hotspot Processing: $i,$rid,$PREFIX $Type"
    if [[ $Type != "T1Baseline" ]]; then
      continue
    fi

    # temporary only process ADNI GO and 2
    if [[ $rid -lt 2000 || $rid -gt 5999 ]]; then
      continue
    fi

    SUBJASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
    SUBJMSTTHKDIR=$SUBJALLDIR/$PREFIX/ASHST1_MTLCORTEX_MSTTHK
    SUBJMSTDATADIR=$SUBJMSTTHKDIR/data
    SUBJMSTREGATLASDIR=$SUBJMSTTHKDIR/RegToAtlases
    SUBJMSTVTREGDIR=$SUBJMSTTHKDIR/RegToInitTemp

    # report status
    echo "   Adding $rid $PREFIX $Type"
    echo "Warp Hotspot: Adding $rid $PREFIX $Type" >> $LOGDIR/warpVThotspot_log.txt

    # submit jobs
    for side in left right; do

      if [[ -d $SUBJMSTTHKDIR/work_${side}/RegToInitTemp/inittemp ]]; then

        #qsub -cwd -o $DUMPDIR -j y \
        #   -q all.q \
        #   -l h_vmem=8.1G,s_vmem=8G \
        #   -N "${JPREFIX}_${id}_${PREFIX}_${side}" $0 \
        #   WarpVTHotspotToBL_sub $id $PREFIX $side $SUBJALLDIR

        pybatch.sh \
          -m "8G" \
          -n 1 \
          -N "${JPREFIX}_${id}_${PREFIX}_${side}" \
          -o $DUMPDIR \
          $0 WarpVTHotspotToBL_sub $id $PREFIX $side $SUBJALLDIR

        # report status
        echo "   Adding $rid $PREFIX $Type"
        echo "   Warp Hotspot: Adding $rid $PREFIX $Type $side" >> $LOGDIR/warpVThotspot_log.txt


      else
        # report status
        echo "   Skipping $rid $PREFIX $Type"
        echo "   Warp Hotspot: Skipping $rid $PREFIX $Type $side" >> $LOGDIR/warpVThotspot_log.txt
      fi


    done

  done

  # Wait for completion
  pybatch.sh \
    -n 1 \
    -w "${JPREFIX}_*"

  #qsub -cwd -o $DUMPDIR -j y \
  #     -q all.q \
  #     -hold_jid "${JPREFIX}_*" -sync y -b y \
  #     sleep 1
}

function WarpVTHotspotToBL_sub()
{
  id=$1
  PREFIX=$2
  side=$3
  SUBJALLDIR=$4
  T1SR=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
  SUBJASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
  SUBJMSTTHKDIR=$SUBJALLDIR/$PREFIX/ASHST1_MTLCORTEX_MSTTHK
  SUBJMSTTHKWORKDIR=$SUBJMSTTHKDIR/work_${side}
  SUBJDATADIR=$SUBJMSTTHKWORKDIR/data
  SUBJREGVTTEMPDIR=$SUBJMSTTHKWORKDIR/RegToInitTemp
  SUBJMEMBERSHIPDIR=$SUBJMSTTHKWORKDIR/membership
  #SUBJREGUTTEMPDIR=$SUBJMSTTHKWORKDIR/RegToUT
  SUBJHOTSPOTDIR=$SUBJMSTTHKDIR/VThotspot06032021
  mkdir -p $SUBJHOTSPOTDIR

  # group
  grp=$(cat $SUBJMEMBERSHIPDIR/autogroup_${side}.txt | head -n 1 | cut -d " " -f 1)
  MSTVTTEMPDIR=$LXROOT/ASHS_T1/pipeline_package/MultiTempThkMTLWithHippo/template/GSTemplate/

  # Reference space (root node in cm-rep space)
  REFSPACE=$SUBJHOTSPOTDIR/refspace_${side}.nii.gz
  if [ ! -f $REFSPACE ]; then
  PAD=60
  c3d $SUBJDATADIR/${id}_${PREFIX}_${side}_seg_orig.nii.gz \
    -pad ${PAD}x${PAD}x${PAD}vox ${PAD}x${PAD}x${PAD}vox 0 \
    -o $TMPDIR/${id}_${PREFIX}_${side}_seg_orig.nii.gz

  MSTUTTEMPDIR=$LXROOT/ASHS_T1/pipeline_package/MultiTempThkMTLWithHippo/template/GSUTemplate/
  $LXROOT/pkg/bin/antsbin/bin/AverageImages 3 \
    $REFSPACE 0 \
    $MSTVTTEMPDIR/InitTemp/template_${grp}/work/iter_04/template_${grp}_seg.nii.gz \
    $MSTVTTEMPDIR/gshoot/template_${grp}/refspace_${grp}.nii.gz \
    $TMPDIR/${id}_${PREFIX}_${side}_seg_orig.nii.gz

  c3d $REFSPACE \
    -thresh 0.0001 inf 1 0 \
    -trim 20vox \
    -resample-mm 0.4x0.4x0.4mm \
    -o $REFSPACE
  fi

  # warp hotspot mesh to subject space
  WARPCHAIN=$(cat $SUBJREGVTTEMPDIR/inittemp/chain_unwarp_to_final_${side}.txt \
    | sed 's/\/data\/jux\/wolk_group/\/home\/lxie/g' \
    | sed 's/\/home\/longxie/\/home\/lxie/g')
  HOTSPOTDIR=$LXROOT/exvivo/exvivo_hotspot_ivtemplate06032021/
  for label in 2 3 4; do
    hotspot=$HOTSPOTDIR/ev_to_ivtemplate_allhotspots_resliced_template${grp}_label${label}.vtk
    TGHSP=$SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_hotspot_label${label}.vtk
    TGHNII=$SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_hotspot_label${label}.nii.gz

    if [[ ! -f $TGHNII ]]; then
    greedy -d 3 \
      -rf $REFSPACE \
      -rs $hotspot $TGHSP \
      -r $WARPCHAIN \
         $SUBJDATADIR/${id}_${PREFIX}_${side}_to_MSTInitTemp_mlaffine.txt

    # convert mesh to nifti
    mesh2img -f -vtk $TGHSP \
      -a 0.3 0.3 0.3 4 \
      $TGHNII
    c3d $T1SR \
      $TGHNII \
      -int 0 -reslice-identity \
      -o $TGHNII
    fi
  done

  c3d $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_hotspot_label2.nii.gz \
    -scale 2 \
    $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_hotspot_label3.nii.gz \
    -scale 3 \
    -add \
    -replace 5 2 \
    $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_hotspot_label4.nii.gz \
    -scale 4 \
    -add \
    -replace 6 2 -replace 7 3 \
    -o $SUBJHOTSPOTDIR/${id}_${PREFIX}_${side}_hotspot.nii.gz
}



######################################################
function AlohaHotspot()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=1075
  N_begin=1068
  #N_begin=2

  # log
  echo "" >> $LOGDIR/ALOHAHotspot_log.txt
  echo $(date) >> $LOGDIR/ALOHAHotspot_log.txt

  # columns   RIDCol=1
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

  # Submit job to copy data
  JPREFIX=ALOHA${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    #SUBJDIR=$DATADIR/$id

    # skip if it is not longitudinal
    echo "ALOHA Processing: $i,$rid,$PREFIX $Type"
    if [[ $Type != "T1Longitudinal" ]]; then
      continue
    fi

    # perform SR
    FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
    FUASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
    FUASHST1LSEG=$FUASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
    FUASHST1RSEG=$FUASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

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
    BLMSTTHKDIR=$BLSUBJDIR/ASHST1_MTLCORTEX_MSTTHK
    BLLHOTSPOT=$BLMSTTHKDIR/hotspot/${id}_${BLPREFIX}_left_hotspot.nii.gz
    BLRHOTSPOT=$BLMSTTHKDIR/hotspot/${id}_${BLPREFIX}_right_hotspot.nii.gz

    # start running ALOHA  forward
    SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaHotspot_ASHST1_${BLPREFIX}
    if [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $BLLHOTSPOT && -f $BLRHOTSPOT && ! -d $SUBJALOHADIR ]]; then
   
      echo "   Adding $rid $PREFIX $Type"
      echo "ALOHA hotspot: Adding $rid $PREFIX $Type" >> $LOGDIR/ALOHAHotspot_log.txt

      qsub -cwd -o $DUMPDIR -j y \
           -p -1023 \
           -q all.q \
           -l h_vmem=5.1G,s_vmem=5G \
           -N "${JPREFIX}_${i}_${rid}_${PREFIX}" \
           $0 AlohaHotspot_sub \
               $BLT1SRIMG \
               $FUT1SRIMG \
               $BLLHOTSPOT \
               $BLRHOTSPOT \
               $SUBJALOHADIR
        sleep 0.1


      #fi
    fi

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}

function AlohaHotspot_sub()
{
  BLT1SRIMG=$1
  YT1SRIMG=$2
  BLLHOTSPOT=$3
  BLRHOTSPOT=$4
  SUBJALOHADIR=$5
  mkdir -p $SUBJALOHADIR

  if [[ ! -f $SUBJALOHADIR/results/volumes_left.txt || ! -f $SUBJALOHADIR/results/volumes_right.txt ]]; then

  c3d $BLLHOTSPOT \
     -thresh 2 4 1 0 \
     -o $SUBJALOHADIR/bl_seg_left.nii.gz

  c3d $BLRHOTSPOT \
     -thresh 2 4 1 0 \
     -o $SUBJALOHADIR/bl_seg_right.nii.gz

  $ALOHA_ROOT/scripts/aloha_main.sh \
    -b $BLT1SRIMG \
    -f $YT1SRIMG \
    -r $SUBJALOHADIR/bl_seg_left.nii.gz \
    -s $SUBJALOHADIR/bl_seg_right.nii.gz \
    -w $SUBJALOHADIR \
    -t 1-4

  fi

  for sub in label2 label3 label4; do

    if [[ ! -f $SUBJALOHADIR/results_${sub}/volumes_left.txt || ! -f $SUBJALOHADIR/results_${sub}/volumes_right.txt ]]; then

      if [[ $sub == "label2" ]]; then
        label=(2 2)
      elif [[ $sub == "label3" ]]; then
        label=(3 3)
      elif [[ $sub == "label4" ]]; then
        label=(4 4)
      fi

      SUBJALOHASUBDIR=$TMPDIR/aloha_${sub}
      mkdir -p $SUBJALOHASUBDIR/results
      ln -sf $SUBJALOHADIR/deformable $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/init $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/global $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/dump $SUBJALOHASUBDIR
      ln -sf $SUBJALOHADIR/final $SUBJALOHASUBDIR

      c3d $BLLHOTSPOT \
        -thresh ${label[0]} ${label[1]} 1 0 \
        -o $SUBJALOHASUBDIR/bl_seg_left.nii.gz

      c3d $BLRHOTSPOT \
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
<<COMMENT
    for side in left right; do

      if [[ ! -f $SUBJALOHADIR/results_${sub}/mean_thickness_${side}.txt || ! -f $SUBJALOHADIR/results_${sub}/median_thickness_${side}.txt ]]; then

        # measure thickness of bl mesh
        cmrep_vskel -Q $PAULROOT/bin/qvoronoi \
          -T $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_thickmap.vtk \
          -p 1.2 -e 6 \
          $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}.vtk \
          $SUBJALOHADIR/results_${sub}/skel_blmptrim_seg_${side}.vtk

        # measure thickness of fu mesh
        cmrep_vskel -Q $PAULROOT/bin/qvoronoi \
          -T $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_warped_to_futrim_om_thickmap.vtk \
          -p 1.2 -e 6 \
          $SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_warped_to_futrim_om.vtk \
          $SUBJALOHADIR/results_${sub}/skel_blmptrim_seg_${side}_warped_to_futrim_om.vtk

        # get thickness
        #$MATLAB_BIN -nojvm -nosplash -nodesktop <<-MATCODE
        #addpath('/home/longxie/ASHS_T1/Application/TAUPET/longitudinal/code');
        #MeasureMeanMedianThickness('$SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_thickmap.vtk','$SUBJALOHADIR/results_${sub}/blmptrim_seg_${side}_warped_to_futrim_om_thickmap.vtk','$SUBJALOHADIR/results_${sub}/mean_thickness_${side}.txt','$SUBJALOHADIR/results_${sub}/median_thickness_${side}.txt');
#MATCODE

     fi

    done
COMMENT
  done

  # clean up
  cp $SUBJALOHADIR/deformable/mprage_global_long_*_omRAS_half.mat $SUBJALOHADIR/
  cp $SUBJALOHADIR/deformable/mp_antsreg3d_*Warp*vec.nii.gz $SUBJALOHADIR/
  rm -rf $SUBJALOHADIR/final $SUBJALOHADIR/global \
         $SUBJALOHADIR/init \
         $SUBJALOHADIR/dump  $SUBJALOHADIR/bl_seg_*.nii.gz
}

######################################################
function AlohaVTHotspot()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  N=$(cat $OUTPUTDIR/checkdata.csv | wc -l)
  #N=1216
  N_begin=1100
  #N_begin=2

  # log
  echo "" >> $LOGDIR/ALOHAVTHotspot_log.txt
  echo $(date) >> $LOGDIR/ALOHAVTHotspot_log.txt

  # columns   RIDCol=1
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

  # Submit job to copy data
  JPREFIX=ALOHA${expid}
  for ((i=${N_begin};i<=${N};i++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

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
    blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
    blscandate=$(ReFormateDate $blscandate)
    Type=$(echo $ROW | cut -f $TYPECol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    #SUBJDIR=$DATADIR/$id

    # skip if it is not longitudinal
    echo "ALOHA Processing: $i,$rid,$PREFIX $Type"
    if [[ $Type != "T1Longitudinal" ]]; then
      continue
    fi

    # perform SR
    FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz
    FUASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
    FUASHST1LSEG=$FUASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
    FUASHST1RSEG=$FUASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

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
    BLMSTTHKDIR=$BLSUBJDIR/ASHST1_MTLCORTEX_MSTTHK
    BLLHOTSPOT=$BLMSTTHKDIR/VThotspot06032021/${id}_${BLPREFIX}_left_hotspot.nii.gz
    BLRHOTSPOT=$BLMSTTHKDIR/VThotspot06032021/${id}_${BLPREFIX}_right_hotspot.nii.gz

    # start running ALOHA  forward
    SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaVTHotspot06032021_ASHST1_${BLPREFIX}

    rm -rf $SUBJALOHADIR

    if [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $BLLHOTSPOT && -f $BLRHOTSPOT && ! -f $SUBJALOHADIR/results_label4/volumes_right.txt ]]; then

      echo "   Adding $rid $PREFIX $Type"
      echo "ALOHA hotspot: Adding $rid $PREFIX $Type" >> $LOGDIR/ALOHAVTHotspot_log.txt

      pybatch.sh \
        -m "5G" \
        -n 1 \
        -N "${JPREFIX}_${i}_${rid}_${PREFIX}" \
        -o $DUMPDIR \
        $0 AlohaHotspot_sub \
               $BLT1SRIMG \
               $FUT1SRIMG \
               $BLLHOTSPOT \
               $BLRHOTSPOT \
               $SUBJALOHADIR
      sleep 0.1

      #qsub -cwd -o $DUMPDIR -j y \
      #     -q all.q \
      #     -l h_vmem=5.1G,s_vmem=5G \
      #     -N "${JPREFIX}_${i}_${rid}_${PREFIX}" \
      #     $0 AlohaHotspot_sub \
      #         $BLT1SRIMG \
      #         $FUT1SRIMG \
      #         $BLLHOTSPOT \
      #         $BLRHOTSPOT \
      #         $SUBJALOHADIR
        sleep 0.1


      #fi
    fi

  done
 
  # Wait for completion
  pybatch.sh \
    -n 1 \
    -w "${JPREFIX}_*"
}


######################################################
function Summarize()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  QC_TXT=/home/lxie/ADNI2018/QC/QC_ALOHA20190613.csv
  #QC_TXT=$OUTPUTDIR/QC_ALOHA20190613.csv
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=10

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp

  # Submit job to copy data
  if [[ 1 == 1 ]]; then
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  PREFIX=SM
  for ((i=2;i<=${N};i++)); do

    echo "$i submitted"
    pybatch.sh \
      -m "4G" \
      -n 1 \
      -N "${PREFIX}_${i}" \
      -o $DUMPDIR \
      $0 Summarize_sub $i $OUTTMPDIR
    sleep 0.1

  done

  # Wait for completion
  ALLJOBS=$(bjobs | wc -l)
  while [ $ALLJOBS -gt 0 ]; do
    sleep 60
    ALLJOBS=$(bjobs | wc -l)
  done
  #pybatch.sh \
  #  -n 1 \
  #  -w "${PREFIX}_*"
  fi

  # get all info
  # header
  blheader=$(cat -A $DEMOG_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  #sed -e 's/\^M\$//g'
  header="$header,ALOHA_Success,DateDiffFromBLScanDate,EXCL_L,EXCL_R,EXCL_MRI"
  longheader=""
  for type in VOL MeanTHK MedianTHK; do
  for side in L R; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  done
  done
  done
  for type in MeanTHKFITTED MedianTHKFITTED; do
  for side in L R; do
  for sub in ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  done
  done
  done
  header="$header$longheader,$blheader,END"
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
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

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
  QC_TXT=/home/lxie/ADNI2018/QC/QC_ALOHA20190613.csv
  QCROW=$(cat -A $QC_TXT | grep $id | grep $scandate | grep T1Longitudinal | sed -e 's/\^M//g' | sed -e 's/\$//g')
  if [[ $QCROW == "" ]]; then
    excl_l=""
    excl_r=""
    excl_mri=""
  else
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

  # get thickness measurements from the fitted mesh
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


  echo "$ROW,$EXIST,$date_diff,$QCINFO$RAWMEASURE$LONGIMEASURE,$BLROW,-1" >> $OUTTMPDIR/$(printf %04d $i)_${id}.csv

  set -e

}

######################################################
function SummarizeHotspot()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  QC_TXT=/home/lxie/ADNI2018/output_runall_10012019/QC_ALOHA20190613.csv
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=1700
  N_begin=2

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp
  #rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # Submit job to copy data
  if [[ 1 == 0 ]]; then
  PREFIX=SM
  for ((i=${N_begin};i<=${N};i++)); do

    echo "$i submitted"
    pybatch.sh \
      -m "4G" \
      -n 1 \
      -N "${PREFIX}_${i}" \
      -o $DUMPDIR \
      $0 SummarizeHotspot_sub $i $OUTTMPDIR
    sleep 0.1


  done

  # Wait for completion
  pybatch.sh \
    -n 1 \
    -w "${JPREFIX}_*"
  fi

  # get all info
  # header
  blheader=$(cat -A $DEMOG_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  #sed -e 's/\^M\$//g'
  header="$header,ALOHA_Success,DateDiffFromBLScanDate,EXCL_L,EXCL_R,EXCL_MRI"
  longheader=""
  for type in VOL MeanTHK MedianTHK; do
  for side in L R; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  done
  done
  done
  for type in MeanTHKFITTED MedianTHKFITTED; do
  for side in L R; do
  for sub in ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  done
  done
  done
  header="$header$longheader,ALOHAHotspot_Success"
  longheader=""
  for type in VOL MeanTHK MedianTHK; do
  for side in L R; do
  for sub in HotspotL2 HotspotL3 HotspotL4; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  done
  done
  done
  header="$header$longheader,$blheader,END"
  echo $header > $OUTPUTDIR/longitudinal_hotspot_06032021.csv

  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/longitudinal_hotspot_06032021.csv
}

function SummarizeHotspot_sub()
{
  i=$1
  OUTTMPDIR=$2
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

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
  #QC_TXT=$OUTPUTDIR/QC_ALOHA20190613.csv
  QC_TXT=/home/lxie/ADNI2018/output_runall_10012019/QC_ALOHA20190613.csv
  QCROW=$(cat -A $QC_TXT | grep $id | grep $scandate | grep T1Longitudinal | sed -e 's/\^M//g' | sed -e 's/\$//g')
  if [[ $QCROW == "" ]]; then
    excl_l=""
    excl_r=""
    excl_mri=""
  else
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

  # get thickness measurements from the fitted mesh
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

  ######################## extract hotspot info ########################
  # check if exist
  BLMSTTHKDIR=$BLSUBJDIR/ASHST1_MTLCORTEX_MSTTHK
  BLLHOTSPOT=$BLMSTTHKDIR/VThotspot06032021/${id}_${BLPREFIX}_left_hotspot.nii.gz
  BLRHOTSPOT=$BLMSTTHKDIR/VThotspot06032021/${id}_${BLPREFIX}_right_hotspot.nii.gz
  EXISTHSP=0
  if [[ $Type != "T1Longitudinal" ]]; then
    EXISTHSP="-1"
  elif [[ -f $FUT1SRIMG && -f $BLT1SRIMG && -f $BLLHOTSPOT && -f $BLRHOTSPOT ]]; then
    EXISTHSP=1
  else
    EXISTHSP=0
  fi
 

  SUBJALOHAHSPDIR=$SUBJALLDIR/${PREFIX}/AlohaVTHotspot06032021_ASHST1_${BLPREFIX}
  # go through all the fields
  HSPRAWMEASURE=""
  HSPLONGIMEASURE=""
  for type in volumes; do
  for side in left right; do
    for sub in label2 label3 label4; do
    SUBJALOHASUBDIR=$SUBJALOHAHSPDIR/results_${sub}
      if [[ $EXISTHSP == 1 ]]; then
        # get volume and thickness
        if [[ -f $SUBJALOHASUBDIR/${type}_${side}.txt ]]; then
          MEA=$(cat $SUBJALOHASUBDIR/${type}_${side}.txt)
          MEA="${MEA// /,}"
          BL=$(echo $MEA | cut -d , -f 1)
          FU=$(echo $MEA | cut -d , -f 2)
          CHANGE=$(echo "scale=7;($BL-$FU)*100*365/$BL/$date_diff" | bc)
          HSPLONGIMEASURE="$HSPLONGIMEASURE,$CHANGE"
        else
          MEA=","
          EXISTHSP=0
          HSPLONGIMEASURE="$HSPLONGIMEASURE,"
        fi
        HSPRAWMEASURE="$HSPRAWMEASURE,$MEA"
      else
        HSPRAWMEASURE="$HSPRAWMEASURE,,"
        HSPLONGIMEASURE="$HSPLONGIMEASURE,"
      fi
    done
  done
  done


  echo "$ROW,$EXIST,$date_diff,$QCINFO$RAWMEASURE$LONGIMEASURE,$EXISTHSP$HSPRAWMEASURE$HSPLONGIMEASURE,$BLROW,-1" >> $OUTTMPDIR/$(printf %04d $i)_${id}.csv

  set -e

}


######################################################
function SummarizeRobin()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  #ROBIN_TXT=/data/jux/rdeflores/scripts/thickness_ADNI/longit/list/list_longit_ID_Date.txt
  ROBIN_TXT=$ANALYSISDIR/list_longit_ID_Date_Robin.txt
  QC_TXT=$OUTPUTDIR/QC_ALOHA20190613.csv
  N=$(cat $ROBIN_TXT | wc -l)
  #N=10

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # Submit job to copy data
  if [[ 1 == 1 ]]; then
  PREFIX=SMR
  for ((i=1;i<=${N};i++)); do

    echo "$i submitted"
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX}_${i}" \
         $0 SummarizeRobin_sub $i $OUTTMPDIR
    sleep 0.05

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  fi

  # get all info
  # header
  blheader=$(cat -A $DEMOG_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
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
  #ROBIN_TXT=/data/jux/rdeflores/scripts/thickness_ADNI/longit/list/list_longit_ID_Date.txt
  ROBIN_TXT=$ANALYSISDIR/list_longit_ID_Date_Robin.txt
  QC_TXT=$OUTPUTDIR/QC_ALOHA20190613.csv

  # get row from Robin spreadsheet
  RBROW=$(cat -A $ROBIN_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g') 
  id=$(echo $RBROW | cut -f 1 -d ",")
  scandate=$(echo $RBROW | cut -f 3 -d ",")
  blscandate=$(echo $RBROW | cut -f 2 -d ",")

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)
  EXCL_L_Col=$(csvcol.sh $QC_TXT EXCL_L)
  EXCL_R_Col=$(csvcol.sh $QC_TXT EXCL_R)
  EXCL_MRI_Col=$(csvcol.sh $QC_TXT EXCL_MRI)

  # search allinfo for this line
  ROW=$(cat -A $ALLINFO_TXT | grep $id | grep $scandate | grep $blscandate | while read line; do
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
  QCROW=$(cat -A $QC_TXT | grep $id | grep $scandate | grep T1Longitudinal | sed -e 's/\^M//g' | sed -e 's/\$//g')
  if [[ $QCROW == "" ]]; then
    excl_l=""
    excl_r=""
    excl_mri=""
  else
    EXCL_L_Col=$(csvcol.sh $QC_TXT EXCL_L)
    EXCL_R_Col=$(csvcol.sh $QC_TXT EXCL_R)
    EXCL_MRI_Col=$(csvcol.sh $QC_TXT EXCL_MRI)
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
function SummarizeMengjin()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  QC_TXT=$OUTPUTDIR/QC_ALOHA20190613.csv
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=2030

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # Submit job to copy data
  if [[ 1 == 1 ]]; then
  PREFIX1=SM
  for ((i=2;i<=${N};i++)); do
    for type in forward reverse; do

      echo "$i $type submitted"
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX1}_${i}_${type}" \
           $0 SummarizeMengjin_sub $i $type $OUTTMPDIR
      sleep 0.05

    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX1}_*" -sync y -b y \
       sleep 1

  fi

  # get all info
  # header
  blheader=$(cat -A $DEMOG_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  #sed -e 's/\^M\$//g'
  header="$header,ALOHA_Success,ALOHA_BLScanDate,ALOHA_FUScanDate,ALOHA_Dir,DateDiffFromBLScanDate,ALOHA_run_type,EXCL_L,EXCL_R,EXCL_MRI"
  longheader=""
  for type in VOL MeanTHK MedianTHK; do
  for side in L R; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  done
  done
  done
  header="$header$longheader,$blheader,END"
  echo $header > $OUTPUTDIR/longitudinal_Mengjin.csv

  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/longitudinal_Mengjin.csv

}

function SummarizeMengjin_sub()
{
  i=$1
  type=$2
  OUTTMPDIR=$3
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  ROW=$(cat -A $ALLINFO_TXT | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

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
  FUASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
  FUASHST1LSEG=$FUASHST1DIR/final/${id}_left_lfseg_heur.nii.gz
  FUASHST1RSEG=$FUASHST1DIR/final/${id}_right_lfseg_heur.nii.gz

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

  # output status
  echo "Checking ASHS ICV: $i,$rid,$PREFIX"
  set +e

  # baseline row
  BLROW=$(cat -A $DEMOG_TXT | grep $id | sed -e 's/\^M//g' | sed -e 's/\$//g')
  if [[ $BLROW == "" ]]; then
    BLROW=$(cat -A $DEMOG_TXT | head -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | awk  -F "," '{for(i=1;i<=NF-1;i++)printf ","}')
  fi

  # aloha baseline and followup
  if [[ $type == "forward" ]]; then
    ALOHABLSCANDATE=$blscandate
    ALOHAFUSCANDATE=$scandate
    ALOHABLT1SRIMG=$BLT1SRIMG
    ALOHAFUT1SRIMG=$FUT1SRIMG
    ALOHABLASHST1LSEG=$BLASHST1LSEG
    ALOHABLASHST1RSEG=$BLASHST1RSEG
    SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
    POSTFIX=0
  else
    ALOHABLSCANDATE=$scandate
    ALOHAFUSCANDATE=$blscandate
    ALOHABLT1SRIMG=$FUT1SRIMG
    ALOHAFUT1SRIMG=$BLT1SRIMG
    ALOHABLASHST1LSEG=$FUASHST1LSEG
    ALOHABLASHST1RSEG=$FUASHST1RSEG
    SUBJALOHADIR=$SUBJALLDIR/${BLPREFIX}/AlohaMTL_ASHST1_${PREFIX}
    POSTFIX=1
  fi

  # check if exist
  EXIST=0
  if [[ $Type != "T1Longitudinal" ]]; then
    EXIST="-1"
    SUBJALOHADIR=""
  elif [[ -f $ALOHAFUT1SRIMG && -f $ALOHABLT1SRIMG && -f $ALOHABLASHST1LSEG && -f $ALOHABLASHST1RSEG ]]; then
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
  if [[ $type == "forward" ]]; then
    QC_TXT=$OUTPUTDIR/QC_ALOHA20190613.csv
  else
    QC_TXT=$OUTPUTDIR/QC_ALOHA20190402_reverse.csv
  fi
  QCROW=$(cat -A $QC_TXT | grep $id | grep $scandate | grep T1Longitudinal | sed -e 's/\^M//g' | sed -e 's/\$//g')
  if [[ $QCROW == "" ]]; then
    excl_l=""
    excl_r=""
    excl_mri=""
  else
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

  echo "$ROW,$EXIST,$ALOHABLSCANDATE,$ALOHAFUSCANDATE,$SUBJALOHADIR,$date_diff,$type,$QCINFO$RAWMEASURE$LONGIMEASURE,$BLROW,-1" >> $OUTTMPDIR/$(printf %04d $i)_${id}_${POSTFIX}.csv

  set -e

}

######################################################
function SummarizeComparisons()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  QC_TXT=$OUTPUTDIR/QC_ALOHA20190613.csv
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=1100
  #N_begin=1054
  N_begin=2

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # Submit job to copy data
  if [[ 1 == 1 ]]; then
  PREFIX=SM
  for ((i=$N_begin;i<=${N};i++)); do

    echo "$i submitted"
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,gpu.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX}_${i}" \
         $0 SummarizeComparisons_sub $i $OUTTMPDIR
    sleep 0.05

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,gpu.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 10

  fi

  # get all info
  # header
  #blheader=$(cat -A $DEMOG_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  #sed -e 's/\^M\$//g'
  header="$header,ALOHA_Success,DateDiffFromBLScanDate,EXCL_L,EXCL_R,EXCL_MRI"
  longheader=""
  for type in VOL MeanTHK MedianTHK; do
  for side in L R; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  done
  done
  done
  for type in MeanTHKFITTED MedianTHKFITTED; do
  for side in L R; do
  for sub in ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_BL,${side}_${sub}_${type}_FU"
    longheader="$longheader,${side}_${sub}_${type}_ChangeAnnualized"
  done
  done
  done
  header="$header$longheader"

  for method in MSTVT MSTUT OLDMTVT OLDMTUT; do
  for side in L R M; do 
  for type in MeanTHK MedianTHK; do
  for sub in ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_${method}"
  done
  done
  done
  done

  for type in VOL MeanTHK MedianTHK; do
  for side in L R M; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_ALOHAMeasures"
  done
  done
  done

  for side in L R M; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC CS OTS AMEN PMEN MISC; do
    header="$header,${side}_${sub}_VOL_ASHST1RAW"
  done
  done

  # output
  header="$header,COMEND"
  echo $header > $OUTPUTDIR/longitudinal_comparisons.csv

  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/longitudinal_comparisons.csv

}

function SummarizeComparisons_sub()
{
  ii=$1
  OUTTMPDIR=$2
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

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

  #######################################################
  # get longitudinal measurements
  SUBJALOHADIR=$SUBJALLDIR/${PREFIX}/AlohaMTL_ASHST1_${BLPREFIX}
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

  # get QC information
  QC_TXT=$OUTPUTDIR/QC_ALOHA20190613.csv
  QCROW=$(cat -A $QC_TXT | grep $id | grep $scandate | grep T1Longitudinal | sed -e 's/\^M//g' | sed -e 's/\$//g')
  if [[ $QCROW == "" ]]; then
    excl_l=""
    excl_r=""
    excl_mri=""
  else
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

  # get thickness measurements from the fitted mesh
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

  #####################################################
  # get MST variant template thickness
  SUBJMSTTHKDIR=$SUBJALLDIR/$PREFIX/ASHST1_MTLCORTEX_MSTTHK
  MSTVTTHK=""
  for side in left right; do
    idside=${id}_${side}
    SUBJMSTVTGEOSHOOTDIR=$SUBJMSTTHKDIR/GeoShoot/$side
    MeanTHK=$(cat -A $SUBJMSTVTGEOSHOOTDIR/${idside}_mean_thickness.vtk | cut -d , -f 3-6 | sed -e 's/\$//g')
    if [[ $MeanTHK == "" ]]; then
      MeanTHK=",,,"
    fi
    MedianTHK=$(cat -A $SUBJMSTVTGEOSHOOTDIR/${idside}_median_thickness.vtk | cut -d , -f 3-6 | sed -e 's/\$//g')
    if [[ $MedianTHK == "" ]]; then
      MedianTHK=",,,"
    fi
    MSTVTTHK="$MSTVTTHK,$MeanTHK,$MedianTHK"
  done
  MSTVTTHK="${MSTVTTHK:1}"
  # compute mean
  for ((i=1;i<=8;i++)); do
    LMEA=$(echo $MSTVTTHK | cut -d, -f $i)
    RMEA=$(echo $MSTVTTHK | cut -d, -f $((i+8)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=10;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    MSTVTTHK="$MSTVTTHK,$MMEA"
  done

  #####################################################
  # get MST unified template thicvkness
  MSTUTTHK=""
  for side in left right; do
    idside=${id}_${side}
    SUBJMSTUTGEOSHOOTDIR=$SUBJMSTTHKDIR/GeoShootUT/$side
    MeanTHK=$(cat -A $SUBJMSTUTGEOSHOOTDIR/${idside}_mean_thickness.vtk | cut -d , -f 3-6 | sed -e 's/\$//g')
    if [[ $MeanTHK == "" ]]; then
      MeanTHK=",,,"
    fi
    MedianTHK=$(cat -A $SUBJMSTUTGEOSHOOTDIR/${idside}_median_thickness.vtk | cut -d , -f 3-6 | sed -e 's/\$//g')
    if [[ $MedianTHK == "" ]]; then
      MedianTHK=",,,"
    fi
    MSTUTTHK="$MSTUTTHK,$MeanTHK,$MedianTHK"
  done
  MSTUTTHK="${MSTUTTHK:1}"
  # compute mean
  for ((i=1;i<=8;i++)); do
    LMEA=$(echo $MSTUTTHK | cut -d, -f $i)
    RMEA=$(echo $MSTUTTHK | cut -d, -f $((i+8)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=10;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    MSTUTTHK="$MSTUTTHK,$MMEA"
  done

  #####################################################
  # get multi-template variant template thickness (the old way)
  SUBJOLDMTTHKDIR=$SUBJALLDIR/$PREFIX/ASHST1_MTLCORTEX_MultiTempThk
  OLDMTVTTHK=""
  for side in left right; do
    SUBJOLDMTVTTHKDIR=$SUBJOLDMTTHKDIR/RegToMultTemp
    MeanTHK=$(cat -A $SUBJOLDMTVTTHKDIR/${idside}_mean_thickness.txt | cut -d , -f 3-6 | sed -e 's/\$//g')
    if [[ $MeanTHK == "" ]]; then
      MeanTHK=",,,"
    fi
    MedianTHK=$(cat -A $SUBJOLDMTVTTHKDIR/${idside}_median_thickness.txt | cut -d , -f 3-6 | sed -e 's/\$//g')
    if [[ $MedianTHK == "" ]]; then
      MedianTHK=",,,"
    fi
    OLDMTVTTHK="$OLDMTVTTHK,$MeanTHK,$MedianTHK"
  done
  OLDMTVTTHK="${OLDMTVTTHK:1}"
  # compute mean
  for ((i=1;i<=8;i++)); do
    LMEA=$(echo $OLDMTVTTHK | cut -d, -f $i)
    RMEA=$(echo $OLDMTVTTHK | cut -d, -f $((i+8)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=10;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    OLDMTVTTHK="$OLDMTVTTHK,$MMEA"
  done

  #####################################################
  # get multi-template unified template thickness (the old way)
  OLDMTUTTHK=""
  for side in left right; do
    SUBJOLDMTUTTHKDIR=$SUBJOLDMTTHKDIR/RegToUniTemp
    MeanTHK=$(cat -A $SUBJOLDMTUTTHKDIR/${idside}_mean_thickness.txt | cut -d , -f 3-6 | sed -e 's/\$//g')
    if [[ $MeanTHK == "" ]]; then
      MeanTHK=",,,"
    fi
    MedianTHK=$(cat -A $SUBJOLDMTUTTHKDIR/${idside}_median_thickness.txt | cut -d , -f 3-6 | sed -e 's/\$//g')
    if [[ $MedianTHK == "" ]]; then
      MedianTHK=",,,"
    fi
    OLDMTUTTHK="$OLDMTUTTHK,$MeanTHK,$MedianTHK"
  done
  OLDMTUTTHK="${OLDMTUTTHK:1}"
  # compute mean
  for ((i=1;i<=8;i++)); do
    LMEA=$(echo $OLDMTUTTHK | cut -d, -f $i)
    RMEA=$(echo $OLDMTUTTHK | cut -d, -f $((i+8)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=10;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    OLDMTUTTHK="$OLDMTUTTHK,$MMEA"
  done

  #####################################################
  # get ALOHA thickness and volume measurement (baseline, look for the folders in followup) 
  ALOHAMEAS=""
  if [[ $Type == "T1Baseline" ]]; then

    # find the cloest longitudinal timepoint row
    ROWLONG=$(cat -A $ALLINFO_TXT | grep $id | grep T1Longitudinal | head -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
    longscandate=$(echo $ROWLONG | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$longscandate" '+%d')
    MM=$(date -d "$longscandate" '+%m')
    YYYY=$(date -d "$longscandate" '+%Y')
    LONGPREFIX="${YYYY}-${MM}-${DD}"
    longscandate=$(ReFormateDate $longscandate)
    longType=$(echo $ROW | cut -f $TYPECol -d ",")
    LONGALOHADIR=$SUBJALLDIR/${LONGPREFIX}/AlohaMTL_ASHST1_${BLPREFIX}

    # extract aloha info
    for type in volumes mean_thickness median_thickness; do
    TMPMEAS=""
    for side in left right; do
      for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
        LONGALOHASUBDIR=$LONGALOHADIR/results_${sub}
        # get volume and thickness
        if [[ -f $LONGALOHASUBDIR/${type}_${side}.txt ]]; then
          MEA=$(cat $LONGALOHASUBDIR/${type}_${side}.txt)
          MEA="${MEA// /,}"
          BL=$(echo $MEA | cut -d , -f 1)
        else
          BL=""
        fi
        TMPMEAS="$TMPMEAS,$BL"
      done
    done
    TMPMEAS="${TMPMEAS:1}"
    # compute mean  
    for ((i=1;i<=7;i++)); do
      LMEA=$(echo $TMPMEAS | cut -d, -f $i)
      RMEA=$(echo $TMPMEAS | cut -d, -f $((i+7)))
      if [[ $LMEA != "" && $RMEA != "" ]]; then
        MMEA=$(echo "scale=10;($LMEA+$RMEA)/2" | bc -l)
      else
        MMEA=""
      fi
      TMPMEAS="$TMPMEAS,$MMEA"
    done

    # output all
    ALOHAMEAS="$ALOHAMEAS,$TMPMEAS"
    done

  elif [[ $Type == "T1Longitudinal" && $EXIST == 1 ]]; then

    for type in volumes mean_thickness median_thickness; do
    TMPMEAS=""
    for side in left right; do
      for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
        SUBJALOHASUBDIR=$SUBJALOHADIR/results_${sub}
        # get volume and thickness
        if [[ -f $SUBJALOHASUBDIR/${type}_${side}.txt ]]; then
          MEA=$(cat $SUBJALOHASUBDIR/${type}_${side}.txt)
          MEA="${MEA// /,}"
          FU=$(echo $MEA | cut -d , -f 2)
        else
          FU=""
        fi
        TMPMEAS="$TMPMEAS,$FU"
      done
    done
    TMPMEAS="${TMPMEAS:1}"
    # compute mean
    for ((i=1;i<=7;i++)); do
      LMEA=$(echo $TMPMEAS | cut -d, -f $i)
      RMEA=$(echo $TMPMEAS | cut -d, -f $((i+7)))
      if [[ $LMEA != "" && $RMEA != "" ]]; then
        MMEA=$(echo "scale=10;($LMEA+$RMEA)/2" | bc -l)
      else
        MMEA=""
      fi
      TMPMEAS="$TMPMEAS,$MMEA"
    done

    # output all
    ALOHAMEAS="$ALOHAMEAS,$TMPMEAS"
    done

  else
    # output empty fields
    for type in volumes mean_thickness median_thickness; do
    for side in left right mean; do
      for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
        ALOHAMEAS="$ALOHAMEAS,"
      done
    done
    done
  fi

  #####################################################
  # get volume measurements
  SUBJASHSDIR=$SUBJALLDIR/$PREFIX/ASHST1
  VOL=""
  for side in left right; do
    idside=${id}_${side}
    VOL_AHIPPO=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Anterior_hippocampus | cut -d ' ' -f 5)
    VOL_PHIPPO=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Posterior_hippocampus | cut -d ' ' -f 5)
    VOL_HIPPO=$(echo "scale=3;$VOL_AHIPPO+$VOL_PHIPPO" | bc -l)
    VOL_PHC=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep " PHC" | cut -d ' ' -f 5 )
    VOL_MISC=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep MISC | cut -d ' ' -f 5)
    VOL_ERC=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep ERC | cut -d ' ' -f 5)
    VOL_BA35=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Br35 | cut -d ' ' -f 5)
    VOL_BA36=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Br36 | cut -d ' ' -f 5)
    VOL_CS=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep ColSul | cut -d ' ' -f 5)
    VOL_OTS=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep OTSul | cut -d ' ' -f 5)
    VOL_AMEN=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep "Meninges " | cut -d ' ' -f 5)
    VOL_PMEN=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep "Meninges_PHC" | cut -d ' ' -f 5)
    if [[ $VOL == "" ]]; then
      VOL="$VOL_AHIPPO,$VOL_PHIPPO,$VOL_HIPPO,$VOL_ERC,$VOL_BA35,$VOL_BA36,$VOL_PHC,$VOL_CS,$VOL_OTS,$VOL_AMEN,$VOL_PMEN,$VOL_MISC"
    else
      VOL="$VOL,$VOL_AHIPPO,$VOL_PHIPPO,$VOL_HIPPO,$VOL_ERC,$VOL_BA35,$VOL_BA36,$VOL_PHC,$VOL_CS,$VOL_OTS,$VOL_AMEN,$VOL_PMEN,$VOL_MISC"
    fi
  done
  # compute mean
  for ((i=1;i<=12;i++)); do
    LMEA=$(echo $VOL | cut -d, -f $i)
    RMEA=$(echo $VOL | cut -d, -f $((i+12)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=3;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    VOL="$VOL,$MMEA"
  done
  

  echo "$ROW,$EXIST,$date_diff,$QCINFO$RAWMEASURE$LONGIMEASURE,$MSTVTTHK,$MSTUTTHK,$OLDMTVTTHK,$OLDMTUTTHK$ALOHAMEAS,$VOL,-1" >> $OUTTMPDIR/$(printf %04d $ii)_${id}.csv

  set -e

}

######################################################
function SummarizeMismatchCohort()
{
  # Load id number
  #ALLINFO_TXT=$ANALYSISDIR/Mismatch_MRIDates.csv
  ALLINFO_TXT=$ANALYSISDIR/Sandy_aplus_04012021.csv
  #ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  #QC_TXT=$OUTPUTDIR/QC_ALOHA20190613.csv
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=1100
  #N_begin=1054
  N_begin=2

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # Submit job to copy data
  if [[ 1 == 1 ]]; then
  PREFIX=SM
  for ((i=$N_begin;i<=${N};i++)); do

    echo "$i submitted"
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX}_${i}" \
         $0 SummarizeMismatchCohort_sub $i $OUTTMPDIR $ALLINFO_TXT
    sleep 0.05

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,gpu.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 10

  fi

  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header="$header,ICV_ASHST1"

  for type in MeanTHK MedianTHK; do
  for side in L R M; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_ALOHAMeasures"
  done
  done
  done

  for side in L R M; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC CS OTS AMEN PMEN MISC; do
    header="$header,${side}_${sub}_VOL_ASHST1RAW"
  done
  done

  #for method in MSTVT MSTUT; do
  for method in MSTVT; do
  for side in L R M; do
  for type in MeanTHK MedianTHK; do
  for sub in ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_${method}"
  done
  done
  done
  done

  # output
  header="$header,COMEND"
  echo $header > $OUTPUTDIR/Mismatch_MRIDates_ASHST1_volume_thickness.csv
  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/Mismatch_MRIDates_ASHST1_volume_thickness.csv
}

function SummarizeMismatchCohort_sub()
{
  ii=$1
  OUTTMPDIR=$2
  ALLINFO_TXT=$3
  #ALLINFO_TXT=$ANALYSISDIR/Mismatch_MRIDates.csv
  #ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT MRISCANDATE)

  # get information
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  SUBJALLDIR=$ALLDATADIR/$id
  FUT1SRIMG=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim_denoised_SR.nii.gz

  set +e

  # get ICV from ASHS T1 pipeline
  SUBJASHSDIR=$SUBJALLDIR/$PREFIX/ASHST1
  ICV=$(cat $SUBJASHSDIR/final/${id}_icv.txt | cut -d " " -f 2)

  #####################################################
  # get ALOHA thickness and volume measurement (baseline, look for the folders in followup)
  ALOHAMEAS=""
  LONGALOHADIRS=($(ls $SUBJALLDIR/$PREFIX | grep AlohaMTL_ASHST1))
  LONGALOHADIR=$SUBJALLDIR/$PREFIX/${LONGALOHADIRS[0]}
  # extract aloha info
  for type in mean_thickness median_thickness; do
  TMPMEAS=""
  for side in left right; do
    for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC; do
      LONGALOHASUBDIR=$LONGALOHADIR/results_${sub}
      # get volume and thickness
      if [[ -f $LONGALOHASUBDIR/${type}_${side}.txt ]]; then
        MEA=$(cat $LONGALOHASUBDIR/${type}_${side}.txt)
        MEA="${MEA// /,}"
        BL=$(echo $MEA | cut -d , -f 1)
      else
        BL=""
      fi
      TMPMEAS="$TMPMEAS,$BL"
    done
  done
  TMPMEAS="${TMPMEAS:1}"
  # compute mean
  for ((i=1;i<=7;i++)); do
    LMEA=$(echo $TMPMEAS | cut -d, -f $i)
    RMEA=$(echo $TMPMEAS | cut -d, -f $((i+7)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=10;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    TMPMEAS="$TMPMEAS,$MMEA"
  done

  # output all
  ALOHAMEAS="$ALOHAMEAS,$TMPMEAS"
  done

  #####################################################
  # get volume measurements
  SUBJASHSDIR=$SUBJALLDIR/$PREFIX/ASHST1
  VOL=""
  for side in left right; do
    idside=${id}_${side}
    VOL_AHIPPO=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Anterior_hippocampus | cut -d ' ' -f 5)
    VOL_PHIPPO=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Posterior_hippocampus | cut -d ' ' -f 5)
    VOL_HIPPO=$(echo "scale=3;$VOL_AHIPPO+$VOL_PHIPPO" | bc -l)
    VOL_PHC=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep " PHC" | cut -d ' ' -f 5 )
    VOL_MISC=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep MISC | cut -d ' ' -f 5)
    VOL_ERC=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep ERC | cut -d ' ' -f 5)
    VOL_BA35=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Br35 | cut -d ' ' -f 5)
    VOL_BA36=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Br36 | cut -d ' ' -f 5)
    VOL_CS=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep ColSul | cut -d ' ' -f 5)
    VOL_OTS=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep OTSul | cut -d ' ' -f 5)
    VOL_AMEN=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep "Meninges " | cut -d ' ' -f 5)
    VOL_PMEN=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep "Meninges_PHC" | cut -d ' ' -f 5)
    if [[ $VOL == "" ]]; then
      VOL="$VOL_AHIPPO,$VOL_PHIPPO,$VOL_HIPPO,$VOL_ERC,$VOL_BA35,$VOL_BA36,$VOL_PHC,$VOL_CS,$VOL_OTS,$VOL_AMEN,$VOL_PMEN,$VOL_MISC"
    else
      VOL="$VOL,$VOL_AHIPPO,$VOL_PHIPPO,$VOL_HIPPO,$VOL_ERC,$VOL_BA35,$VOL_BA36,$VOL_PHC,$VOL_CS,$VOL_OTS,$VOL_AMEN,$VOL_PMEN,$VOL_MISC"
    fi
  done
  # compute mean
  for ((i=1;i<=12;i++)); do
    LMEA=$(echo $VOL | cut -d, -f $i)
    RMEA=$(echo $VOL | cut -d, -f $((i+12)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=3;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    VOL="$VOL,$MMEA"
  done


  #####################################################
  # get MST variant template thickness
  SUBJMSTTHKDIR=$SUBJALLDIR/$PREFIX/ASHST1_MTLCORTEX_MSTTHK
  MSTVTTHK=""
  for side in left right; do
    MEANTHKFN=$SUBJMSTTHKDIR/${id}_${PREFIX}_${side}_thickness.csv
    SUBJMSTVTGEOSHOOTDIR=$SUBJMSTTHKDIR/GeoShoot/$side
    MeanTHK=""
    if [[ -f $MEANTHKFN ]]; then
      MeanTHK=$(cat -A $MEANTHKFN | head -n 2 | tail -n 1 | cut -d , -f 5-8 | sed -e 's/\$//g')
    elif [[ -f $SUBJMSTVTGEOSHOOTDIR/${id}_${side}_mean_thickness.vtk ]]; then
      MeanTHK=$(cat -A $SUBJMSTVTGEOSHOOTDIR/${id}_${side}_mean_thickness.vtk | cut -d , -f 3-6 | sed -e 's/\$//g')
    fi
    if [[ $MeanTHK == "" ]]; then
      MeanTHK=",,,"
    fi

    MEDIANTHKFN=$SUBJMSTTHKDIR/${id}_${PREFIX}_${side}_thickness.csv
    MedianTHK=""
    if [[ -f $MEDIANTHKFN ]]; then
      MedianTHK=$(cat -A $MEDIANTHKFN | head -n 2 | tail -n 1 | cut -d , -f 5-8 | sed -e 's/\$//g')
    elif [[ -f $SUBJMSTVTGEOSHOOTDIR/${id}_${side}_median_thickness.vtk ]]; then
      MedianTHK=$(cat -A $SUBJMSTVTGEOSHOOTDIR/${id}_${side}_median_thickness.vtk | cut -d , -f 3-6 | sed -e 's/\$//g')
    fi
    if [[ $MedianTHK == "" ]]; then
      MedianTHK=",,,"
    fi

    MSTVTTHK="$MSTVTTHK,$MeanTHK,$MedianTHK"
  done
  MSTVTTHK="${MSTVTTHK:1}"
  # compute mean
  for ((i=1;i<=8;i++)); do
    LMEA=$(echo $MSTVTTHK | cut -d, -f $i)
    RMEA=$(echo $MSTVTTHK | cut -d, -f $((i+8)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=10;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    MSTVTTHK="$MSTVTTHK,$MMEA"
  done


  echo "$ROW,$ICV$ALOHAMEAS,$VOL,$MSTVTTHK,-1" >> $OUTTMPDIR/$(printf %04d $ii)_${id}.csv

  set -e

}

######################################################
function SummarizeBLMengjin()
{
  # Load id number
  SUBJLIST_TXT=$ANALYSISDIR/subjects_2GO_Mengjin2020Paper.csv
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  BLINFO_TXT=$ANALYSISDIR/ADNI_T1_AllInfo_MTL_20181025.csv
  N=$(cat $SUBJLIST_TXT | wc -l)
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
           -q all.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX1}_${i}" \
           $0 SummarizeBLMengjin_sub $i $OUTTMPDIR
      sleep 0.05

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX1}_*" -sync y -b y \
       sleep 1

  fi

  # get all info
  # header
  blheader=$(cat -A $BLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header="$blheader,ICV_ASHS"

  for side in L R M; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC CS OTS AMEN PMEN MISC; do
    header="$header,${side}_${sub}_VOL_ASHST1RAW"
  done
  done

  #for method in MSTVT MSTUT; do
  for method in MSTVT; do
  for side in L R M; do
  for type in MeanTHK MedianTHK; do
  for sub in ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_${method}"
  done
  done
  done
  done

  #header="$header,QC_MRI,QC_ICV,QC_CORTEX_LEFT,QC_CORTEX_R,QC_CORTEX_BOTH,QC_HIPPO_L,QC_HIPPO_R,QC_HIPPO_BOTH,END"
  header="$header,END"
  echo $header > $OUTPUTDIR/baseline_Mengjin2020Paper.csv

  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/baseline_Mengjin2020Paper.csv
}

function SummarizeBLMengjin_sub()
{
  ii=$1
  OUTTMPDIR=$2
  SUBJLIST_TXT=$ANALYSISDIR/subjects_2GO_Mengjin2020Paper.csv
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  BLINFO_TXT=$ANALYSISDIR/ADNI_T1_AllInfo_MTL_20181025.csv

  # get id
  id=$(cat $SUBJLIST_TXT | head -n $ii | tail -n 1)
  
  # get qc info
  #QCMRICol=$(csvcol.sh $BLINFO_TXT QC_MRI)
  #QCICVCol=$(csvcol.sh $BLINFO_TXT QC_ICV)
  #QCCORTEXLCol=$(csvcol.sh $BLINFO_TXT QC_MTLCORTEX_LEFT)
  #QCCORTEXRCol=$(csvcol.sh $BLINFO_TXT QC_MTLCORTEX_RIGHT)
  #QCCORTEXBCol=$(csvcol.sh $BLINFO_TXT QC_MTLCORTEX_BOTH)
  #QCHIPPOLCol=$(csvcol.sh $BLINFO_TXT QC_HIPPO_LEFT)
  #QCHIPPORCol=$(csvcol.sh $BLINFO_TXT QC_HIPPO_RIGHT)
  #QCHIPPOBCol=$(csvcol.sh $BLINFO_TXT QC_HIPPO_BOTH)
  
  BLROW=$(cat -A $BLINFO_TXT | grep $id | sed -e 's/\^M//g' | sed -e 's/\$//g')
  #QCMRI=$(echo $BLROW | cut -f $QCMRICol -d,)
  #QCICV=$(echo $BLROW | cut -f $QCICVCol -d,)
  #QCCORTEXL=$(echo $BLROW | cut -f $QCCORTEXLCol -d,)
  #QCCORTEXR=$(echo $BLROW | cut -f $QCCORTEXBCol -d,)
  #QCCORTEXB=$(echo $BLROW | cut -f $QCCORTEXRCol -d,)
  #QCHIPPOL=$(echo $BLROW | cut -f $QCHIPPOLCol -d,)
  #QCHIPPOR=$(echo $BLROW | cut -f $QCHIPPORCol -d,)
  #QCHIPPOB=$(echo $BLROW | cut -f $QCHIPPOBCol -d,)

  # get baseline scandate
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  ROWS=($(cat -A $ALLINFO_TXT | grep $id))
  ROW=${ROWS[0]}
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

  # get ICV from ASHS T1 pipeline
  SUBJALLDIR=$ALLDATADIR/$id
  SUBJASHSICVDIR=$SUBJALLDIR/$BLPREFIX/ASHSICV
  ICV=$(cat $SUBJASHSICVDIR/final/${id}_left_corr_nogray_volumes.txt | cut -d " " -f 5)

  #####################################################
  # get volume measurements
  SUBJASHSDIR=$SUBJALLDIR/$BLPREFIX/ASHST1
  VOL=""
  for side in left right; do
    idside=${id}_${side}
    VOL_AHIPPO=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Anterior_hippocampus | cut -d ' ' -f 5)
    VOL_PHIPPO=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Posterior_hippocampus | cut -d ' ' -f 5)
    VOL_HIPPO=$(echo "scale=3;$VOL_AHIPPO+$VOL_PHIPPO" | bc -l)
    VOL_PHC=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep " PHC" | cut -d ' ' -f 5 )
    VOL_MISC=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep MISC | cut -d ' ' -f 5)
    VOL_ERC=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep ERC | cut -d ' ' -f 5)
    VOL_BA35=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Br35 | cut -d ' ' -f 5)
    VOL_BA36=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep Br36 | cut -d ' ' -f 5)
    VOL_CS=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep ColSul | cut -d ' ' -f 5)
    VOL_OTS=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep OTSul | cut -d ' ' -f 5)
    VOL_AMEN=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep "Meninges " | cut -d ' ' -f 5)
    VOL_PMEN=$(cat $SUBJASHSDIR/final/${idside}_heur_volumes.txt | grep "Meninges_PHC" | cut -d ' ' -f 5)
    if [[ $VOL == "" ]]; then
      VOL="$VOL_AHIPPO,$VOL_PHIPPO,$VOL_HIPPO,$VOL_ERC,$VOL_BA35,$VOL_BA36,$VOL_PHC,$VOL_CS,$VOL_OTS,$VOL_AMEN,$VOL_PMEN,$VOL_MISC"
    else
      VOL="$VOL,$VOL_AHIPPO,$VOL_PHIPPO,$VOL_HIPPO,$VOL_ERC,$VOL_BA35,$VOL_BA36,$VOL_PHC,$VOL_CS,$VOL_OTS,$VOL_AMEN,$VOL_PMEN,$VOL_MISC"
    fi
  done
  # compute mean
  for ((i=1;i<=12;i++)); do
    LMEA=$(echo $VOL | cut -d, -f $i)
    RMEA=$(echo $VOL | cut -d, -f $((i+12)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=3;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    VOL="$VOL,$MMEA"
  done
 
  #####################################################
  # get MST variant template thickness
  SUBJMSTTHKDIR=$SUBJALLDIR/$BLPREFIX/ASHST1_MTLCORTEX_MSTTHK
  MSTVTTHK=""
  for side in left right; do
    MEANTHKFN=$SUBJMSTTHKDIR/${id}_${PREFIX}_${side}_thickness.csv
    SUBJMSTVTGEOSHOOTDIR=$SUBJMSTTHKDIR/GeoShoot/$side
    MeanTHK=""
    if [[ -f $MEANTHKFN ]]; then
      MeanTHK=$(cat -A $MEANTHKFN | head -n 2 | tail -n 1 | cut -d , -f 5-8 | sed -e 's/\$//g')
    elif [[ -f $SUBJMSTVTGEOSHOOTDIR/${id}_${side}_mean_thickness.vtk ]]; then
      MeanTHK=$(cat -A $SUBJMSTVTGEOSHOOTDIR/${id}_${side}_mean_thickness.vtk | cut -d , -f 3-6 | sed -e 's/\$//g')
    fi
    if [[ $MeanTHK == "" ]]; then
      MeanTHK=",,,"
    fi

    MEDIANTHKFN=$SUBJMSTTHKDIR/${id}_${PREFIX}_${side}_thickness.csv
    MedianTHK=""
    if [[ -f $MEDIANTHKFN ]]; then
      MedianTHK=$(cat -A $MEDIANTHKFN | head -n 2 | tail -n 1 | cut -d , -f 5-8 | sed -e 's/\$//g')
    elif [[ -f $SUBJMSTVTGEOSHOOTDIR/${id}_${side}_median_thickness.vtk ]]; then
      MedianTHK=$(cat -A $SUBJMSTVTGEOSHOOTDIR/${id}_${side}_median_thickness.vtk | cut -d , -f 3-6 | sed -e 's/\$//g')
    fi
    if [[ $MedianTHK == "" ]]; then
      MedianTHK=",,,"
    fi

    MSTVTTHK="$MSTVTTHK,$MeanTHK,$MedianTHK"
  done
  MSTVTTHK="${MSTVTTHK:1}"
  # compute mean
  for ((i=1;i<=8;i++)); do
    LMEA=$(echo $MSTVTTHK | cut -d, -f $i)
    RMEA=$(echo $MSTVTTHK | cut -d, -f $((i+8)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=10;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    MSTVTTHK="$MSTVTTHK,$MMEA"
  done

  #echo "$BLROW,$ICV,$VOL,$MSTVTTHK,$QCMRI,$QCICV,$QCCORTEXL,$QCCORTEXR,$QCCORTEXB,$QCHIPPOL,$QCHIPPOR,$QCHIPPOB,-1" >> $OUTTMPDIR/$(printf %04d $ii)_${id}.csv
  echo "$BLROW,$ICV,$VOL,$MSTVTTHK,-1" >> $OUTTMPDIR/$(printf %04d $ii)_${id}.csv

  set -e

}


######################################################
function MSTMTUTThicknessForMisMatchCohort()
{
  # Load id number
  ALLINFO_TXT=$ANALYSISDIR/Mismatch_MRIDates.csv
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=10
  N_begin=2

  echo "" >> $LOGDIR/MisMatchMSTThkASHST1_log.txt
  echo $(date) >> $LOGDIR/MisMatchMSTThkASHST1_log.txt

  # columns
  #RIDCol=1
  #IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  #Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  #ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  #BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  #TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT MRISCANDATE)

  # Submit job to copy data
  JPREFIX=MSTTHK${expid}
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
    SUBJALLDIR=$ALLDATADIR/$id
    #SUBJDIR=$DATADIR/$id

    # skip if it is longitudinal
    echo "MST thickness ASHS T1 Processing: $i,$rid,$PREFIX $Type"
    SUBJASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
    SUBJMSTTHKDIR=$SUBJALLDIR/$PREFIX/ASHST1_MTLCORTEX_MSTTHK

    # submit jobs
    for side in left right; do
      AUTOSEG=$SUBJASHST1DIR/final/${id}_${side}_lfseg_heur.nii.gz
      if [[ -f $AUTOSEG ]]; then
        # report status
        echo "   Adding $rid $PREFIX $side"
        echo "MST thickness ASHS T1: Adding $rid $PREFIX $side" >> $LOGDIR/MisMatchMSTThkASHST1_log.txt
      #qsub -cwd -o $DUMPDIR -j y \
      #     -p -1023 \
      #     -q all.q \
      #     -l h_vmem=8.1G,s_vmem=8G \
      #     -N "${JPREFIX}_${id}_${PREFIX}_${side}" $0 \
      #     MSTMTUTThickness_sub $id $PREFIX $side $AUTOSEG $SUBJMSTTHKDIR
      else
        echo "skip!"
      fi
    done

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}


######################################################
function CleanUpSpace()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=1510
  #N_begin=1500
  N_begin=1000

  # put date in log
  LOGFILE=$LOGDIR/cleanup_space_log.txt
  echo "" >> $LOGFILE
  echo $(date) >> $LOGFILE

  # columns
  ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)

  # Submit job to copy data
  PREFIX=SM
  for ((ii=$N_begin;ii<=${N};ii++)); do

    ALLINFO_TXT=$OUTPUTDIR/checkdata.csv
    ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
    
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

    echo "ASHS T1 Processing: $i,$rid,$PREFIX"
    if [[ $Type != "T1Baseline" && $Type != "T1Longitudinal" ]]; then
      continue
    fi

    # clean up step 1: convert the T1 SR image to short datatype
    if [[ -f $FUT1SRIMG ]]; then
      datatype=$(c3d $FUT1SRIMG -info-full | grep datatype \
                     | cut -d = -f 2 | cut -d " " -f 2)
      if [[ $datatype -gt 8 ]]; then
        echo "    Changing T1 SR datatype for $rid $PREFIX"
        echo "Clean up space: Changing T1 SR datatype for $rid $PREFIX" >> $LOGFILE
        c3d $FUT1SRIMG -type short -o $FUT1SRIMG
      fi
    fi

    # clean up step 2: remove the old multi-template thickness folder
    SUBJOLDMTTHKDIR=$SUBJALLDIR/$PREFIX/ASHST1_MTLCORTEX_MultiTempThk
    if [[ -d $SUBJOLDMTTHKDIR ]]; then
      echo "    Removing ASHST1_MTLCORTEX_MultiTempThk for $rid $PREFIX"
      echo "Clean up space: Removing ASHST1_MTLCORTEX_MultiTempThk for $rid $PREFIX" >> $LOGFILE
      rm -rf $SUBJOLDMTTHKDIR
    fi

    
    # cleanup MST multi-template thickness
    SUBJMSTTHKDIR=$SUBJALLDIR/$PREFIX/ASHST1_MTLCORTEX_MSTTHK
    SUBJMSTDATADIR=$SUBJMSTTHKDIR/data
    SUBJMSTREGATLASDIR=$SUBJMSTTHKDIR/RegToAtlases
    SUBJMSTVTREGDIR=$SUBJMSTTHKDIR/RegToInitTemp
    SUBJMSTUTREGDIR=$SUBJMSTTHKDIR/RegToUT
    if [[ -d $SUBJMSTDATADIR || -d $SUBJMSTREGATLASDIR || -d $SUBJMSTVTREGDIR || -d $SUBJMSTUTREGDIR ]]; then
      echo "    Removing intermediate files in the MST folder for $rid $PREFIX"
      echo "Clean up space: Removing intermediate files in the MST folder for $rid $PREFIX" >> $LOGFILE
      rm -rf $SUBJMSTDATADIR $SUBJMSTREGATLASDIR $SUBJMSTVTREGDIR $SUBJMSTUTREGDIR
    fi

  done
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
    BLScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT BLSCANDATE)
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
           -q all.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${i}_${rid}_${PREFIX}" \
           $0 MTLCortexRegionalThk_sub $i \

    fi
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
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
  BLScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT BLSCANDATE)
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
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
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
           -q all.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${id}_${side}" \
           $0 ComputeRegionalAtropyRate_sub $id $side
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
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
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
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
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
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
          -q all.q \
          -l h_vmem=4.1G,s_vmem=4G \
          -N "${PREFIX}_${exp}_${side}_${grp}_Atrophy" \
          $0 RegionalStats_meshglm_sub $exp $side $grp
        sleep 0.1
      done
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
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
    BLScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT BLSCANDATE)
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
           -q all.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${i}_${rid}_${PREFIX}" \
           $0 MTLCortexRegionalThkMultiTemp_sub $i \

    fi
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
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
  BLScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT BLSCANDATE)
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
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
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
           -q all.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${JPREFIX}_${id}_${side}" \
           $0 ComputeRegionalAtropyRateMultiTemp_sub $id $side
    done
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
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
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
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
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
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
          -q all.q \
          -l h_vmem=4.1G,s_vmem=4G \
          -N "${PREFIX}_${exp}_${side}_Atrophy" \
          $0 UTRegionalStats_meshglm_sub $exp $side
        sleep 0.1
      done
    done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
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
         -q all.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX1}_${i}" \
         $0 SummarizeMultiTemp_sub $i $OUTTMPDIR
    sleep 0.01

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
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
  BLScanDateCol=$(csvcol_tab.sh $ALLINFO_TXT BLSCANDATE)
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
function VoxelwiseLongiRobin()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/longitudinal_Robin.csv
  #N=$(cat $ALLINFO_TXT | wc -l)
  #N_begin=2
  #N=1064
  rm -f $OUTPUTDIR/log.txt

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)
  SUCCESS=$(csvcol.sh $ALLINFO_TXT ALOHA_Success)
  DateDiffCol=$(csvcol.sh $ALLINFO_TXT DateDiffFromBLScanDate)

  # IDs
  IDs=($(cat $ALLINFO_TXT | awk -F, '{print $2}' | uniq))
  N=${#IDs[*]}
  #N=2
  N_begin=0

  # Submit job for each subject
  JPREFIX=VTHK${expid}
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
    qsub -cwd -o $DUMPDIR -j y \
         -p -1023 \
         -q all.q \
         -l h_vmem=8.1G,s_vmem=8G \
         -N "${JPREFIX}_${id}" \
         $0 VoxelwiseLongiRobin_sub $id

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${JPREFIX}_*" -sync y -b y \
       sleep 1
}

function VoxelwiseLongiRobin_sub()
{
  id=$1
  side=$2
  ALLINFO_TXT=$OUTPUTDIR/longitudinal_Robin.csv
  SUBJALLDIR=$ALLDATADIR/$id
  VWLONGIOUTDIR=$SUBJALLDIR/T1_WB_longSingleSubjectTemplate/VoxelwiseLongitudinal
  rm -rf $VWLONGIOUTDIR
  mkdir -p $VWLONGIOUTDIR

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
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

  #########################################
  # gather scans and information for processing
  # get baseline information
  case=${cases[0]}
  rid=$(echo $case | cut -f $RIDCol -d ",")
  blscandate=$(echo $case | cut -f $BLScanDateCol -d ",")
  blscandate=$(ReFormateDate $blscandate)
  DD=$(date -d "$blscandate" '+%d')
  MM=$(date -d "$blscandate" '+%m')
  YYYY=$(date -d "$blscandate" '+%Y')
  BLPREFIX="${YYYY}-${MM}-${DD}"
  BLTHKSCAN=$(find $SUBJALLDIR/${BLPREFIX}_${id}_T1w*/${BLPREFIX}_${id}_T1w*CorticalThickness.nii.gz) 
  NORMBLTHKSCAN=$VWLONGIOUTDIR/${BLPREFIX}_${id}_T1w_trimCorticalThicknessNormalizedToTemplate_smooth${smooth}.nii.gz
  echo "$id,$BLPREFIX,0,$BLTHKSCAN,$NORMBLTHKSCAN,-1" > $VWLONGIOUTDIR/timepoint_info.csv

  # get all the meshes
  for ((i=0;i<${#cases[*]};i++)); do

    case=${cases[i]}
    date_diff=$(echo $case | cut -f $DateDiffCol -d ",")
    scandate=$(echo $case | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    THKSCAN=$(find $SUBJALLDIR/${PREFIX}_${id}_T1w_*/${PREFIX}_${id}_T1w*CorticalThickness.nii.gz)
    NORMTHKSCAN=$VWLONGIOUTDIR/${PREFIX}_${id}_T1w_trimCorticalThicknessNormalizedToTemplate_smooth${smooth}.nii.gz
    echo "$id,$PREFIX,$date_diff,$THKSCAN,$NORMTHKSCAN,-1" >> $VWLONGIOUTDIR/timepoint_info.csv

  done

  #########################################
  # warp thickness to OASIS template (not SST)
  cat $VWLONGIOUTDIR/timepoint_info.csv |
    while read line; do
      
      # info
      PREFIX=$(echo $line | cut -d, -f 2)
      THKSCAN=$(echo $line | cut -d, -f 4)
      THKSCANDIR=$(dirname $THKSCAN)
      NORMTHKSCAN=$(echo $line | cut -d, -f 5)
      SSTDIR=$SUBJALLDIR/T1_WB_longSingleSubjectTemplate

      # use antsApplyTransforms to warp thickness map to OASIS template (code borrowed rom antsCorticalThickness.sh)
      antsApplyTransforms -d 3 \
        -i $THKSCAN \
        -o $NORMTHKSCAN \
        -r $TEMPLATEDIR/T_template0.nii.gz \
        -n Gaussian \
        --float 1 \
        -t $SSTDIR/T_templateSubjectToTemplate1Warp.nii.gz \
        -t $SSTDIR/T_templateSubjectToTemplate0GenericAffine.mat \
        -t $THKSCANDIR/${PREFIX}_${id}_T1w*SubjectToTemplate1Warp.nii.gz \
        -t $THKSCANDIR/${PREFIX}_${id}_T1w*SubjectToTemplate0GenericAffine.mat

      c3d $NORMTHKSCAN -smooth 1vox -o $NORMTHKSCAN

    done 
  

  #########################################
  # call MATLAB script to compute regional atropy map
  OUTPREFIX=$VWLONGIOUTDIR/${id}_voxelwise_atrophyrate_${ylthresh}-${yhthresh}yrs_smooth${smooth}_leastsquare
  MASK=$TEMPLATEDIR/T_template0_gmmask.nii.gz
  robust=0
  $MATLAB_BIN -nosplash -nodesktop <<-MATCODE
    addpath('/home/longxie/ADNI2018/matlabcode');
    compute_voxelwise_annualized_atrophy_rate('$VWLONGIOUTDIR/timepoint_info.csv','$MASK','$OUTPREFIX',$robust);
MATCODE

  OUTPREFIX=$VWLONGIOUTDIR/${id}_voxelwise_atrophyrate_${ylthresh}-${yhthresh}yrs_smooth${smooth}_robust
  robust=1
  $MATLAB_BIN -nosplash -nodesktop <<-MATCODE
    addpath('/home/longxie/ADNI2018/matlabcode');
    compute_voxelwise_annualized_atrophy_rate('$VWLONGIOUTDIR/timepoint_info.csv','$MASK','$OUTPREFIX',$robust);
MATCODE
}

######################################################
function CheckVoxelwiseResult()
{
  # Load id number
  ALLINFO_TXT=$OUTPUTDIR/longitudinal_Robin.csv
  #N=$(cat $ALLINFO_TXT | wc -l)
  #N_begin=2
  #N=1064
  OUTCSV=$OUTPUTDIR/check_voxelwise_LongiRobin_${ylthresh}-${yhthresh}yrs_smooth${smooth}.csv
  echo "ID,NTP,MaxDateDiff,AllDateDiff,File1,File2,File3,File4" > $OUTCSV


  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  Vis2Col=$(csvcol.sh $ALLINFO_TXT VISCODE2)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSCANDATE)
  TYPECol=$(csvcol.sh $ALLINFO_TXT AnalysisType)
  SUCCESS=$(csvcol.sh $ALLINFO_TXT ALOHA_Success)
  DateDiffCol=$(csvcol.sh $ALLINFO_TXT DateDiffFromBLScanDate)

  # IDs
  IDs=($(cat $ALLINFO_TXT | awk -F, '{print $2}' | uniq))
  N=${#IDs[*]}
  #N=2
  N_begin=0

  # Submit job for each subject
  JPREFIX=VTHK${expid}
  for ((i=${N_begin};i<${N};i++)); do

    id=${IDs[$i]}

    # get cases
    set +e
    fcases=($(cat $ALLINFO_TXT | grep $id | grep T1Longitudinal | sed -e 's/\ //g'))
    set -e
    # if no case is found, skip
    if [[ ${#fcases[*]} == 0 ]]; then
      echo "$id,0,,,,,," >> $OUTCSV
      continue
    fi

    # check how many qualified cases
    cases=""
    idx=0
    all_date_diff=""
    max_date_diff=""
    lthresh=$(printf "%.0f" $(echo "scale=2;${ylthresh}*365" | bc))
    hthresh=$(printf "%.0f" $(echo "scale=2;${yhthresh}*365" | bc))
    for ((j=0;j<${#fcases[*]};j++)); do
      fcase=${fcases[j]}
      date_diff=$(echo $fcase | cut -f $DateDiffCol -d ",")
      success=$(echo $fcase | cut -f $SUCCESS -d ",")
      if [[ $date_diff -le $hthresh && $date_diff -ge $lthresh && $success -eq 1 ]]; then
        cases[${idx}]=$fcase
        max_date_diff=$date_diff
        all_date_diff="${all_date_diff};${date_diff}"
        idx=$((idx+1))
      fi
    done

    # if no qualified cases, skip
    if [[ $cases == "" ]]; then
      echo "$id,0,$max_date_diff,$all_date_diff,,,," >> $OUTCSV
      continue
    fi

    # output information
    ntp=${#cases[*]}
    ntp=$((ntp+1))
    
    # check whether the file has been generated
    SUBJALLDIR=$ALLDATADIR/$id
    VWLONGIOUTDIR=$SUBJALLDIR/T1_WB_longSingleSubjectTemplate/VoxelwiseLongitudinal
    set +e
    OUTPREFIX=$VWLONGIOUTDIR/${id}_voxelwise_atrophyrate_${ylthresh}-${yhthresh}yrs_smooth${smooth}_leastsquare
    FILES="$(find ${OUTPREFIX}_ChangeAnnualized.nii.gz),$(find ${OUTPREFIX}_PercChangeAnnualized.nii.gz)"
    OUTPREFIX=$VWLONGIOUTDIR/${id}_voxelwise_atrophyrate_${ylthresh}-${yhthresh}yrs_smooth${smooth}_robust
    FILES="$FILES,$(find ${OUTPREFIX}_ChangeAnnualized.nii.gz),$(find ${OUTPREFIX}_PercChangeAnnualized.nii.gz)"
    set -e    

    # output
    echo "$id,$((idx+1)),$max_date_diff,$all_date_diff,$FILES" >> $OUTCSV

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
