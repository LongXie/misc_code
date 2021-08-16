#//bin/bash
#$ -S /bin/bash
set -e
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
ROOT=/data/picsl/longxie/ADNI2018/
CODEDIR=$ROOT/scripts
ANALYSISDIR=$ROOT/analysis_input
#ALLINFO_TXT=$ANALYSISDIR/MRI3TListWithNIFTIPath_11262018.tsv
ALLINFO_TXT=$ANALYSISDIR/ADNI_T2MRI_3T_TextureAnalysis_09092019.csv
#ALLINFO_TXT=$ANALYSISDIR/MRILIST_T1T2_012_S_4128.tsv
ALLDATADIR=$ROOT/dataset
LOCALDATADIR=$ROOT/dataset
OUTPUTDIR=$ROOT/output_ADNIT2Texture_09092019
LOGDIR=$OUTPUTDIR/log
DUMPDIR=$OUTPUTDIR/dump
mkdir -p $LOGDIR

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

#############################################
function main()
{
  reset_dir

  # perform ASHS T2 segmentation if necessary
  #ASHST2SegSUB
  #MakeASHST2QCPNG
  #CleanUpASHST2
  
  #CopyASHST2ToDatasetDir


  # summarize ASHS information (both T1 and T2)
  #SummarizeASHS

  # Sample patches from the selected timepoints (the last timepoint of each subject)
  #SampleDBPatch

  # Sample isotropic patches
  #SampleIsoDBPatch

  # generate image dataset for learning
  #GenerateDLDataset

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
    nYYYY=$(date -d "$indate" '+%Y')
    outdate="${YYYY}-${MM}-${DD}"

  fi

  echo $outdate
}

######################################################
function ASHST2SegSUB()
{
  # Load id number
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=2
  ATLAS=/home/pauly/wolk/atlas2016/ashs01/final
  Nrun=20
  echo "" >> $LOGDIR/ASHST2_log.txt
  echo "$(date)" >> $LOGDIR/ASHST2_log.txt

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE_T2MRI_3T_List)
  PrimaryLastCol=$(csvcol.sh $ALLINFO_TXT PrimaryLast)

  # Submit job to copy data
  PREFIX1=ASHST2${expid}
  for ((i=2;i<=${N};i++)); do

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
    last=$(echo $ROW | cut -f $PrimaryLastCol -d ",")
    SUBJALLDIR=$ALLDATADIR/$id
    T1=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T1w_trim.nii.gz
    T2=$SUBJALLDIR/$PREFIX/${PREFIX}_${id}_T2w.nii.gz
    SUBJASHST2DIR=$SUBJALLDIR/${PREFIX}/sfsegnibtend
    OUTDIR=$SUBJASHST2DIR

    # start segmentation
    echo "ASHS T2 Processing: $i,$rid,$PREFIX"
    if [[ $last == 0 ]]; then
      continue
    fi
    #SUBJLOCALDIR=$LOCALDATADIR/$id
    #OUTDIR=$SUBJLOCALDIR/${PREFIX}/sfsegnibtend

    if [[ -f $T1 && -f $T2 && \
       ! -f $SUBJASHST2DIR/final/${id}_left_lfseg_corr_nogray.nii.gz  ]]; then

      echo "   Adding $rid $PREFIX "
      echo "ASHS T2: Adding $rid $PREFIX ($(date)) " >> $LOGDIR/ASHST2_log.txt

      mkdir -p $OUTDIR
      $ASHS_ROOT/bin/ashs_main.sh \
        -a $ATLAS -d -T -I $id -g $T1 -f $T2 \
        -Q -s 1-7 \
        -z $CODEDIR/ashs-fast-z.sh \
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
function MakeASHST2QCPNG()
{
  # Load id number
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=10
  #N_begin=1054
  N_begin=2

  # output directory
  QCDIR=$OUTPUTDIR/qc
  rm -rf $QCDIR
  mkdir -p $QCDIR

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE_T2MRI_3T_List)
  PrimaryLastCol=$(csvcol.sh $ALLINFO_TXT PrimaryLast)
  idxall=0
  idx=0

  # Submit job to copy data
  PREFIXJ=MPNG
  CMD=""
  for ((ii=$N_begin;ii<=${N};ii++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    last=$(echo $ROW | cut -f $PrimaryLastCol -d ",")

    # directories
    SUBJALLDIR=$ALLDATADIR/$id
    SUBJASHST2DIR=$SUBJALLDIR/${PREFIX}/sfsegnibtend
    SUBJLOCALDIR=$LOCALDATADIR/$id
    LOCALTPDIR=$SUBJLOCALDIR/$PREFIX
    if [[ ! -d $SUBJASHST2DIR ]]; then
      SUBJASHST2DIR=$LOCALTPDIR/sfsegnibtend
    fi
    LSEG=$SUBJASHST2DIR/final/${id}_left_lfseg_corr_nogray.nii.gz
    RSEG=$SUBJASHST2DIR/final/${id}_right_lfseg_corr_nogray.nii.gz

    if [[ ! -f $LSEG || ! -f $RSEG ]]; then

      echo "Segmentation files from $id are missing"

    else

      # montage command line
      CMD="$CMD -label ${id}_${PREFIX}_left $SUBJASHST2DIR/qa/qa_seg_bootstrap_corr_usegray_left_qa.png -geometry +10+10 -label ${id}_${PREFIX}_right $SUBJASHST2DIR/qa/qa_seg_bootstrap_corr_usegray_right_qa.png -geometry +10+10"
      idx=$((idx+1))

    fi
    echo $idx

    if [ $idx -gt 100 ]; then
      echo "$idxall   $idx  "
      montage $CMD $QCDIR/qa_${idxall}.png
      CMD=""
      idx=0
      idxall=$((idxall+1))
    fi

  done

  if [[ $CMD != "" ]]; then
    echo "$idxall   $idx  "
    montage $CMD $QCDIR/qa_${idxall}.png
  fi
}

##################################################
function CopyASHST2ToDatasetDir()
{
  # Load id number
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=2
  #N_begin=1054
  N_begin=2
  echo "" >> $LOGDIR/copyASHST2.txt
  echo "$(date)" >> $LOGDIR/copyASHST2.txt

  # output directory
  QCDIR=$OUTPUTDIR/qc
  rm -rf $QCDIR
  mkdir -p $QCDIR

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE_T2MRI_3T_List)
  PrimaryLastCol=$(csvcol.sh $ALLINFO_TXT PrimaryLast)
  idxall=0
  idx=0
 
  # Submit job to copy data
  PREFIXJ=MPNG
  CMD=""
  for ((ii=$N_begin;ii<=${N};ii++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    last=$(echo $ROW | cut -f $PrimaryLastCol -d ",")

    # status
    echo "Copying output to dataset directory: $i, $rid, $PREFIX"

    # directories
    SUBJALLDIR=$ALLDATADIR/$id
    SUBJASHST2DIR=$SUBJALLDIR/${PREFIX}/sfsegnibtend
    SUBJT2=$SUBJALLDIR/${PREFIX}/${PREFIX}_${id}_T2w.nii.gz
    SUBJLOCALDIR=$LOCALDATADIR/$id
    SUBJLOCALASHST2DIR=$SUBJLOCALDIR/$PREFIX/sfsegnibtend
    SUBJSDDIR=/home/srdas/wd/ADNI23/$id/
    SUBJSDASHST2DIR=$SUBJSDDIR/$PREFIX/sfsegnibtend

    # check whether it has been processed
    ASHST2DIR=$SUBJASHST2DIR
    LSEG=$ASHST2DIR/final/${id}_left_lfseg_corr_nogray.nii.gz
    RSEG=$ASHST2DIR/final/${id}_right_lfseg_corr_nogray.nii.gz
    if [[ -f $LSEG && -f $RSEG ]]; then
      echo "    $ii, $rid, $PREFIX is  already processed in dataset directory"
      continue
    fi

    # check Sandy's directory
    ASHST2DIR=$SUBJSDASHST2DIR
    LSEG=$ASHST2DIR/final/${id}_left_lfseg_corr_nogray.nii.gz
    RSEG=$ASHST2DIR/final/${id}_right_lfseg_corr_nogray.nii.gz
    #DIFF=""
    if [[ -f $LSEG && -f $RSEG ]]; then
      set +e 
      SDT2=$ASHST2DIR/tse.nii.gz

      RAWDIFF=$(c3d $SDT2 $SUBJT2 -scale -1 -add -voxel-sum | cut -d " " -f 3)
      if [[ $RAWDIFF == 0 ]]; then
        msg="    $ii, $rid, $PREFIX copying data from Sandy's directory."
        echo $msg
        echo $msg >> $LOGDIR/copyASHST2.txt
        mkdir -p $SUBJASHST2DIR/
        cp -r $ASHST2DIR/final/ $SUBJASHST2DIR/
        cp -r $ASHST2DIR/qa $SUBJASHST2DIR/
        cp -r $ASHST2DIR/*native_chunk* $SUBJASHST2DIR/
        echo "Copied from $ASHST2DIR/" > $SUBJASHST2DIR/log.txt
        continue
      fi

      # for those that are not exactly the same, do a more detail check
      orient_code=$(c3d $SUBJT2 -info | cut -d ';' -f 5 | cut -d ' ' -f 5)
      if [[ $orient_code == "Oblique," ]]; then
        orient_code=$(c3d $SUBJT2 -info | cut -d ';' -f 5 | cut -d ' ' -f 8)
      fi
      DIFF=($(c3d $SUBJT2 -dup $SDT2 -swapdim ${orient_code} \
                  -copy-transform -scale -1 -add -info \
                  | grep 'range' | sed -e "s/.*\[//g" -e "s/,/ /g" -e "s/\].*//g"))
      set -e
      if [[ $DIFF == "" ]]; then
        msg="    $ii, $rid, $PREFIX data available from Sandy's directory but different T2 scan was used.($SDT2, $SUBJT2, ${DIFF[*]})"
        echo $msg
        echo $msg >> $LOGDIR/copyASHST2.txt
      elif [[ $(echo "${DIFF[0]}>-5" | bc) == 1 && $(echo "${DIFF[1]}<5" | bc) == 1 ]]; then
        msg="    $ii, $rid, $PREFIX found small header difference, converting data from Sandy's directory (${DIFF[*]}). "
        echo $msg
        echo $msg >> $LOGDIR/copyASHST2.txt
        mkdir -p $SUBJASHST2DIR/
        cp -r $ASHST2DIR/final/ $SUBJASHST2DIR/
        cp -r $ASHST2DIR/qa $SUBJASHST2DIR/
        cp -r $ASHST2DIR/*native_chunk* $SUBJASHST2DIR/
        echo "Converted from $ASHST2DIR/" > $SUBJASHST2DIR/log.txt
        for S in left right; do
          for T in heur corr_nogray corr_usegray; do
            SEG=$SUBJASHST2DIR/final/${id}_${S}_lfseg_${T}.nii.gz
            c3d $SUBJT2 $SEG -swapdim ${orient_code} -copy-transform -o $SEG
          done
          c3d $SUBJT2 -as T2 $SUBJASHST2DIR/tse_native_chunk_${S}.nii.gz \
              -int 0 -reslice-identity -thresh 1 inf 1 0 -trim 0x0x0vox \
              -push T2 -int 0 -reslice-identity \
              -o $SUBJASHST2DIR/tse_native_chunk_${S}.nii.gz          
        done
        continue
      else
        msg="    $ii, $rid, $PREFIX data available from Sandy's directory but different T2 scan was used.($SDT2, $SUBJT2, ${DIFF[*]})"
        echo $msg
        echo $msg >> $LOGDIR/copyASHST2.txt
      fi
    fi 

    # check my local directory
    ASHST2DIR=$SUBJLOCALASHST2DIR
    LSEG=$ASHST2DIR/final/${id}_left_lfseg_corr_nogray.nii.gz
    RSEG=$ASHST2DIR/final/${id}_right_lfseg_corr_nogray.nii.gz
    if [[ -f $LSEG && -f $RSEG ]]; then
      msg="    $ii, $rid, $PREFIX copying data from my local directory"
      echo $msg
      echo $msg >> $LOGDIR/copyASHST2.txt
      #mkdir -p $SUBJASHST2DIR/
      #cp -r $ASHST2DIR/final/ $SUBJASHST2DIR/
      #cp -r $ASHST2DIR/qa $SUBJASHST2DIR/
      #cp -r $ASHST2DIR/*native_chunk* $SUBJASHST2DIR/
      #echo "Copied from $ASHST2DIR/" > $SUBJASHST2DIR/log.tx
    else
      msg="    $ii, $rid, $PREFIX segmentation data not found"
      echo $msg
      echo $msg >> $LOGDIR/copyASHST2.txt
    fi

  done
}

######################################################
function SampleDBPatch()
{
  # Load id number
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=10
  #N_begin=1054
  N_begin=2

  # output directory
  OUTDIR=$OUTPUTDIR/patches_last
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR $OUTDIR
  mkdir -p $OUTTMPDIR
  mkdir -p $OUTDIR

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE_T2MRI_3T_List)
  PrimaryLastCol=$(csvcol.sh $ALLINFO_TXT PrimaryLast)

  # Submit job to copy data
  PREFIXJ=SP
  for ((ii=$N_begin;ii<=${N};ii++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    last=$(echo $ROW | cut -f $PrimaryLastCol -d ",")

    echo "Samping patches: $ii,$rid,$PREFIX"
    if [[ $last == "1" ]]; then
      echo "    submitting $ii,$rid,$PREFIX"
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIXJ}_${ii}" $0 \
           SampleDBPatch_sub $ii $OUTDIR $OUTTMPDIR
      sleep 0.05
    else
      echo "    skip."
    fi
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q,gpu.q \
       -hold_jid "${PREFIXJ}_*" -sync y -b y \
       sleep 10

  # get all info
  # header
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header="$header,ASHST2Success,QC"
  echo $header > $OUTPUTDIR/DBPatch_dataset.csv
  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/DBPatch_dataset.csv
}

function SampleDBPatch_sub()
{
  ii=$1
  OUTDIR=$2
  OUTTMPDIR=$3
  ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE_T2MRI_3T_List)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSANDATE)
  PrimaryLastCol=$(csvcol.sh $ALLINFO_TXT PrimaryLast)

  # get information
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  last=$(echo $ROW | cut -f $PrimaryLastCol -d ",")

  # directories
  SUBJALLDIR=$ALLDATADIR/$id
  SUBJASHST2DIR=$SUBJALLDIR/${PREFIX}/sfsegnibtend
  SUBJTPDIR=$SUBJALLDIR/${PREFIX}
  #SUBJLOCALDIR=$LOCALDATADIR/$id
  #LOCALTPDIR=$SUBJLOCALDIR/$PREFIX
  #LOCALPATCHDIR=$LOCALTPDIR/DBPatches/
  LOCALPATCHDIR=$TMPDIR
  mkdir -p $LOCALPATCHDIR
  #if [[ ! -d $SUBJASHST2DIR ]]; then
  #  SUBJASHST2DIR=$LOCALTPDIR/sfsegnibtend
  #else
  #  mkdir -p $LOCALTPDIR
  #  ln -sf $SUBJALLDIR/${PREFIX}/sfsegnibtend $LOCALTPDIR/
  #fi
  
  # link some scans
  #ln -sf $SUBJALLDIR/${PREFIX}/${PREFIX}_${id}_T1w_trim.nii.gz $LOCALTPDIR/
  #ln -sf $SUBJALLDIR/${PREFIX}/${PREFIX}_${id}_T2w.nii.gz $LOCALTPDIR/

  success=1
  for side in left right; do

    SEG=$SUBJASHST2DIR/final/${id}_${side}_lfseg_corr_nogray.nii.gz
    if [[ -f $SEG ]]; then

      # generate the mask for sampling region
      w=1
      mask=$LOCALPATCHDIR/${id}_${PREFIX}_${side}_mask_darkband_width${w}.nii.gz
      if [[ ! -f $mask ]]; then
      c3d $SEG -as SEG \
        -thresh 3 3 1 0 -popas DG \
        -push SEG -replace 1 999 2 99 4 999 7 999 8 999 \
        -thresh 999 999 1 0 -popas NONDG \
        -push DG -dilate 1 ${w}x${w}x0vox \
        -push NONDG -multiply -popas RIM1 \
        -push NONDG -dilate 1 ${w}x${w}x0vox \
        -push DG -multiply -push RIM1 -add \
        -push SEG -thresh 10 12 1 0 \
        -dilate 1 100x100x0vox \
        -multiply \
        -o ${mask}
      
      fi

      # perform sampling
      ps1=11
      ps2=2
      OUTDAT=$OUTDIR/${id}_${PREFIX}_${side}_PS${ps1}-${ps1}-${ps2}_NS1_NC2.dat
      NP=$(c3d $SUBJTPDIR/${PREFIX}_${id}_T2w.nii.gz $SEG \
             ${mask} \
             -xp $OUTDAT \
             ${ps1}x${ps1}x${ps2} 30 | cut -d : -f 2)
      OUTDATFINAL=$OUTDIR/${id}_${PREFIX}_${side}_PS${ps1}-${ps1}-${ps2}_NS1_NC2_NP${NP}.dat
      mv $OUTDAT $OUTDATFINAL

    else
      success=0
    fi
  done

  echo "$ROW,$success,0" > $OUTTMPDIR/$(printf %04d $ii)_${id}_${PREFIX}.csv
}

######################################################
function SampleIsoDBPatch()
{
  # Load id number
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=4
  #N_begin=1054
  N_begin=2

  # output directory
  OUTDIR=$OUTPUTDIR/isopatches_last
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR $OUTDIR
  mkdir -p $OUTTMPDIR
  mkdir -p $OUTDIR

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE_T2MRI_3T_List)
  PrimaryLastCol=$(csvcol.sh $ALLINFO_TXT PrimaryLast)

  # Submit job to copy data
  PREFIXJ=SP
  for ((ii=$N_begin;ii<=${N};ii++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    last=$(echo $ROW | cut -f $PrimaryLastCol -d ",")

    echo "Samping patches: $ii,$rid,$PREFIX"
    if [[ $last == "1" ]]; then
      echo "    submitting $ii,$rid,$PREFIX"
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIXJ}_${ii}" $0 \
           SampleIsoDBPatch_sub $ii $OUTDIR $OUTTMPDIR
      sleep 0.05
    else
      echo "    skip."
    fi
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q,gpu.q \
       -hold_jid "${PREFIXJ}_*" -sync y -b y \
       sleep 10

  # get all info
  # header
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header="$header,ASHST2Success,QC"
  echo $header > $OUTPUTDIR/IsoDBPatch_dataset.csv
  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/IsoDBPatch_dataset.csv
}

function SampleIsoDBPatch_sub()
{
  ii=$1
  OUTDIR=$2
  OUTTMPDIR=$3
  ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE_T2MRI_3T_List)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSANDATE)
  PrimaryLastCol=$(csvcol.sh $ALLINFO_TXT PrimaryLast)

  # get information
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  last=$(echo $ROW | cut -f $PrimaryLastCol -d ",")

  # directories
  SUBJALLDIR=$ALLDATADIR/$id
  SUBJASHST2DIR=$SUBJALLDIR/${PREFIX}/sfsegnibtend
  SUBJTPDIR=$SUBJALLDIR/${PREFIX}
  #SUBJLOCALDIR=$LOCALDATADIR/$id
  #LOCALTPDIR=$SUBJLOCALDIR/$PREFIX
  #LOCALPATCHDIR=$LOCALTPDIR/DBPatches/
  LOCALPATCHDIR=$TMPDIR
  mkdir -p $LOCALPATCHDIR
  #if [[ ! -d $SUBJASHST2DIR ]]; then
  #  SUBJASHST2DIR=$LOCALTPDIR/sfsegnibtend
  #else
  #  mkdir -p $LOCALTPDIR
  #  ln -sf $SUBJALLDIR/${PREFIX}/sfsegnibtend $LOCALTPDIR/
  #fi

  # link some scans
  #ln -sf $SUBJALLDIR/${PREFIX}/${PREFIX}_${id}_T1w_trim.nii.gz $LOCALTPDIR/
  #ln -sf $SUBJALLDIR/${PREFIX}/${PREFIX}_${id}_T2w.nii.gz $LOCALTPDIR/

  #TMPDIR=$SUBJASHST2DIR/final/

  success=1
  for side in left right; do

    SEG=$SUBJASHST2DIR/final/${id}_${side}_lfseg_corr_nogray.nii.gz
    if [[ -f $SEG ]]; then

      # generate the mask for sampling region
      w=1
      mask=$TMPDIR/${id}_${PREFIX}_${side}_mask_darkband_width${w}.nii.gz
      if [[ ! -f $mask ]]; then
      c3d $SEG -as SEG \
        -thresh 3 3 1 0 -popas DG \
        -push SEG -replace 1 999 2 99 4 999 7 999 8 999 \
        -thresh 999 999 1 0 -popas NONDG \
        -push DG -dilate 1 ${w}x${w}x0vox \
        -push NONDG -multiply -popas RIM1 \
        -push NONDG -dilate 1 ${w}x${w}x0vox \
        -push DG -multiply -push RIM1 -add \
        -push SEG -thresh 10 12 1 0 \
        -dilate 1 100x100x0vox \
        -multiply \
        -o ${mask}

      fi

      # resample T2 image to isotropic
      c3d $SUBJTPDIR/${PREFIX}_${id}_T2w.nii.gz \
        -resample 100x100x500% \
        -o $TMPDIR/${PREFIX}_${id}_T2w_iso.nii.gz
      c3d $TMPDIR/${PREFIX}_${id}_T2w_iso.nii.gz \
        $mask \
        -reslice-identity \
        -thresh 0.5 inf 1 0 \
        -o $mask

      # perform sampling
      ps1=11
      ps2=11
      OUTDAT=$OUTDIR/${id}_${PREFIX}_${side}_PS${ps1}-${ps1}-${ps2}_NS1_NC1.dat
      NP=$(c3d $TMPDIR/${PREFIX}_${id}_T2w_iso.nii.gz \
             ${mask} \
             -xp $OUTDAT \
             ${ps1}x${ps1}x${ps2} 120 | cut -d : -f 2)
      OUTDATFINAL=$OUTDIR/${id}_${PREFIX}_${side}_PS${ps1}-${ps1}-${ps2}_NS1_NC1_NP${NP}.dat
      mv $OUTDAT $OUTDATFINAL

    else
      success=0
    fi
  done

  echo "$ROW,$success,0" > $OUTTMPDIR/$(printf %04d $ii)_${id}_${PREFIX}.csv
}

######################################################
function GenerateDLDataset()
{
  # Load id number
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=4
  #N=4
  #N_begin=1054
  N_begin=2

  # output directory
  OUTDIR=$OUTPUTDIR/MTLPatchDataset
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR $OUTDIR
  mkdir -p $OUTTMPDIR
  mkdir -p $OUTDIR

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE_T2MRI_3T_List)
  PrimaryLastCol=$(csvcol.sh $ALLINFO_TXT PrimaryLast)

  # Submit job to copy data
  PREFIXJ=SP
  for ((ii=$N_begin;ii<=${N};ii++)); do

    ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
    rid=$(echo $ROW | cut -f $RIDCol -d ",")
    id=$(echo $ROW | cut -f $IDCol -d ",")
    scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
    DD=$(date -d "$scandate" '+%d')
    MM=$(date -d "$scandate" '+%m')
    YYYY=$(date -d "$scandate" '+%Y')
    PREFIX="${YYYY}-${MM}-${DD}"
    scandate=$(ReFormateDate $scandate)
    last=$(echo $ROW | cut -f $PrimaryLastCol -d ",")

    echo "Samping patches: $ii,$rid,$PREFIX"
    if [[ $last == "1" ]]; then
      for side in left right; do
        echo "    submitting $ii,$rid,$PREFIX,$side"
        SUBJOUTDIR=$OUTDIR/${id}_${side}
        qsub -cwd -o $DUMPDIR -j y \
             -q all.q,basic.q \
             -l h_vmem=4.1G,s_vmem=4G \
             -N "${PREFIXJ}_${ii}_${side}" $0 \
             GenerateDLDataset_sub $ii $side $SUBJOUTDIR $OUTTMPDIR
        sleep 0.05
      done
    else
      echo "    skip."
    fi
  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q,gpu.q \
       -hold_jid "${PREFIXJ}_*" -sync y -b y \
       sleep 10

  # get all info
  # header
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header="$header,ASHST2Success,QC"
  echo $header > $OUTPUTDIR/MTLPatch_dataset.csv
  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/MTLPatch_dataset.csv
}

function GenerateDLDataset_sub()
{
  ii=$1
  side=$2
  OUTDIR=$3
  OUTTMPDIR=$4
  mkdir -p $OUTDIR
  ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE_T2MRI_3T_List)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSANDATE)
  PrimaryLastCol=$(csvcol.sh $ALLINFO_TXT PrimaryLast)

  # get information
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  last=$(echo $ROW | cut -f $PrimaryLastCol -d ",")

  # directories
  SUBJALLDIR=$ALLDATADIR/$id
  SUBJASHST2DIR=$SUBJALLDIR/${PREFIX}/sfsegnibtend
  SUBJTPDIR=$SUBJALLDIR/${PREFIX}
  #TEMPDIR=$OUTDIR/${id}_${side}
  TEMPDIR=$TMPDIR
  mkdir -p $TEMPDIR

  success=1
  SEGORIG=$SUBJASHST2DIR/final/${id}_${side}_lfseg_corr_nogray.nii.gz
  if [[ -f $SEGORIG ]]; then

    # crop the segmentation
    SEG=$OUTDIR/${id}_${PREFIX}_${side}_lfseg_corr_nogray.nii.gz
    c3d $SEGORIG -trim 20x20x4vox -o $SEG
    
    # generate the mask for sampling region
    w=1
    mask=$OUTDIR/${id}_${PREFIX}_${side}_samplemask_width${w}.nii.gz
    c3d $SEG -as SEG \
      -thresh 3 3 1 0 -popas DG \
      -push SEG -replace 1 999 2 99 4 999 7 999 8 999 \
      -thresh 999 999 1 0 -popas NONDG \
      -push DG -dilate 1 ${w}x${w}x0vox \
      -push NONDG -multiply -popas RIM1 \
      -push NONDG -dilate 1 ${w}x${w}x0vox \
      -push DG -multiply -push RIM1 -add \
      -push SEG -thresh 10 10 2 0 -add \
      -push SEG -thresh 11 11 3 0 -add \
      -o ${mask}

    # resample T2 image to the cropped space
    T2=$OUTDIR/${id}_${PREFIX}_${side}_T2w.nii.gz
    c3d $SEG $SUBJTPDIR/${PREFIX}_${id}_T2w.nii.gz \
      -int 0 -reslice-identity \
      -o $T2
    c3d $SUBJTPDIR/${PREFIX}_${id}_T2w.nii.gz -resample 100x100x500% \
      -region 20x20x0% 60x60x100% \
      -type short -o $TEMPDIR/tse_iso.nii.gz
    greedy -d 3 -a -dof 6 -m MI -n 100x100x10 \
      -i $TEMPDIR/tse_iso.nii.gz $SUBJTPDIR/${PREFIX}_${id}_T1w.nii.gz \
      -ia-identity \
      -o $TEMPDIR/greedy_t2_to_t1.mat

    # register the T1 to T2 cropped space
    T1=$OUTDIR/${id}_${PREFIX}_${side}_T1w_reg.nii.gz
    greedy -d 3 \
      -rf $T2 \
      -rm $SUBJTPDIR/${PREFIX}_${id}_T1w.nii.gz \
          $TEMPDIR/${id}_${PREFIX}_${side}_mprage_to_tse_init.nii.gz\
      -r $TEMPDIR/greedy_t2_to_t1.mat
    greedy -d 3 \
      -a -dof 12 -m MI -n 20 \
      -i $T2 $TEMPDIR/${id}_${PREFIX}_${side}_mprage_to_tse_init.nii.gz \
      -o $TEMPDIR/${id}_${PREFIX}_${side}_t2_to_t1.mat
    greedy -d 3 \
      -rf $T2 \
      -rm $SUBJTPDIR/${PREFIX}_${id}_T1w.nii.gz \
          $T1 \
      -r $TEMPDIR/${id}_${PREFIX}_${side}_t2_to_t1.mat \
         $TEMPDIR/greedy_t2_to_t1.mat


  else
    success=0
  fi

  echo "$ROW,$success,0" > $OUTTMPDIR/$(printf %04d $ii)_${id}_${PREFIX}.csv
}







######################################################
function SummarizeASHS()
{
  # Load id number
  N=$(cat $ALLINFO_TXT | wc -l)
  #N=5
  #N_begin=1054
  N_begin=2

  # output directory
  OUTTMPDIR=$OUTPUTDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR

  # Submit job to copy data
  if [[ 1 == 1 ]]; then
  PREFIX=SM
  for ((ii=$N_begin;ii<=${N};ii++)); do

    echo "$i submitted"
    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q,gpu.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX}_${ii}" \
         $0 SummarizeASHS_sub $ii $OUTTMPDIR
    sleep 0.05

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q,gpu.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 10

  fi

  # get all info
  # header
  #blheader=$(cat -A $DEMOG_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  header=$(cat -A $ALLINFO_TXT | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  #sed -e 's/\^M\$//g'
  header="$header,ICV"

  for side in L R M; do
  for sub in AHippo PHippo Hippo ERC BA35 BA36 PHC CS OTS AMEN PMEN MISC; do
    header="$header,${side}_${sub}_VOL_ASHST1"
  done
  done

  for method in MSTVT; do
  for side in L R M; do
  for type in MeanTHK MedianTHK; do
  for sub in ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_${type}_${method}"
  done
  done
  done
  done

  for side in L R M; do
  for sub in CA1 CA2 DG CA3 SUB ERC BA35 BA36 PHC; do
    header="$header,${side}_${sub}_NS_ASHST2,${side}_${sub}_VOL_ASHST2"
  done
  done

  # output
  echo $header > $OUTPUTDIR/summarize_ASHS.csv

  cat $OUTTMPDIR/*.csv >> $OUTPUTDIR/summarize_ASHS.csv
}

function SummarizeASHS_sub()
{
  ii=$1
  OUTTMPDIR=$2
  ROW=$(cat -A $ALLINFO_TXT | head -n $ii | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')

  # columns
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO_TXT ID)
  ScanDateCol=$(csvcol.sh $ALLINFO_TXT SCANDATE_T2MRI_3T_List)
  BLScanDateCol=$(csvcol.sh $ALLINFO_TXT BLSANDATE)

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

  # baseline information
  blscandate=$(echo $ROW | cut -f $BLScanDateCol -d ",")
  blscandate=$(ReFormateDate $blscandate)
  DD=$(date -d "$blscandate" '+%d')
  MM=$(date -d "$blscandate" '+%m')
  YYYY=$(date -d "$blscandate" '+%Y')
  BLPREFIX="${YYYY}-${MM}-${DD}"
  BLSUBJDIR=$SUBJALLDIR/$BLPREFIX

  set +e

  #####################################################
  # get ICV
  SUBJASHSICVDIR=$BLSUBJDIR/ASHSICV
  ICV=$(cat $SUBJASHSICVDIR/final/${id}_left_corr_nogray_volumes.txt | awk -F' ' '{ print $5 }')

  #####################################################
  # get ASHS T1 volume
  SUBJASHST1DIR=$SUBJALLDIR/$PREFIX/ASHST1
  T1VOL=""
  for side in left right; do
    idside=${id}_${side}
    VOL_AHIPPO=$(cat $SUBJASHST1DIR/final/${idside}_heur_volumes.txt | grep Anterior_hippocampus | cut -d ' ' -f 5)
    VOL_PHIPPO=$(cat $SUBJASHST1DIR/final/${idside}_heur_volumes.txt | grep Posterior_hippocampus | cut -d ' ' -f 5)
    VOL_HIPPO=$(echo "scale=3;$VOL_AHIPPO+$VOL_PHIPPO" | bc -l)
    VOL_PHC=$(cat $SUBJASHST1DIR/final/${idside}_heur_volumes.txt | grep " PHC" | cut -d ' ' -f 5 )
    VOL_MISC=$(cat $SUBJASHST1DIR/final/${idside}_heur_volumes.txt | grep MISC | cut -d ' ' -f 5)
    VOL_ERC=$(cat $SUBJASHST1DIR/final/${idside}_heur_volumes.txt | grep ERC | cut -d ' ' -f 5)
    VOL_BA35=$(cat $SUBJASHST1DIR/final/${idside}_heur_volumes.txt | grep Br35 | cut -d ' ' -f 5)
    VOL_BA36=$(cat $SUBJASHST1DIR/final/${idside}_heur_volumes.txt | grep Br36 | cut -d ' ' -f 5)
    VOL_CS=$(cat $SUBJASHST1DIR/final/${idside}_heur_volumes.txt | grep ColSul | cut -d ' ' -f 5)
    VOL_OTS=$(cat $SUBJASHST1DIR/final/${idside}_heur_volumes.txt | grep OTSul | cut -d ' ' -f 5)
    VOL_AMEN=$(cat $SUBJASHST1DIR/final/${idside}_heur_volumes.txt | grep "Meninges " | cut -d ' ' -f 5)
    VOL_PMEN=$(cat $SUBJASHST1DIR/final/${idside}_heur_volumes.txt | grep "Meninges_PHC" | cut -d ' ' -f 5)
    T1VOL="$T1VOL,$VOL_AHIPPO,$VOL_PHIPPO,$VOL_HIPPO,$VOL_ERC,$VOL_BA35,$VOL_BA36,$VOL_PHC,$VOL_CS,$VOL_OTS,$VOL_AMEN,$VOL_PMEN,$VOL_MISC"
  done
  T1VOL=$(echo $T1VOL | sed 's/.\(.*\)/\1/')
  # compute mean
  ALLMEA=$T1VOL
  NL=12
  for ((i=1;i<=$NL;i++)); do
    LMEA=$(echo $ALLMEA | cut -d, -f $i)
    RMEA=$(echo $ALLMEA | cut -d, -f $((i+$NL)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=3;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    ALLMEA="$ALLMEA,$MMEA"
  done
  T1VOL=$ALLMEA

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
  ALLMEA=$MSTVTTHK
  NL=8
  for ((i=1;i<=$NL;i++)); do
    LMEA=$(echo $ALLMEA | cut -d, -f $i)
    RMEA=$(echo $ALLMEA | cut -d, -f $((i+$NL)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=3;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    ALLMEA="$ALLMEA,$MMEA"
  done
  MSTVTTHK=$ALLMEA

  #####################################################
  # get ASHS T2 volume and number of slides
  SUBJASHST2DIR=$SUBJALLDIR/$PREFIX/sfsegnibtend
  if [[ ! -d $SUBJASHST2DIR ]]; then
    SUBJASHST2DIR=$LOCALDATADIR/$id/$PREFIX/sfsegnibtend
  fi
  T2VOL=""
  for side in left right; do
    idside=${id}_${side}
    for sub in CA1 CA2 DG CA3 SUB ERC BA35 BA36 PHC; do
      MEA=$(cat $SUBJASHST2DIR/final/${idside}_corr_nogray_volumes.txt | grep "$sub" | awk '{printf("%d,%4.4f",$4,$5)}')
      if [[ $MEA == "" ]]; then
        MEA=","
      fi
      T2VOL="$T2VOL,$MEA"
    done
  done
  T2VOL=$(echo $T2VOL | sed 's/.\(.*\)/\1/')
  # compute mean
  ALLMEA=$T2VOL
  NL=18
  for ((i=1;i<=$NL;i++)); do
    LMEA=$(echo $ALLMEA | cut -d, -f $i)
    RMEA=$(echo $ALLMEA | cut -d, -f $((i+$NL)))
    if [[ $LMEA != "" && $RMEA != "" ]]; then
      MMEA=$(echo "scale=3;($LMEA+$RMEA)/2" | bc -l)
    else
      MMEA=""
    fi
    ALLMEA="$ALLMEA,$MMEA"
  done
  T2VOL=$ALLMEA

  # output
  echo "$ROW,$ICV,$T1VOL,$MSTVTTHK,$T2VOL" >> $OUTTMPDIR/$(printf %04d $ii)_${id}_${PREFIX}.csv 
  
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


