#/bin/bash
#$ -S /bin/bash
set -e -x

#set -x -e

##############################################
# Setup environment
# Software PATH
LXROOT=/home/lxie
ANTSPATH=$LXROOT/pkg/bin/antsbin/bin
C3DPATH=$LXROOT/pkg/c3d_tool/bin
FSLROOT=$LXROOT/pkg/fsl
FSLPATH=$FSLROOT/bin
PAULROOT=/project/hippogang_2/pauly/
export ASHS_ROOT=/project/hippogang_2/longxie/pkg/ashs/ashs-fast
export PATH=$PATH:$ASHS_ROOT/bin
#MATLAB_BIN=/share/apps/matlab/R2017a/bin/matlab
export ALOHA_ROOT=/project/hippogang_2/longxie/pkg/aloha
SRPATH=$LXROOT/pkg/PatchSuperResolution/build_release

##############################################
# Directories
ROOT=$LXROOT/ADNI2018/
CODEDIR=$ROOT/scripts
#ASHSTEMPLATEDIR=/home/longxie/ASHS_T1/ASHSexp/exp008/headtailatlas/final
ANALYSISDIR=$ROOT/analysis_input
RAWDATADIR=$ROOT/rawdata
RAWDICOMDIR=/project/wolk_1/PUBLIC/dicom
RAWNIFTIDIR=/project/wolk_1/PUBLIC/nifti
DATADIR=$ROOT/dataset
DUMPDIR=$ROOT/dump_organize

# spreadsheets
OLDDATE="02202021"
#OLDT1T2TSV=$ROOT/info/T1Longitudinal_allTP/MRI3TListWithNIFTIPath_06132019.tsv
OLDT1T2TSV=$ROOT/RefreshT1T2NIFTI_${OLDDATE}/MRI3TListWithNIFTIPath_${OLDDATE}.tsv
NEWDATE="08022021"
OUTDIR=$ROOT/RefreshT1T2NIFTI_${NEWDATE}
LOGDIR=$OUTDIR/runlog
IMAGEUIDCSV=$ANALYSISDIR/MRILIST_T1T2_${NEWDATE}.csv
SCANNERFILE=/project/hippogang_1/srdas/wd/ADNI2/ADNI2_vs_ADNI3_MRI_Scanner_edited_withADNI2scannerchangeforT2.csv

#FULLLIST_TXT=$ANALYSISDIR/ADNI_T1_AllInfo_MTL_20181025.csv
#MRIMETA=$ANALYSISDIR/MRI3META_20181026.csv
#BLINFODIR=$ROOT/info/T1Baseline
#BLINFO_TXT=$BLINFODIR/checkdata_afterQC_04062018.csv
#BLNIFTIINFO_TXT=$BLINFODIR/subjectlist_04072018.csv
#LONGINFODIR=$ROOT/info/T1Longitudinal_allTP


# columns
RIDCol1="1"
IDCol1=$(csvcol.sh $IMAGEUIDCSV ID)
VISCODE2Col1=$(csvcol.sh $IMAGEUIDCSV VISCODE2)
SCANDATECol1=$(csvcol.sh $IMAGEUIDCSV SMARTDATE)
IMGIDT1Col1=$(csvcol.sh $IMAGEUIDCSV IMAGUID_T1)
IMGIDT2Col1=$(csvcol.sh $IMAGEUIDCSV IMAGUID_T2)
T1ISACCECol1=$(csvcol.sh $IMAGEUIDCSV T1ISACCE)
NT1NONEACCECol1=$(csvcol.sh $IMAGEUIDCSV NT1NONEACCE)
IMAGEUID_T1NONEACCECol1=$(csvcol.sh $IMAGEUIDCSV IMAGEUID_T1NONEACCE)
NT1ACCECol1=$(csvcol.sh $IMAGEUIDCSV NT1ACCE)
IMAGEUID_T1ACCECol1=$(csvcol.sh $IMAGEUIDCSV IMAGEUID_T1ACCE)
NT2Col1=$(csvcol.sh $IMAGEUIDCSV NT2)
IMAGEUID_T2ALLCol1=$(csvcol.sh $IMAGEUIDCSV IMAGEUID_T2ALL)

# old spreadsheet
FINALT1NIFTICol2=$(csvcol_tab.sh $OLDT1T2TSV FINALT1NIFTI)
FINALT2NIFTICol2=$(csvcol_tab.sh $OLDT1T2TSV FINALT2NIFTI)
IMGIDT1Col2=$(csvcol_tab.sh $OLDT1T2TSV IMAGEUID_T1)
T1ISACCECol2=$(csvcol_tab.sh $OLDT1T2TSV T1ISACCE)
IMGIDT2Col2=$(csvcol_tab.sh $OLDT1T2TSV IMAGEUID_T2)
BLSCANDATECol2=$(csvcol_tab.sh $OLDT1T2TSV BLSCANDATE)

# other
SCANSiteCol=$(csvcol.sh $SCANNERFILE SiteID)
Vendor3Col=$(csvcol.sh $SCANNERFILE MRI3)
Vendor2Col=$(csvcol.sh $SCANNERFILE MRI2)
Model3Col=$(csvcol.sh $SCANNERFILE Model3)
Model2Col=$(csvcol.sh $SCANNERFILE Model2)


#RIDCol2="1"
#IDCol2=$(csvcol_tab.sh $OLDT1T2TSV ID)
#SCANDATECol1=$(csvcol_tab.sh $OLDT1T2TSV SCANDATE)
#IMGIDT1Col2=$(csvcol_tab.sh $OLDT1T2TSV IMAGUID_T1)
#IMGIDT2Col2=$(csvcol_tab.sh $OLDT1T2TSV IMAGUID_T2)
#T1ISACCECol2=$(csvcol_tab.sh $OLDT1T2TSV T1ISACCE)
#NT1NONEACCECol2=$(csvcol_tab.sh $OLDT1T2TSV NT1NONEACCE)
#IMAGEUID_T1NONEACCECol2=$(csvcol_tab.sh $OLDT1T2TSV IMAGEUID_T1NONEACCE)
#NT1ACCECol2=$(csvcol_tab.sh $OLDT1T2TSV NT1ACCE)
#IMAGEUID_T1ACCECol2=$(csvcol_tab.sh $OLDT1T2TSV IMAGEUID_T1ACCE)
#NT2Col2=$(csvcol_tab.sh $OLDT1T2TSV NT2)
#IMAGEUID_T2ALLCol2=$(csvcol.sh $OLDT1T2TSV IMAGEUID_T2ALL)



#############################################
function main()
{
  reset_dir

  # refresh the T1T2 MRI spreadsheet
  #RefreshT1T2Spreadsheet

  # convert NIFTI
  ConvertALLNIFTI


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

#############################################
function RefreshT1T2Spreadsheet()
{
  # Load id number
  N=$(cat $IMAGEUIDCSV | wc -l)
  #N=100
  Nbegin=2
  #N=10
  OUTTMPDIR=$OUTDIR/tmp

  # link data and get orientation
  if [[ 1 == 0 ]]; then
  PREFIX=CD
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR
  mkdir -p $LOGDIR
  for ((i=${Nbegin};i<=${N};i++)); do

    pybatch.sh \
      -m "2G" \
      -n 1 \
      -N "${PREFIX}_${i}" \
      -o $DUMPDIR \
      $0 RefreshT1T2Spreadsheet_sub $i $OUTTMPDIR

    sleep 0.5

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
  fi

  # organize output
  header=$(cat -A $OLDT1T2TSV | head -n 1 | tail -n 1 | sed -e 's/\$//g' | sed -e 's/\^I/,/g' | cut -d , -f 1-7)
  echo "$header,IMAGEUID_T1,T1ISACCE,FINALT1NIFTI,IMAGEUID_T2,FINALT2NIFTI,NT1NONEACCE,IMAGEUID_T1NONEACCE,NT1ACCE,IMAGEUID_T1ACCE,NT2,IMAGEUID_T2ALL,Vendor2,Vendor3,Model2,Model3,BLSCANDATE" | sed 's/,/\t/g' > $OUTDIR/SearchScans.tsv
  cat $OUTTMPDIR/*.tsv >> $OUTDIR/SearchScans.tsv
}

function RefreshT1T2Spreadsheet_sub()
{
  i=$1
  OUTTMPDIR=$2

  # get rid and id
  ROW=$(cat -A $IMAGEUIDCSV | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g'  | sed "s/\[//g" | sed "s/\]//g" | sed -e "s/',\ '/;/g" | sed -e "s/;0//g" | sed -e "s/'//g"  | sed -e "s/\"//g")
  rid=$(echo $ROW | cut -f $RIDCol1 -d ",")
  id=$(echo $ROW | cut -f $IDCol1 -d ",")
  viscode2=$(echo $ROW | cut -f $VISCODE2Col1 -d ",")
  scandate=$(echo $ROW | cut -f $SCANDATECol1 -d ",")
  scandate=$(ReFormateDate $scandate)
  OUTROW=$(echo $ROW | cut -d , -f 1-7)
  SUBJDICOMDIR=$RAWDICOMDIR/$id
  SUBJNIFTIDIR=$RAWNIFTIDIR/$id

  # check whether id valide
  OUTTSV=$OUTTMPDIR/$(printf %05d $i)_${rid}.tsv
  if [[ $id == "" ]]; then
    set +e
    id=$(ls $RAWDICOMDIR | grep _S_$rid)
    set -e
  fi
  echo "$rid" >> $OUTTMPDIR/$(printf %05d $i)_${rid}.tsv
  if [[ $id == "" ]]; then
    exit
  fi

  # check if this scan has been processed before
  set +e
  line=$(cat -A $OLDT1T2TSV | grep $id | grep $viscode2 | grep $scandate | sed -e 's/\$//g' | sed -e 's/\^I/,/g') 
  OLDT1NIFTI=$(echo $line | cut -f $FINALT1NIFTICol2 -d ',')
  OLDT2NIFTI=$(echo $line | cut -f $FINALT2NIFTICol2 -d ',')

  ###############################################
  # check T1 first
  # if exist reuse, if not search for stuff
  if [[ $OLDT1NIFTI != "" ]]; then
    ImageIDT1=$(echo $line | cut -f $IMGIDT1Col2 -d ',')
    T1ISACCE=$(echo $line | cut -f $T1ISACCECol2 -d ',')
    T1NIFTI=$(echo $OLDT1NIFTI | sed  's/\/data\/tesla\-data\/PUBLIC\/ADNI2018/\/project\/wolk\_1\/PUBLIC/g')
  else
    T1NIFTI=""
    ImageIDT1=""
    CHANGET1=""
    T1ISACCE=""
    sep=$(mktemp -u  | awk -F "/" '{print $NF}')

    ImageIDT1=$(echo $ROW | cut -f $IMGIDT1Col1 -d ",")
    if [[ $ImageIDT1 != "" ]]; then
      T1DICOM=($(find $SUBJDICOMDIR/*/${scandate}*/*/*I${ImageIDT1}.dcm))
      T1ISACCE=$(echo $ROW | cut -f $T1ISACCECol1 -d ",")
      CHANGET1="0"

      # if can not find the one that is indicated, search for the other ones
      if [[ ${#T1DICOM[*]} == 0 ]]; then

        # first search none accelarated
        NT1NONEACCE=$(echo $ROW | cut -f $NT1NONEACCECol1 -d ",")
        if [[ $NT1NONEACCE -gt 1 ]]; then

          IMAGEUIDS=($(echo $ROW | cut -f $IMAGEUID_T1NONEACCECol1 -d "," | sed -e 's/\;/ /g'))

          for ((j=0;j<${#IMAGEUIDS[*]};j++)); do
            ImageIDT1TMP=${IMAGEUIDS[j]}
            if [[ $ImageIDT1TMP != $ImageIDT1 ]]; then
              T1DICOM=($(find $SUBJDICOMDIR/*/${scandate}*/*/*I${ImageIDT1TMP}.dcm))
              if [[ ${#T1DICOM[*]} != 0 ]]; then
                CHANGET1="1"
                ImageIDT1=$ImageIDT1TMP
                T1ISACCE="N"
                break
              fi
            fi
          done
        fi
      fi

      # search accelarated if necessary
      if [[ ${#T1DICOM[*]} == 0 ]]; then

        # could not find another one in the none acce, try to find one in the acce
        NT1ACCE=$(echo $ROW | cut -f $NT1ACCECol1 -d ",")
        if [[ $NT1ACCE -gt 0 ]]; then

          IMAGEUIDS=($(echo $ROW | cut -f $IMAGEUID_T1ACCECol1 -d "," | sed -e 's/\;/ /g' ))

          for ((j=0;j<${#IMAGEUIDS[*]};j++)); do
            ImageIDT1TMP=${IMAGEUIDS[j]}
            if [[ $ImageIDT1TMP != $ImageIDT1 ]]; then
              T1DICOM=($(find $SUBJDICOMDIR/*/${scandate}*/*/*I${ImageIDT1TMP}.dcm))
              if [[ ${#T1DICOM[*]} != 0 ]]; then
                CHANGET1="1"
                ImageIDT1=$ImageIDT1TMP
                T1ISACCE="Y"
                break
              fi
            fi
          done
        fi
      fi

      # get the NIFTI name
      if [[ ${#T1DICOM[*]} != 0 ]]; then
        dcm=$(c3d -dicom-series-list ${T1DICOM[0]} | sed -n '2,$p' | sed -e "s/\t/$sep/g")
        ser=$(echo $dcm | awk -F "$sep" '{print $1}')
        descr=$(echo $dcm | awk -F "$sep" '{print $4}')
        T1NIFTIDIR=$(dirname ${T1DICOM[0]} | sed 's/dicom/nifti/g')
        T1NIFTI=$(echo $T1NIFTIDIR/$( echo ${scandate}_ser$(printf %02d $ser)_${descr}.nii.gz | sed -e 's/ /_/g' | sed -e "s/\//_/g" ))
      else
        ImageIDT1=""
        CHANGET1=""
        T1ISACCE=""
        T1NIFTI=""
      fi

      # output status, if found new one, output to log
      if [[ $T1NIFTI != "" ]]; then
        msg="Found T1 scan for $id $scandate: $T1NIFTI"
        echo $msg
        echo $msg >> $LOGDIR/SearchNewScans.txt
      fi
    fi
  fi
  OUTROW="$OUTROW,$ImageIDT1,$T1ISACCE,$T1NIFTI"

  ###############################################
  # check T2
  # if exist reuse, if not search for stuff
  if [[ $OLDT2NIFTI != "" ]]; then
    ImageIDT2=$(echo $line | cut -f $IMGIDT2Col2 -d ',')
    T2NIFTI=$(echo $OLDT2NIFTI | sed  's/\/data\/tesla\-data\/PUBLIC\/ADNI2018/\/project\/wolk\_1\/PUBLIC/g')
  else
    T2NIFTI=""
    ImageIDT2=""
    CHANGET2=""

    ImageIDT2=$(echo $ROW | cut -f $IMGIDT2Col1 -d ",")
    if [[ $ImageIDT2 != "" ]]; then
      T2DICOM=($(find $SUBJDICOMDIR/*/${scandate}*/*/*I${ImageIDT2}.dcm))
      CHANGET2="0"

      # if can not find the one that is indicated, search for the other ones
      if [[ ${#T2DICOM[*]} == 0 ]]; then
        # first search none accelarated
        NT2=$(echo $ROW | cut -f $NT2Col1 -d ",")
        if [[ $NT2 -gt 1 ]]; then
          IMAGEUIDS=($(echo $ROW | cut -f $IMAGEUID_T2ALLCol1 -d "," | sed -e 's/\[//g' | sed -e 's/\]//g' | sed -e 's/,/ /g'))
          for ((j=0;j<${#IMAGEUIDS[*]};j++)); do
            ImageIDT2TMP=${IMAGEUIDS[j]}
            if [[ $ImageIDT2TMP != $ImageIDT2 ]]; then
              T2DICOM=($(find $SUBJDICOMDIR/*/${scandate}*/*/*I${ImageIDT2TMP}.dcm))
              if [[ ${#T2DICOM[*]} != 0 ]]; then
                CHANGET2="1"
                ImageIDT2=$ImageIDT2TMP
                break
              fi
            fi
          done
        fi
      fi
    fi

    if [[ ${#T2DICOM[*]} != 0 ]]; then
      dcm=$(c3d -dicom-series-list ${T2DICOM[0]} | sed -n '2,$p' | sed -e "s/\t/$sep/g")
      ser=$(echo $dcm | awk -F "$sep" '{print $1}')
      descr=$(echo $dcm | awk -F "$sep" '{print $4}')
      T2NIFTIDIR=$(dirname ${T2DICOM[0]} | sed 's/dicom/nifti/g')
      T2NIFTI=$(echo $T2NIFTIDIR/$( echo ${scandate}_ser$(printf %02d $ser)_${descr}.nii.gz | sed -e 's/ /_/g' | sed -e "s/\//_/g" ))
    else
      CHANGET2=""
      ImageIDT2=""
    fi
  fi
  OUTROW="$OUTROW,$ImageIDT2,$T2NIFTI"
  set -e

  # get the rest info
  NT1NONEACCE=$(echo $ROW | cut -f $NT1NONEACCECol1 -d ",")
  IMAGEUID_T1NONEACCE=$(echo $ROW | cut -f $IMAGEUID_T1NONEACCECol1 -d ",")
  NT1ACCE=$(echo $ROW | cut -f $NT1ACCECol1 -d ",")
  IMAGEUID_T1ACCE=$(echo $ROW | cut -f $IMAGEUID_T1ACCECol1 -d ",")
  NT2=$(echo $ROW | cut -f $NT2Col1 -d ",")
  IMAGEUID_T2ALL=$(echo $ROW | cut -f $IMAGEUID_T2ALLCol1 -d ",")
  OUTROW="$OUTROW,$NT1NONEACCE,$IMAGEUID_T1NONEACCE,$NT1ACCE,$IMAGEUID_T1ACCE,$NT2,$IMAGEUID_T2ALL"

  # get site, phas, viscode, viscode2, model2, model3, vendor2, vendor3
  Site=$(echo $id | cut -f 1 -d "_")
  searchSite=$(echo $Site | sed -e "s/^0//g" | sed -e "s/^0//g")
  set +e
  Scanline=$(cat $SCANNERFILE | grep ^${searchSite}, )
  set -e
  if [ "$Scanline" != "" ]; then
    Model2=$(echo $Scanline | cut -f $Model2Col -d "," | sed -e 's/ //g')
    Model3=$(echo $Scanline | cut -f $Model3Col -d "," | sed -e 's/ //g')
    Vendor2=$(echo $Scanline | cut -f $Vendor2Col -d "," | sed -e 's/ //g')
    Vendor3=$(echo $Scanline | cut -f $Vendor3Col -d "," | sed -e 's/ //g')
  fi
  OUTROW="$OUTROW,$Model2,$Model3,$Vendor2,$Vendor3"

  # determine basline scan date
  set +e
  BLScanDate=$(echo $line | cut -f $BLSCANDATECol2 -d ',')
  if [[ $BLScanDate == "" ]]; then
    # search for the earliest scan date
    BLScanDate=$(cat -A $IMAGEUIDCSV | grep $id | head -n 1 | cut -f $SCANDATECol1 -d ",")
  fi
  BLScanDate=$(ReFormateDate $BLScanDate)
  set -e
  OUTROW="$OUTROW,$BLScanDate"

  # output
  echo "$OUTROW" | sed 's/,/\t/g' > $OUTTMPDIR/$(printf %05d $i)_${rid}.tsv  
}

#############################################
function ConvertALLNIFTI()
{
  # Load id number
  TSV=$OUTDIR/SearchScans.tsv  
  N=$(cat $TSV | wc -l)
  Nbegin=2
  #N=20
  OUTTMPDIR=$OUTDIR/tmp

  # link data and get orientation
  if [[ 1 == 1 ]]; then
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR
  mkdir -p $LOGDIR

  PREFIX=CAN
  for ((i=${Nbegin};i<=${N};i++)); do

    pybatch.sh \
      -m "2G" \
      -n 1 \
      -N "${PREFIX}_${i}" \
      -o $DUMPDIR \
      $0 ConvertALLNIFTI_sub $i $OUTTMPDIR
    sleep 0.5

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

  fi

  # organize output
  echo "RID,ID,T1NIFTI,T1SUCCESS,T2NIFTI,T2SUCCESS" > $OUTDIR/log_convert_nifti.csv
  cat $OUTTMPDIR/*.csv >> $OUTDIR/log_convert_nifti.csv

}

function ConvertALLNIFTI_sub()
{
  i=$1
  OUTTMPDIR=$2
  TSV=$OUTDIR/SearchScans.tsv
  OUTCSV=$OUTTMPDIR/$(printf %05d $i).csv
  

  # find column
  #RIDCol1=$(csvcol.sh $IMAGEUIDCSV RID)
  RIDCol1="1"
  IDCol1=$(csvcol_tab.sh $TSV ID)
  SCANDATECol1=$(csvcol_tab.sh $TSV SCANDATE)
  IMGIDT1Col1=$(csvcol_tab.sh $TSV IMAGEUID_T1)
  IMGIDT2Col1=$(csvcol_tab.sh $TSV IMAGEUID_T2)
  T1NIFTICol1=$(csvcol_tab.sh $TSV FINALT1NIFTI)
  T2NIFTICol1=$(csvcol_tab.sh $TSV FINALT2NIFTI)

  # get ROWs
  ROW=$(cat -A $TSV | head -n $1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g' | sed -e 's/\^I/,/g')
  id=$(echo $ROW | cut -f $IDCol1 -d ",")
  rid=$(echo $ROW | cut -f $RIDCol1 -d ",")

  # output basic information
  OUTCSV=$OUTTMPDIR/$(printf %05d $i)_${rid}.csv
  echo "$rid,$id,,,," > $OUTCSV

  # get other info
  scandate=$(echo $ROW | cut -f $SCANDATECol1 -d ",")
  scandate=$(ReFormateDate $scandate)
  IMAGEIDT1=$(echo $ROW | cut -f $IMGIDT1Col1 -d ",")
  IMAGEIDT2=$(echo $ROW | cut -f $IMGIDT2Col1 -d ",")
  T1NIFTI=$(echo $ROW | cut -f $T1NIFTICol1 -d ",")
  T2NIFTI=$(echo $ROW | cut -f $T2NIFTICol1 -d ",")
  SUBJDICOMDIR=$RAWDICOMDIR/$id
  SUBJNIFTIDIR=$RAWNIFTIDIR/$id
  sep=$(mktemp -u  | awk -F "/" '{print $NF}')

  set +e
  # convert T1
  T1STATUS=""
  if [[ $T1NIFTI != "" ]]; then
  if [[ ! -f $T1NIFTI ]]; then

    # get DICOM and convert to NIFTI
    T1DICOM=($(find $SUBJDICOMDIR/*/${scandate}*/*/*I${IMAGEIDT1}.dcm))
    if [[ ${#T1DICOM[*]} != 0 ]]; then
      dcm=$(c3d -dicom-series-list ${T1DICOM[0]} | sed -n '2,$p' | sed -e "s/\t/$sep/g")
      ser=$(echo $dcm | awk -F "$sep" '{print $1}')
      serid=$(echo $dcm | awk -F "$sep" '{print $NF}');
      T1NIFTIDIR=$(dirname ${T1DICOM[0]} | sed 's/dicom/nifti/g')
      mkdir -p $T1NIFTIDIR
      c3d -dicom-series-read "${T1DICOM[0]}" "${serid}" \
        -o $T1NIFTI

      if [[ -f $T1NIFTI ]]; then
        T1STATUS="CONVERED"
      else
        T1STATUS="FAILED"
      fi
    else
      T1STATUS="DICOMNOTEXIST"
    fi
  else
    T1STATUS="EXIST"
  fi
  fi

  # convert T2
  T2STATUS=""
  if [[ $T2NIFTI != "" ]]; then
  if [[ ! -f $T2NIFTI ]]; then

    # get DICOM and convert to NIFTI
    T2DICOM=($(find $SUBJDICOMDIR/*/${scandate}*/*/*I${IMAGEIDT2}.dcm))
    if [[ ${#T2DICOM[*]} != 0 ]]; then
      dcm=$(c3d -dicom-series-list ${T2DICOM[0]} | sed -n '2,$p' | sed -e "s/\t/$sep/g")
      ser=$(echo $dcm | awk -F "$sep" '{print $1}')
      serid=$(echo $dcm | awk -F "$sep" '{print $NF}');
      T2NIFTIDIR=$(dirname ${T2DICOM[0]} | sed 's/dicom/nifti/g')
      mkdir -p $T2NIFTIDIR
      c3d -dicom-series-read "${T2DICOM[0]}" "${serid}" \
        -o $T2NIFTI

      if [[ -f $T2NIFTI ]]; then
        T2STATUS="CONVERED"
      else
        T2STATUS="FAILED"
      fi
    else
      T2STATUS="DICOMNOTEXIST"
    fi
  else
    T2STATUS="EXIST"
  fi
  fi
  set -e
  
  # output
  echo "$rid,$id,$T1NIFTI,$T1STATUS,$T2NIFTI,$T2STATUS" > $OUTCSV
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

  cmd=$1
  shift
  $cmd $@

fi
