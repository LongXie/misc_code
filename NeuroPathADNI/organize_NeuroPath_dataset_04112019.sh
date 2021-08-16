#/bin/bash
#$ -S /bin/bash
set -e 
set -x -e

##############################################
# Setup environment
# Software PATH
ANTSPATH=/data/picsl/longxie/pkg/antsbin/bin
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

##############################################
# Directories
ROOT=/home/longxie/NeuroPathADNI/
CODEDIR=$ROOT/scripts
ANALYSISDIR=$ROOT/analysis_input
FULLLIST_TXT=$ANALYSISDIR/NEUROPATH_04_12_18.csv
SCANLIST_TXT=$ANALYSISDIR/NeuroPath_04_12_18_4_05_2019.csv
RAWDATADIR=/data/jux/PUBLIC/ADNI2018/NeuroPath_dataset
RAWDICOMDIR=$RAWDATADIR/dicom
RAWNIFTIDIR=$RAWDATADIR/nifti
DATADIR=$ROOT/dataset
DUMPDIR=$ROOT/dump_organize
OUTORGANIZEDIR=$ROOT/output_organize

# find all the column numbers
dateCol=$(csvcol.sh $MRIMETA EXAMDATE)
RIDCol=$(csvcol.sh $MRIMETA RID)
SiteCol=$(csvcol.sh $MRIMETA SITEID)
PhaseCol=$(csvcol.sh $MRIMETA PHASE)

#############################################
function main()
{
  reset_dir

  # check whether data exist
  #FilterScanInfo

  # convert to nifti
  ConvertNIFTI


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
function FilterScanInfo()
{
  # Load id number
  N=$(cat $FULLLIST_TXT | wc -l)
  #N=100
  rm -rf $OUTORGANIZEDIR/tmp
  INFOTMPDIR=$OUTORGANIZEDIR/tmp
  mkdir -p $INFOTMPDIR

  # link data and get orientation
  PREFIX=FSI
  for ((i=2;i<=${N};i++)); do

    qsub -cwd -o $DUMPDIR -j y \
         -q all.q,basic.q \
         -l h_vmem=2.1G,s_vmem=2G \
         -N "${PREFIX}_${i}" \
         $0 FilterScanInfo_sub $i $INFOTMPDIR
    sleep 0.1

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # organize output
  echo "RID,ID,SCANDATE,Group,Sex,Age,Sequence,ImageUID,NUniqueScandate,SelectedScans,BLSCANDATE,DATEDIFFTOBL" > $OUTORGANIZEDIR/filterdata.csv
  cat $INFOTMPDIR/*.csv >> $OUTORGANIZEDIR/filterdata.csv

}

function FilterScanInfo_sub()
{
  i=$1
  longi=$(printf %03d $i)
  INFOTMPDIR=$2
  rm -f $INFOTMPDIR/${i}_${rid}.csv

  # get rid and id
  ROW=$(cat $FULLLIST_TXT | head -n $i | tail -n 1)
  rid=$(echo $ROW | awk -F "," '{print $1}')
  longrid=$(printf %04d $rid)
  id=$(ls $RAWDICOMDIR | grep $longrid)

  if [[ $id == "" ]]; then
    echo "$rid" >> $INFOTMPDIR/${i}_${rid}.csv
    continue
  fi

  # summarize information
  SCANDATES=($(cat $SCANLIST_TXT | grep $id | grep -E 'RAGE|SPGR' \
    | while read line; do
      scandate=$(echo $line | awk -F "," '{print $10}' | sed -e 's/"//g')
      scandate=$(ReFormateDate $scandate)
      echo $scandate
    done))

  UNIQSCANDATES=($(echo ${SCANDATES[@]} | tr ' ' '\n' | sort -u | tr '\n' ' '))
  NUNIQUESCANDATES=${#UNIQSCANDATES[*]}
  blscandate=${UNIQSCANDATES[0]}

  #echo ${UNIQSCANDATES[*]}

  # select one scan for each scan date (give preference to non accelarated scans)
  FINALSELIDS=""
  MISMATCHSTR=("SAGITTAL|GRAPPA|grappa|Acce|acce|_P2_|SENSE" "wefwefewfwefwefe")
  for ((ii=0;ii<$NUNIQUESCANDATES;ii++)); do
  for ((jj=0;jj<${#MISMATCHSTR[*]};jj++)); do
    SELIDS=($(cat $SCANLIST_TXT | grep $id | grep -E 'RAGE|SPGR' | grep -vE "${MISMATCHSTR[jj]}" \
      | while read line; do
        scandate=$(echo $line | awk -F "," '{print $10}' | sed -e 's/"//g')
        scandate=$(ReFormateDate $scandate)
        if [[ $scandate == ${UNIQSCANDATES[ii]} ]]; then
          imageid=$(echo $line | awk -F "," '{print $1}' | sed -e 's/"//g')
          echo $imageid
        fi
      done))

    #echo "$id, $jj, ${SELIDS[*]}"

    if [[ $SELIDS == "" ]]; then
      echo "$id $scandate is looking for accelarated scans"
    else
      SORTEDIDS=($(echo ${SELIDS[@]} | tr ' ' '\n' | sort -r | tr '\n' ' '))
      FINALSELIDS[$ii]=${SORTEDIDS[0]}
      break
    fi

  done

  if [[ $SELIDS == "" ]]; then
    echo "No image id is found for $id $scandate"
    echo "No image id is found for $id $scandate" >> $INFOTMPDIR/${i}_${rid}.csv
    exit
  fi

  done

  #echo ${FINALSELIDS[*]}

  # check the dicom info file for scan infos
  if [[ 1 == 1 ]]; then
  cat $SCANLIST_TXT | grep $id | grep -E 'RAGE|SPGR' \
    | while read line; do
      scandate=$(echo $line | awk -F "," '{print $10}' | sed -e 's/"//g')
      scandate=$(ReFormateDate $scandate)
      demog=$(echo $line | awk -F "," '{print $3","$4","$5}' | sed -e 's/"//g')
      seq=$(echo $line | awk -F "," '{print $8}')
      imageid=$(echo $line | awk -F "," '{print $1}' | sed -e 's/"//g')
      date_diff=$(( \
        ($(date -d $scandate +%s) - \
        $(date -d $blscandate +%s) )/(60*60*24) ))
      # check if this scan is selected
      selected=""
      for ((ii=0;ii<${#FINALSELIDS[*]};ii++)); do
        if [[ $imageid == ${FINALSELIDS[ii]} ]]; then
          selected="selected"
          break
        fi
      done
      echo "$rid,$id,$scandate,$demog,$seq,$imageid,$NUNIQUESCANDATES,$selected,$blscandate,$date_diff" >> $INFOTMPDIR/${longi}_${rid}.csv
    done
  fi
}


#############################################
function ConvertNIFTI()
{
  # Load id number
  ALLINFO=$OUTORGANIZEDIR/filterdata.csv
  N=$(cat $ALLINFO | wc -l)
  #N=3

  rm -rf $OUTORGANIZEDIR/tmp_nifti
  INFOTMPDIR=$OUTORGANIZEDIR/tmp_nifti
  mkdir -p $INFOTMPDIR
  #rm -rf $OUTORGANIZEDIR/log.csv
  
  #echo "RID,ID,ScanDate,NNDDirs,NDDirsNames,NDDateFolders" > $LONGINFODIR/log.csv

  # link data and get orientation
  for ((i=2;i<=${N};i++)); do
    
    # check if selected
    ROW=$(cat $ALLINFO | head -n $i | tail -n 1)
    selected=$(echo $ROW | awk -F "," '{print $10}')

    # submit jobs to convert to NIFTI
    if [[ $selected == "selected" ]]; then
      PREFIX=CN
      qsub -cwd -o $DUMPDIR -j y \
           -q all.q,basic.q,gpu.q \
           -l h_vmem=4.1G,s_vmem=4G \
           -N "${PREFIX}_${i}" \
           $0 ConvertNIFTI_sub $i $INFOTMPDIR
      sleep 0.1
    fi

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q,basic.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # organize output
  header=$(cat -A $ALLINFO | head -n 1 | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
  echo "$header,FINALT1NIFTI" > $OUTORGANIZEDIR/checkdata.csv
  cat $INFOTMPDIR/*.csv >> $OUTORGANIZEDIR/checkdata.csv
}

function ConvertNIFTI_sub()
{
  i=$1
  longi=$(printf %03d $i)
  INFOTMPDIR=$2

  # get columns
  ALLINFO=$OUTORGANIZEDIR/filterdata.csv
  RIDCol=1
  IDCol=$(csvcol.sh $ALLINFO ID)
  ScanDateCol=$(csvcol.sh $ALLINFO SCANDATE)
  IMGUIDCol=$(csvcol.sh $ALLINFO ImageUID)

  # get the information
  ROW=$(cat -A $ALLINFO | head -n $i | tail -n 1 | sed -e 's/\^M//g' | sed -e 's/\$//g')
  rid=$(echo $ROW | cut -f $RIDCol -d ",")
  id=$(echo $ROW | cut -f $IDCol -d ",")
  scandate=$(echo $ROW | cut -f $ScanDateCol -d ",")
  DD=$(date -d "$scandate" '+%d')
  MM=$(date -d "$scandate" '+%m')
  YYYY=$(date -d "$scandate" '+%Y')
  PREFIX1="${YYYY}-${MM}-${DD}"
  scandate=$(ReFormateDate $scandate)
  IMGUID=$(echo $ROW | cut -f $IMGUIDCol -d ",")
  echo "$ROW,None" > $INFOTMPDIR/${longi}_${id}_T1.csv

  # T1 scan
  sep=$(mktemp -u  | awk -F "/" '{print $NF}')
  SUBJDICOMDIR=$RAWDICOMDIR/$id
  SUBJNIFTIDIR=$RAWNIFTIDIR/$id/$dicomfolder
  T1DICOM=($(find $SUBJDICOMDIR/*/${scandate}*/*/*I${IMGUID}.dcm))
  if [[ $T1DICOM == "" ]]; then
    T1NIFTI=""
  else
    line=$(c3d -dicom-series-list ${T1DICOM[0]} | sed -n '2,$p' | sed -e "s/\t/$sep/g")
    ser=$(echo $line | awk -F "$sep" '{print $1}')
    descr=$(echo $line | awk -F "$sep" '{print $4}')
    #serid=$(echo $line | awk -F "$sep" '{print $NF}');
    T1NIFTIDIR=$(dirname ${T1DICOM[0]} | sed 's/dicom/nifti/g')
    T1NIFTI=$(echo $T1NIFTIDIR/$( echo ${scandate}_ser$(printf %02d $ser)_${descr}.nii.gz | sed -e 's/ /_/g' | sed -e "s/\//_/g" ))
    mkdir -p $T1NIFTIDIR
    c3d -dicom-series-read "${T1DICOM[0]}" "${serid}" \
      -o $T1NIFTI
    if [[ ! -f $T1NIFTI ]]; then
      T1NIFTI=""
    fi
  fi

  echo "$ROW,$T1NIFTI" > $INFOTMPDIR/${longi}_${id}_T1.csv


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

