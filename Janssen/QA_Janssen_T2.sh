#/bin/bash
#$ -S /bin/bash
set -e


FILE=$1
STARTLINE=$2
if [[ $STARTLINE == "" ]]; then
  STARTLINE=1
fi

TMPDIR=~/.tmpdir

Nline=$(cat $FILE | wc -l)

for ((i=$STARTLINE;i<=$Nline;i++)); do

line=$(cat $FILE | head -n $i | tail -n 1)
ID=$(echo $line | cut -f1 -d,)
if [[ $ID == "ID" ]]; then
  continue
fi
#DATE=$(echo $line | cut -f2 -d,)

SUBDIR=/data/jux/wolk_group/janssen/CPSS_MRI_vols_UPenn/UPenn/$ID/
ASHSDIR=$SUBDIR/ASHST2/


#echo $line
echo "Line $i: $ID $DATE"

for side in left right; do
echo "    ${side} side"

  mkdir -p $TMPDIR

#if [[ ! -f $ASHSDIR/final/${ID}_${side}_lfseg_corr_usegray_native_chunk.nii.gz ]]; then
  /data/picsl/pauly/bin/c3d $ASHSDIR/tse_native_chunk_${side}.nii.gz \
    $ASHSDIR/final/${ID}_${side}_lfseg_corr_nogray_CortexCleaned.nii.gz \
    -int 0 -reslice-identity \
    -o $TMPDIR/${ID}_${side}_lfseg_corr_nogray_CortexCleaned_native_chunk.nii.gz
#fi

/share/apps/itksnap/itksnap-3.8.0-beta-20181019-Linux-x86_64-qt4/bin/itksnap -g $ASHSDIR/tse_native_chunk_${side}.nii.gz \
        -l /data/picsl/pauly/wolk/abc_prisma/exp03/final/snap/snaplabels.txt \
        -s $TMPDIR/${ID}_${side}_lfseg_corr_nogray_CortexCleaned_native_chunk.nii.gz

    rm -rf $TMPDIR

done

echo ""

done
