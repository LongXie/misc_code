#/bin/bash
#$ -S /bin/bash
set +e


ROOT=/home/longxie/ADNI2018/
SDROOT=/home/srdas/wd/ADNI23
INFOTXT=/home/srdas/wd/ADNI23/ashst2.txt
OUTDIR=$ROOT/ASHST2Seg_CongShan_ADNI

FILES=($(cat $INFOTXT))
N=${#FILES[*]}
#N=2

for ((i=0;i<${N};i++)); do

  file=${FILES[i]}
  id=$(echo $file | cut -d '/' -f 1)
  scandate=$(echo $file | cut -d '/' -f 2)
  ASHSDIR=$SDROOT/$(dirname $file)
  SUBJOUTDIR=$OUTDIR/$id/$scandate

  echo "copying $id $scandate"
  mkdir -p $SUBJOUTDIR

  if [[ ! -d $SUBJOUTDIR/bootstrap/fusion ]]; then
  # copy
  cp $ASHSDIR/final/${id}_*_lfseg_corr_nogray.nii.gz \
    $SUBJOUTDIR
  mkdir -p $SUBJOUTDIR/bootstrap/fusion
  cp $ASHSDIR/bootstrap/fusion/*nogray*.nii.gz \
    $SUBJOUTDIR/bootstrap/fusion

  fi

done

  



