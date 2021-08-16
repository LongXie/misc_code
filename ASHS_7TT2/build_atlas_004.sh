#!/bin/bash
#$ -S /bin/bash
set -e
#set -x -e

##############################################################################
# INIT
#ANTSPATH=/home/avants/bin/ants
ANTSPATH=/data/picsl/longxie/pkg/bin/antsbin/bin/
#ANTSPATH=/data/picsl/pauly/bin/ants/
NEWANTSDIR=/share/apps/ANTs/2014-06-23/build/bin/
C3DPATH=/data/picsl/longxie/pkg/c3d_tool/bin
PKGDIR=/data/picsl/longxie/pkg/
MATLAB_BIN=/share/apps/matlab/R2017a/bin/matlab
GREEDY_BIN=/data/picsl/pauly/bin/
export PATH=$ANTSPATH:$C3DPATH:$PATH
RCODEDIR=/home/longxie/ASHS_T1/scripts
#export ASHS_ROOT=/data/picsl/pauly/wolk/ashs-fast
export ASHS_ROOT=/data/picsl/longxie/pkg/ashs/ashs-fast

# TMPDIR
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp
fi

runid=4

ROOT=/home/longxie/ASHS_7TT2
WORKDIR=$ROOT/work
RAWDATADIR=$ROOT/rawdata/
RAWSEGDIR=$RAWDATADIR/ABC_7T_Atlas_final
RAWATLASDIR=$RAWDATADIR/atlas_withflip
ANALYSISDIR=$ROOT/analysis_input
CODEDIR=$ROOT/scripts
RUNDIR=$WORKDIR/run_$(printf %03d $runid)
EVALDIR=$RUNDIR/evaluation
DUMPDIR=$WORKDIR/dump

ATLASCSV=$ANALYSISDIR/7TAtlasIDs.csv




##############################################################################
function main()
{
  reset_dir

  # find and flip the scans, create manifest file
  #PrepareAtlas

  # create atlas and cross validation
  #CreateAtlas
 
  # compute cross-validation Dice score
  #Evaluate

  # test the pipeline using David Berron's data
  #TestDavid
  #TestDavid3TT2

  # evaluate atlas overlap
  EvaluateDavidDifferentASHSMethods

}

##########################################################
function PrepareAtlas()
{
  N=$(cat $ATLASCSV | wc -l)

  # initialize tmp dir
  OUTTMPDIR=$WORKDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR
  mkdir -p $RAWATLASDIR

  # submit jobs to find scans and flip them
  PREFIX=PA
  for ((i=2;i<=${N};i++)); do

    qsub -cwd -o $DUMPDIR -j y \
         -q all.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX}_${i}" \
         $0 PrepareAtlas_sub $i $OUTTMPDIR
    sleep 0.03

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # get the manifest file
  cat $OUTTMPDIR/*.txt > $RAWATLASDIR/manifest_flip.txt
}

function PrepareAtlas_sub()
{
  i=$1
  OUTTMPDIR=$2

  # get scan info
  aid=$(cat -A $ATLASCSV | head -n $i | tail -n 1 | cut -d, -f 1)
  id=$(echo $aid | cut -d _ -f 1)
  scandate=$(echo $aid | cut -d _ -f 2)

  # find T1, T2 and segmentations
  SCANDATADIR=/data/jux/wolk_group/Prisma3T/relong
  T1=$(ls $SCANDATADIR/$id/MRI3T/*/processed/*mprage_trim.nii.gz | head -n 1)
  #T1=$(ls $SCANDATADIR/$id/MRI7T/$scandate/processed/T1_*.nii.gz)
  T2=$(ls $SCANDATADIR/$id/MRI7T/$scandate/processed/T2_*.nii.gz)
  LSEG=$(ls $RAWSEGDIR/$aid/${id}_left_LW5.nii.gz)
  RSEG=$(ls $RAWSEGDIR/$aid/${id}_right_LW5.nii.gz) 

  # copy files
  mkdir -p $RAWATLASDIR/orig_$aid
  cp $T1 $RAWATLASDIR/orig_$aid/orig_T1.nii.gz
  cp $T2 $RAWATLASDIR/orig_$aid/orig_T2.nii.gz
  cp $LSEG $RAWATLASDIR/orig_$aid/orig_left_LW5.nii.gz
  cp $RSEG $RAWATLASDIR/orig_$aid/orig_right_LW5.nii.gz

  # output a line to manifest file
  echo "orig_${aid} $RAWATLASDIR/orig_$aid/orig_T1.nii.gz $RAWATLASDIR/orig_$aid/orig_T2.nii.gz $RAWATLASDIR/orig_$aid/orig_left_LW5.nii.gz $RAWATLASDIR/orig_$aid/orig_right_LW5.nii.gz" > \
    $OUTTMPDIR/$(printf %04d $i)_${aid}.txt
 
  # flip the scans
  mkdir -p $RAWATLASDIR/flip_$aid
  /data/picsl/pauly/bin/c3d_affine_tool -sform $T2 -tran 447 0 0 -scale -1 1 1 \
    -sform $T2 -inv -mult -mult -mult \
    -o $RAWATLASDIR/flip_$aid/flip_${aid}.mat
  c3d $T1 -pad 30x30x30vox 30x30x30vox 0 \
    -dup -reslice-matrix $RAWATLASDIR/flip_$aid/flip_${aid}.mat \
    -as FT1 \
    -as FT1 \
    -thresh 0.0001 inf 1 0 -trim 0x0x0vox \
    -push FT1 -int 0 -reslice-identity \
    -o $RAWATLASDIR/flip_$aid/flip_T1.nii.gz
  c3d $T2 -dup -reslice-matrix $RAWATLASDIR/flip_$aid/flip_${aid}.mat \
    -o $RAWATLASDIR/flip_$aid/flip_T2.nii.gz
  c3d $LSEG -dup -int 0 -reslice-matrix $RAWATLASDIR/flip_$aid/flip_${aid}.mat \
    -o $RAWATLASDIR/flip_$aid/flip_right_LW5.nii.gz
  c3d $RSEG -dup -int 0 -reslice-matrix $RAWATLASDIR/flip_$aid/flip_${aid}.mat \
    -o $RAWATLASDIR/flip_$aid/flip_left_LW5.nii.gz
 
  # output the flip info to manifest file
  echo "flip_${aid} $RAWATLASDIR/flip_$aid/flip_T1.nii.gz $RAWATLASDIR/flip_$aid/flip_T2.nii.gz $RAWATLASDIR/flip_$aid/flip_left_LW5.nii.gz $RAWATLASDIR/flip_$aid/flip_right_LW5.nii.gz" >> $OUTTMPDIR/$(printf %04d $i)_${aid}.txt 
}

##########################################################
function CreateAtlas()
{
  # make run directory
  mkdir -p $RUNDIR

  # make xval_flip.txt
  #cat $RAWATLASDIR/manifest_flip.txt | awk '{print $1}' > $RUNDIR/xval_flip.txt
    

  # run ASHS train
  $ASHS_ROOT/bin/ashs_train.sh \
    -D $RAWATLASDIR/manifest_flip.txt \
    -L $ANALYSISDIR/7TAtlas_snaplabel.txt \
    -w $RUNDIR \
    -C $ANALYSISDIR/ashs_user_config_7T_3TT1.sh \
    -x $WORKDIR/xval_flip.txt \
    -z $ANALYSISDIR/ashs-fast-train-z-7T.sh \
    -s 1-7 \
    | tee -a $RUNDIR/ashs_train.log
}

##########################################################
function Evaluate()
{
  N=$(cat $ATLASCSV | wc -l)

  # initialize tmp dir
  OUTTMPDIR=$EVALDIR/tmp
  rm -rf $OUTTMPDIR
  mkdir -p $OUTTMPDIR
  mkdir -p $RAWATLASDIR

  # submit jobs to find scans and flip them
  PREFIX=EV
  for ((i=2;i<=${N};i++)); do

    qsub -cwd -o $DUMPDIR -j y \
         -q all.q \
         -l h_vmem=4.1G,s_vmem=4G \
         -N "${PREFIX}_${i}" \
         $0 Evaluate_sub $i $OUTTMPDIR
    sleep 0.03

  done

  # Wait for completion
  qsub -cwd -o $DUMPDIR -j y \
       -q all.q \
       -hold_jid "${PREFIX}_*" -sync y -b y \
       sleep 1

  # get the manifest file
  for type in heur corr_nogray corr_usegray; do
    echo "ID,SIDE,TYPE,CA1,CA2,DG,CA3,Tail,SUB,ERC,BA35,BA36,PHC,MISC,Cysts,HS,CS,END" \
      > $EVALDIR/overlap_${type}.csv
    cat $OUTTMPDIR/*_${type}.csv >> $EVALDIR/overlap_${type}.csv
  done
}

function Evaluate_sub()
{
  i=$1
  OUTTMPDIR=$2

  # get scan info
  MANIFESTFN=$RAWATLASDIR/manifest_flip.txt
  ROW=$(cat $MANIFESTFN | head -n $i | tail -n 1)
  id=$(echo $ROW | cut -d " " -f 1)
  
  # xval dir
  XVALDIR=$RUNDIR/xval
  for type in heur corr_nogray corr_usegray; do
    for side in left right; do

      outstr="${id},${side},${type}"
      OVALFN=$(find $XVALDIR/*/test/*${id}/final/*${id}_${side}_${type}_overlap.txt)
      tmpstr=($(cat $OVALFN | awk '{print $4}'))
      tmpstr=$(echo ${tmpstr[*]} | sed 's/ /,/g')
      outstr="$outstr,$tmpstr,-1"
      echo $outstr >> $OUTTMPDIR/$(printf %04d $i)_${id}_${type}.csv

    done
  done
}

##########################################################
function TestDavid()
{
  ATLAS=$RUNDIR/final
  DAVIDDATADIR=$RUNDIR/testDavid

  for id in $(ls $DAVIDDATADIR ); do

    echo "Processing $id"

    T1=$DAVIDDATADIR/$id/${id}_t1_orig.nii
    T2=$DAVIDDATADIR/$id/${id}_t2_orig.nii
    ASHSRUNDIR=$DAVIDDATADIR/$id/ashs_${id}_ashsfast
    
    $ASHS_ROOT/bin/ashs_main.sh \
        -a $ATLAS -d -T -I $id -g $T1 -f $T2 \
        -G \
        -Q -s 1-7 \
        -z $CODEDIR/ashs-fast-z.sh \
        -w $ASHSRUNDIR &
      sleep 0.1
    


  done

} 

##########################################################
function TestDavid3TT2()
{
  #ATLAS=$RUNDIR/final
  ATLAS=/data/picsl/pauly/wolk/abc_prisma/exp03/final
  DAVIDDATADIR=$RUNDIR/testDavid

  for id in $(ls $DAVIDDATADIR ); do

    echo "Processing $id"

    T1=$DAVIDDATADIR/$id/${id}_t1_orig.nii
    T2=$DAVIDDATADIR/$id/${id}_t2_orig.nii
    ASHSRUNDIR=$DAVIDDATADIR/$id/ashs_${id}_ashsfast_3TT2_nomask

    $ASHS_ROOT/bin/ashs_main.sh \
        -a $ATLAS -d -T -I $id -g $T1 -f $T2 \
        -Q -s 1-7 \
        -z $CODEDIR/ashs-fast-z.sh \
        -w $ASHSRUNDIR &
      sleep 0.1

#        -G \


  done

}

##########################################################
function EvaluateDavidDifferentASHSMethods()
{
  DAVIDDATADIR=$RUNDIR/testDavid

  for id in $(ls $DAVIDDATADIR ); do

    echo "Processing $id"

    T1=$DAVIDDATADIR/$id/${id}_t1_orig.nii
    T2=$DAVIDDATADIR/$id/${id}_t2_orig.nii

    for type in ashsfast ashsfast_3TT2 ashsfast_3TT2_nomask; do

      ASHSRUNDIR=$DAVIDDATADIR/$id/ashs_${id}_${type}
      BSDIR=$ASHSRUNDIR/bootstrap
      OUTTXT=$BSDIR/overlap_atlas.csv
      rm -f $OUTTXT

      for side in left right; do
        for atlas in $(ls $BSDIR | grep tseg_${side}_train); do

          OVL=$(c3d $BSDIR/fusion/lfseg_corr_nogray_${side}.nii.gz \
              $BSDIR/$atlas/atlas_to_native_segvote.nii.gz \
              -foreach -replace 6 0 7 0 13 0 14 0 15 0 16 0 -endfor \
              -label-overlap \
              | awk '{print $3}' | awk '{printf("%s ", $1)}' | awk '{print $3}')
          echo ${id},${type},${side},${atlas},${OVL}
          echo "${id},${type},${side},${atlas},${OVL}" >> $OUTTXT

        done
      done
    done
  done

}

##########################################################
function reset_dir()
{
  rm -rf $DUMPDIR
  mkdir -p $DUMPDIR
}

##########################################################
if [[ $# -lt 2 ]]; then

  main

else

  cmd=$1
  shift
  $cmd $@

fi







