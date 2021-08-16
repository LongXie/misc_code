#/bin/bash
#$ -S /bin/bash

ASHS_WORK=$1
STAGE=$2

case $STAGE in

  1)
  echo " -q all.q -l h_vmem=12.1G,s_vmem=12G";;

  2)
  echo " -q all.q -l h_vmem=4.1G,s_vmem=4G";;

  3)
  echo " -q all.q -l h_vmem=6.1G,s_vmem=6G";;

  4)
  echo " -q all.q -l h_vmem=4.1G,s_vmem=4G";;

  5)
  echo " -q all.q -l h_vmem=6.1G,s_vmem=6G";;

  6)
  echo " -q all.q -l h_vmem=4.1G,s_vmem=4G";;

  7)
  echo " -q all.q -l h_vmem=4.1G,s_vmem=4G";;

esac
