#/bin/bash
#$ -S /bin/bash

ASHS_WORK=$1
STAGE=$2

Q="-q all.q"
QHM="-q all.q,himem.q"

case $STAGE in

  1)
  echo " $Q  -l h_vmem=10.1G,s_vmem=10G -pe serial 1 ";;

  2)
  echo " $Q  -l h_vmem=10.1G,s_vmem=10G";;

  3)
  echo " $Q  -l h_vmem=10.1G,s_vmem=10G";;

  4)
  echo " $Q  -l h_vmem=4.1G,s_vmem=4G";;

  5)
  echo " $Q  -l h_vmem=6.1G,s_vmem=6G";;

  6)
  echo " $Q  -l h_vmem=4.1G,s_vmem=4G";;

  7)
  echo " $Q  -l h_vmem=4.1G,s_vmem=4G";;

esac
