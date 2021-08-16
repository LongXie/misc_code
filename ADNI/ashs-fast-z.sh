#/bin/bash
#$ -S /bin/bash

ASHS_WORK=$1
STAGE=$2

Q=""

case $STAGE in

  1)
  echo " $Q -M 6G";;

  2)
  echo " $Q -M 4G";;

  3)
  echo " $Q -M 6G";;

  4)
  echo " $Q -M 4G";;

  5)
  echo " $Q -M 6G";;

  6)
  echo " $Q -M 10G";;

  7)
  echo " $Q -M 10G";;

esac
