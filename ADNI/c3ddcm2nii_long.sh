#/bin/bash
#$ -S /bin/bash
set -e

if [ $# -lt 2 ]; then
  echo "Usage: $0 source_dir dest_dir [prefix]"
  echo "Traverses source_dir and extracts nifti images from all dicom series encountered"
  echo "Warning: Duplicate scans with the same series number and descriptions will be overwritten"
  exit 1
fi
if [ $# == 3 ]; then
  prefix=${3}_
fi

for dirname in $(find $1 -type d ); do
  #echo Looking in directory $dirname
  folder=$(echo $dirname | awk -F "/" '{print $NF}')
  outdir=${2}/$folder
  #mkdir -p $outdir
  sep=$(mktemp -u  | awk -F "/" '{print $NF}')
  ~srdas/bin/c3d -dicom-series-list "$dirname" | sed -n '2,$p' | sed -e "s/\t/$sep/g" | while read line; do 
    #echo "Extracting series $ser $dim $descr.."
    ser=$(echo $line | awk -F "$sep" '{print $1}'); 
    descr=$(echo $line | awk -F "$sep" '{print $4}'); 
    serid=$(echo $line | awk -F "$sep" '{print $NF}'); 
    dim=$(echo $line | awk -F "$sep" '{print $2}');
    fname=$( echo ser$(printf %02d $ser)_${descr}.nii.gz | sed -e 's/ /_/g' | sed -e "s/\//_/g" );
    # ~/bin/c3d -dicom-series-read "$dirname" "$serid" -o "${2}/$fname"
    mkdir -p $outdir;
    ~srdas/bin/c3d -dicom-series-read "$dirname" "$serid" -o "$outdir/$( echo ${prefix}ser$(printf %02d $ser)_${descr}.nii.gz | sed -e 's/ /_/g')";
    echo $outdir/$( echo ${prefix}ser$(printf %02d $ser)_${descr}.nii.gz | sed -e 's/ /_/g')
  done
done

