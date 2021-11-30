#! /bin/bash

printf "This script finds peaks within n bases of the TSS. Its input is a Homer-annotated region file 
(specifically, a BED file processed using Homer annotatePeaks) and your desired value for n.

It generates two files per value of n: a filtered version of the annotated region file containing
only the peaks within n bases of a TSS, and a file containing genes with peaks and the "

read -p "Enter the path to your annotated peak file: " -r PEAK_DIR

read -p "Enter your desired distance(s) from TSS. If you input n, then this script will find peaks in the region -n to +n from the TSS. You may input multiple values separated by commas: " -r DIST_INPUT

DISTANCES=$(echo $DIST_INPUT | tr "," "\n")

FILE=${PEAK_DIR##*/}
printf "Filename is: $FILE \n"
BASE=${FILE%.*}
printf "Basename is $BASE \n"
DIR=${PEAK_DIR%/*}
printf "Input and output directory is $DIR \n"

cd $DIR

printf 'These distances from the TSS will be used: \n'
echo $DISTANCES

for DIST in $DISTANCES
do
  printf "Processing peaks ± $DIST base pairs from the TSS...  \n"
  OUTFILE=${BASE}.${DIST}bp.tsv
  NEG=-${DIST}
  awk -F "\t" -v dist=$DIST -v neg=$NEG 'BEGIN{FS="\t"; OFS=FS} NR==1; NR > 1{ if ($10 <= dist && $10 >= neg) { print } }' $PEAK_DIR > ${DIR}/${OUTFILE}
  awk -F "\t" 'BEGIN{FS="\t"; OFS=FS} { a[$16]++ } END { for (n in a) print n, a[n] }' ${DIR}/${OUTFILE} > ${DIR}/${BASE}.${DIST}bp.genecounts.tsv
done


