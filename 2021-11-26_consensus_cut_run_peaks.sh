#! /bin/bash

printf 'This script takes SEACR peak files in the IGV directory
of the nf-core CUT&RUN pipeline output and processes them
to generate consensus peak sets, filter these by motif
presence, then annotate them with Homer. 

It requires bedtools, bedmap, and Homer to be installed 
and in your path. Please ensure this is the case before
running! \n\n'

printf 'Enter the IGV directory containing your peak files: '

read IGV_DIR

printf 'Enter the path to a bed file with mapped motifs: '

read MOTIF_DIR

printf 'Enter the genome name (e.g. mm10, hg38): '

read SPECIES

# export IGV_DIR=/Volumes/A340-Projekte/Bryce/betti/cut_run_2/cut_run_mouse_02-1_wt
# export MOTIF_DIR=/Volumes/A340-Projekte/Bryce/betti/cut_run_2/motif/AAAGTTW_mm10.bed
# export SPECIES=mm10

printf 'Thank you. Processing: \n\n\n\n'

cd $IGV_DIR 

mkdir final_peaks
mkdir final_peaks/1.cat
mkdir final_peaks/2.bedmap_0.25
mkdir final_peaks/3.flat
mkdir final_peaks/4.motif_intersect
mkdir final_peaks/5.stats

printf 'Step 1 of 4: Concatenating peak files... \n\n\n\n'

cat ./*.bed.stringent.bed > ./final_peaks/1.cat/cat_unsorted.bed
bedtools sort -i ./final_peaks/1.cat/cat_unsorted.bed > ./final_peaks/1.cat/cat.bed

rm ./final_peaks/1.cat/cat_unsorted.bed

printf 'Step 2 of 4: Merging replicates with 50 percent overlap in two or more peaks... \n\n\n\n'

bedmap --count --echo-map-range --fraction-both 0.5 --delim '\t' ./final_peaks/1.cat/cat.bed \
    | awk '$1>1' - \
    | cut -f2- - \
    | sort-bed - \
    | uniq - \
    > ./final_peaks/2.bedmap_0.25/cat.merge.bed

bedtools merge -i  ./final_peaks/2.bedmap_0.25/cat.merge.bed > ./final_peaks/3.flat/cat.flat.bed

printf 'Step 3 of 4: Filtering peaks by motif presence... \n\n\n\n'

bedtools intersect -a  ./final_peaks/3.flat/cat.flat.bed  -b $MOTIF_DIR  -wa > ./final_peaks/4.motif_intersect/cat.flat.motif_wa.unsorted.bed

sort -k1,1 -k2,2n -k3,3n -u ./final_peaks/4.motif_intersect/cat.flat.motif_wa.unsorted.bed > ./final_peaks/4.motif_intersect/cat.flat.motif_wa.bed

rm ./final_peaks/4.motif_intersect/cat.flat.motif_wa.unsorted.bed

printf 'Step 4 of 4: Annotating peaks... \n\n\n\n'

annotatePeaks.pl ./final_peaks/4.motif_intersect/cat.flat.motif_wa.bed  $SPECIES -annStats ./final_peaks/5.stats/cat.flat.motif_wa.anno_stats.tsv > ./final_peaks/5.stats/cat.flat.motif_wa.anno.tsv

printf 'Done, enjoy your day :) \n\n\n\n'
