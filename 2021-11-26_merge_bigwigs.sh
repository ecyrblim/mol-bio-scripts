#! /bin/bash

printf 'This script takes replicate bigWig files and merges them.

Please ensure that: 
- All files to be merged are in one directory.
- All .bigWig files are present.

This script requires samtools, deeptools, and the UCSC tools (these 
can be found at http://hgdownload.soe.ucsc.edu/admin/exe/)
to be installed and in your path. 

Please ensure this is the case before running! \n\n'

printf 'Enter the directory containing your bigwig files: '
read BIGWIG_DIR

printf 'Enter a pattern common to the bigWig files to be merged, 
but unique relative to other bigWigs that should not be merged: '
read PATTERN

printf 'Enter the genome name, e.g. mm10, hg38: '
read SPECIES

cd $BIGWIG_DIR

printf 'The following files will be merged: \n\n'
ls -lah *$PATTERN*.bigWig

if [ $SPECIES == mm10 ]
then
    CHROM=https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/chromAlias.txt.gz
    SIZES=http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
elif [ $SPECIES == hg38 ]
then
    CHROM=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromAlias.txt.gz
    SIZES=http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
else
    printf 'Genome unsupported, please provide the UCSC URL to the
chromAlias.txt.gz file for your genome here: '
    read CHROM
    printf 'Provide the UCSC URL to the chrom.sizes file for your
genome here: '
    read SIZES
fi

printf '\nStep 1 of 2: Generating merged bedGraph file... \n'

bigWigMerge *$PATTERN*.bigWig $PATTERN.merged.bg

chromToUcsc -i $PATTERN.merged.bg -o $PATTERN.ucsc.bg -a $CHROM

bedSort $PATTERN.ucsc.bg $PATTERN.ucsc.bg 

printf 'Step 2 of 2: Coverting merged bedGraph file to bigWig... \n'

bedGraphToBigWig $PATTERN.ucsc.bg $SIZES $PATTERN.merged.bigWig

rm $PATTERN.merged.bg
rm $PATTERN.ucsc.bg

printf 'Done, enjoy your day :)'
