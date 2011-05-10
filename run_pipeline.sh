#!/bin/bash
# Usage:
# bash run_pipeline.sh data/dir_with_SAM_files num_threads

if [ -z $1 ]
then
    echo "Usage: bash run_pipeline.sh data/dir_with_SAM_files [num_threads]"
    exit 1
fi

if [ -z $2 ]
then
    num_threads="1"
else
    num_threads=$2
fi

rootdir=$(basename $1)
printf "Processing all SAM files in '%s'\n" $rootdir
printf "Using %s thread(s)\n\n" $num_threads

if [ ! -f results ]
then
    echo "Making directory 'results/'"
    mkdir -p results
fi

# Check for all bam, sorted.bam, and sorted.bam.bai files
barcodes=$(find data/$rootdir/*.sam | xargs -n1 basename | sed s/.sam//)
coversion_needed=0

for bc in $barcodes
do
    found_ext=$(find data/$rootdir/$bc* | xargs -n1 basename | sort | sed 's/[^\\.]*.//' | tr '\n' ';')
    if [ ! $found_ext == "bam;sam;sorted.bam;sorted.bam.bai;" ]
    then
        coversion_needed=1
    fi
done

if [ $coversion_needed -eq 1 ]
then
    # Take all SAM files and convert them to BAM
    echo "Converting all SAM files to BAM files."
    find data/$rootdir/*.sam | sed s/.sam// | \
        xargs -n1 -P $num_threads -I{} samtools view -h -S -b -o {}.bam {}.sam

   # Take all BAM files and sort
    echo "Sorting all BAM files."
    find data/$rootdir/*.bam | sed s/.bam// | \
        xargs -n1 -P $num_threads -I{} samtools sort {}.bam {}.sorted
    
   # Index all BAM files
    echo "Indexing all BAM files."
    find data/$rootdir/*.sorted.bam | sed s/.sorted.bam// | \
        xargs -n1 -P $num_threads -I{} samtools index {}.sorted.bam
else
    printf "All barcoded files had the correct extensions indicating \nthey were converted to sorted and index BAM files.\n\n"
fi


if [ ! -f results/$rootdir/ ]
then
    printf "Making directory 'results/%s'" $rootdir
    mkdir -p results/$rootdir
fi

# Run analysis pipeline on all BAM files
echo "Beginning data mining analysis of all sorted and indexed BAM files."
find data/$rootdir/*.sorted.bam | xargs -n1 -P $num_threads -I{} Rscript R/tlminer.R {} $rootdir 2> results/$rootdir/log.txt

# Compress results for post-analysis
#find . -name "*-output" | perl -ne 'chomp; print $_ . "/results/ "' | xargs tar -czf results.tar.gz
