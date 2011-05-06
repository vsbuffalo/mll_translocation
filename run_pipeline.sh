#!/bin/bash

# Take all SAM files and convert them to BAM
echo *.sam | sed s/.sam/\\n/g | xargs -n1 -P10 -I{} samtools view -S -b -o {}.bam {}.sam

# Take all BAM files and sort
echo *.bam | sed s/.bam/\\n/g | xargs -n1 -P10 -I{} samtools sort {}.bam {}.sorted

# Index all BAM files
echo *.sorted.bam | sed s/.sorted.bam/\\n/g | xargs -n1 -P10 -I{} samtools index {}.sorted.bam

# Run analysis pipeline on all *.db SQLite databases
echo *.sorted.bam | sed s/.sorted.bam/\\n/g | xargs -n1 -P10 Rscript tlminer.R

# Compress results for post-analysis
find . -name "*-output" | perl -ne 'chomp; print $_ . "/results/ "' | xargs tar -czf results.tar.gz
