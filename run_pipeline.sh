#!/bin/bash

# Run analysis pipeline on all *.db SQLite databases
ls *.db | xargs -n1 -P10 Rscript analysis.R 

# Compress results for post-analysis
find . -name "*-output" | perl -ne 'chomp; print $_ . "/results/ "' | xargs tar -czf results.tar.gz
