MLL Translocation Miner
=======================

* Running tlminer.R

`tlminer.R` can be run in two ways: via command line, and
interactively. Either way, it requires a BAM file that has been sorted
and indexed.

** Running via command line

The helper script `run_pipeline.sh` is probably the best way to run
`tlminer.R`. It takes two arguments:

    bash run_pipeline.sh dat/dir_with_SAM_files [num_threads]

However, you can also run `tlminer.R` with `Rscript`:

    Rscript tlminer.R sorted_bam_file.bam groupdir

where `groupdir` is a directory inside of results (i.e. "original_run"
or "barcode_run_0x").

** Running interactively

This is the best way to run `tlminer.R` when you want to manually
check some results or debug the pipeline.


First, open up R. Then set:

    mapfile = "data/samples/your_sorted_bam.bam"

and then:

    source("R/tlminer.R")


`tlminer.R` looks for whether this value exists in the global
environment before checking for it in the arguments. Setting it will
run this process on that file.

