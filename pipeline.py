"""
MLL translocation pipeline

Input: SAM file from paired-end reads mapped to the entire human
genome. High coverage reads must be quality-controlled by adaptive
quality trimming and sequence adapter contaminant trimming.

"""

import pysam
import os
import sys
import re
import subprocess
from optparse import OptionParser
import logging

import find_odd_mates

LEVELS = {'debug': logging.DEBUG,
          'info': logging.INFO,
          'warning': logging.WARNING,
          'error': logging.ERROR,
          'critical': logging.CRITICAL}


def check_dir(direc, parent=None):
    """
    Makes a directory and logs it. 
    """
    if parent is not None:
        direc = os.path.join(parent, direc)
    if not os.path.exists(direc):
        logging.info("Making directory '%s'." % direc)
        os.mkdir(direc)
    else:
        logging.info("Directory '%s' exists, rewriting old files!" % direc)
    return direc

def get_rearrangements(splitdir, statsdir, filename):
    """
    Gather rearrangements files involving chr11 from wc -l output.
    """
    cmd = "wc -l %s/* | grep chr11 | sort -rn" % splitdir
    
    proc = subprocess.Popen(cmd,
                            shell=True,
                            stdout=subprocess.PIPE)
    
    stdout_value = proc.communicate()[0]
    rearrangements = list()
    ra_file = "%s/%s" % (statsdir, filename)
    with open(ra_file, 'w') as f:
        for line in stdout_value.split('\n'):
            line = line.strip()
            if len(line) == 0:
                continue
            count, fn = re.split(r' +', line)
            fn = os.path.basename(fn)
            chunks = re.split('[-\.]', fn)
            loc = '-'.join([c for c in chunks if 'chr' in c])
            rearrangements.append(loc)
            f.write("%s\t%s\n" % (loc, count))
    return ra_file
    
def find_split_mates(samfile, outdir, no_unmapped=False):
    """
    Wrapper for finding split-mates sub-pipeline.
    """
    logging.info("Finding split-mates...")
    m = find_odd_mates.PairedReads(samfile)
    m.gather_pairs()
    m.find_odd_pairs()
    logging.info("Writing split-mates files...")
    m.output_pairs(outdir=outdir)
    if not no_unmapped:
        logging.info("Finding paired-end reads with one unmapped mate...")
        m.gather_unmapped_ends(outdir=outdir)

    


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-r", "--ref", dest="ref",
                      help="reference template (default: mll.fasta)",
                      default="mll.fasta")

    (options, args) = parser.parse_args()

    ### Setup ###
    # Number of rearrangements to consider
    N = 4
    
    if len(args) != 1:
        parser.error("Incorrect number of arguments supplied.")

    ## Set up logging - this can be an option
    level_name = 'info'
    level = LEVELS.get(level_name, logging.NOTSET)
    logging.basicConfig(level=level)

    ## Grab basename from filename - useful for barcoded data
    basename = re.split(r'[-\.]', os.path.basename(args[0]))[0]
    logging.info("Using '%s' as basename." % basename)    
    
    ## Check input SAM file
    if not os.path.exists(args[0]):
        parser.error("File '%s' does not exist." % args[0])

    ## Create output directory
    output_dir = check_dir(basename + "-output")

    ## Check reference template
    if not os.path.exists(options.ref):
        parser.error("Reference template FASTA file '%s' does not exist." % options.ref)

    ### Analysis ###

    ## Make stats directory
    stats_dir = check_dir("stats", output_dir)

    ## Part 1: Split the SAM file into reads mapped to different
    ## locations, and reads in which one mate mapped, but another didn't.
    split_mates_dir = check_dir("split-mates", output_dir)

    if False:
        find_split_mates(args[0], split_mates_dir)

    ## Part 2: Run Samtools to find reliable mappings, to get clearer
    ## picture of reads mapped to different chromosomes (_not_ the
    ## singles file).

    # before mapping quality pruning
    ra_file = get_rearrangements(split_mates_dir, stats_dir, "rearrangement-counts.txt")

    # run samtools quality pruning on SAM file
    logging.info("Using Samtools to remove 0-mapping-quality reads in '%s'." % args[0])
    reliable_sam_file = basename + "-reliable.sam"
    cmd = "samtools view -S -h -q1 %s > %s" % (args[0], reliable_sam_file)
    logging.info("Running command: '%s'." % cmd)
    subprocess.call(cmd, shell=True)

    # Re-run find_split_mates sub-pipeline on mapping-quality-trimmed
    # data Note: We can't run the -singles (one mate mapped, the other
    # unmapped) file with mapping quality subsetting.
    reliable_split_dir = check_dir("split-mates-reliable", output_dir)
    find_split_mates(reliable_sam_file, reliable_split_dir, no_unmapped=True)
    
    # Re-running finding split-mate pipeline on quality-pruned data
    qra_file = get_rearrangements(reliable_split_dir, stats_dir,
                                  "reliable-rearrangement-counts.txt")

    ## Part 3: In "reliable" split-mapped (with chr11) reads, extract positions.
    rearrangements = dict()
    with open(qra_file) as f:
        for line in f:
            line = line.strip()
            ra, count = line.split('\t')
            rearrangements[ra] = count

    
    
