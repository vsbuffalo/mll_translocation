"""
pipeline.py

This uses all the default directory names for the modules' functions
called.
"""

import pysam
import os
import sys
import re
import subprocess
from optparse import OptionParser

import find_odd_mates
import find_fusion
import divide_singles

if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()

    output_dir = "output"
    N = 6
    mll_template = "mll.fasta"
    if not os.path.exist(mll_template):
        parser.error("Could not MLL template file '%s'" % mll_template)

    if len(args) != 1:
        parser.error("Incorrect number of arguments.")
    elif not os.path.exists(args[0]):
        parser.error("File '%s' does not exist." % args[0])

    basename = os.path.basename(args[0])
    basename = basename.split('.')[0]

    ## Run find_odd_pairs on BARCODE.sam
    outdir = basename + "-odds-and-singles"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    ## side-effects: writes files.
    m = find_odd_mates.PairedReads(args[0])
    m.gather_pairs()
    m.find_odd_pairs()
    m.output_pairs(outdir=outdir)
    m.gather_unmapped_ends(outdir=outdir)
    

    ## Gather rearrangements files involving chr11
    odd_mates_dir = "%s/%s-chr*" % (outdir, basename)
    cmd = "wc -l %s | grep chr11 | sort -r" % odd_mates_dir
    
    proc = subprocess.Popen(cmd,
                            shell=True,
                            stdout=subprocess.PIPE)
    
    stdout_value = proc.communicate()[0]
    rearrangements = list()
    with open("%s/%s-rearrangement-counts.txt", 'w') as f:
        for line in stdout_value.split('\n'):
            line = line.strip()
            count, fn = re.split(r' +', line)
            chunks = re.split('[-\.]', fn)
            loc = '-'.join(chunks[4:6])
            rearrangement.append(loc)
            f.write("%s\t%s\n" % loc, count)

    top_candidates = rearrangements[:N]

    ## Call to R script to gather statistics
    # here

    ## Convert single file entries into FASTA files
    unmapped_dir = "%s-unmapped" % basename
    if not os.path.exists(unmapped_dir):
        os.mkdir(unmapped_dir)
    singles_filename = "%s/%s-singles.txt" % (outdir, basename)
    divide_singles.build_location_dict(singles_filename, unmapped_dir)
    

    ## Read singles file and run mapping
    mapping_dir = "%s-singles-mapping" % basename
    if not os.path.exists(mapping_dir):
        os.mkdir(mapping_dir)

    # get the FASTA files of those in top_candidates
    unmapped_fasta_files = os.listdir(mapping_dir)
    search_str = ';;;'.join(top_candidates)
    top_unmapped_fasta_files = [f for f in unmapped_fasta_files if f.split('-')[1] in search_str]
    
    # run BWA indexing
    subprocess.call("bwa index -a is %s" % mll_template, shell=True)
    
    # run BWA mapping on each of the top rearrangement candidates
    bwa_cmd = "bwa bwasw database.fasta %s > %s_aln.sam"
    for fasta_file in top_unmapped_fasta_files:
        ffbn = fasta_file.split('-')[1]
        subprocess.call(bwa_cmd % (fasta_file, ffbn))

    


    
