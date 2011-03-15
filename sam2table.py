"""
Extract fields from a SAM file of alignments from BWA samse.
"""
from optparse import OptionParser
import pysam
import os
import sys
import csv

def make_cigar(cigarlist):
    """
    Turn pysam's CIGAR list back to into a string.
    """
    if cigarlist is None:
        return ""
    cigar = "MIDNSHP=X"
    return ''.join([str(len) + cigar[op] for op, len in cigarlist])

def mapped_reads_table(samfile, tablefile):
    """
    
    """
    tf = open(tablefile, 'w')
    table_writer = csv.writer(tf, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
    with pysam.Samfile(samfile, 'r') as sf:
        first = True
        for read in sf:
            if first:
                table_writer.writerow(['name', 'cigar', 'strand', 'mapq', 'pos', 'rname', 'seq'])
                first = False
            strand = ('forward', 'reverse')[read.is_reverse]
            if not read.is_unmapped:
                table_writer.writerow([read.qname, make_cigar(read.cigar), strand, read.mapq, read.pos, sf.getrname(read.rname), read.seq])
    


if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()

    basename = os.path.basename(args[0]).split('.')[0]
    mapped_reads_table(args[0], os.path.join(os.path.dirname(args[0]), basename+'.txt'))
