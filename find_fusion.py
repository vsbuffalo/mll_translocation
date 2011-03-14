"""
find-fusion.py

find-fusion.py takes a SAM file and finds all reads satisfying the
following conditions:

 - The CIGAR string is of the format xMyS where x and y are integers.
 - The read is mapped and mapped to the forward strand.

The input should be a SAM file of reads mapped to the MLL template
*with* BWA bwasw (which soft-clips well-matching reads). These reads
should the unmapped mates of reads mapped to the translocation
partners.
"""
from optparse import OptionParser
import pysam
import os
import sys

def extract_fusion_candidates(filename, outdir="", statsdir=""):
    """
    A fusion candidate is a read mapped to the forward strand with a
    CIGAR string in the format xMyS. This function gathers such reads
    and writes them to file, and returns statistics on the SAM file.
    """
    unmapped = 0
    length_not_2 = 0
    not_xMyS = 0
    reverse = 0
    total = 0
    correct = 0
    mapped = 0
    basename = os.path.basename(filename).split('.')[0]
    #fasta_file = open(os.path.join(outdir, "%s-hybrids.fasta" % basename), 'w')
    with open(os.path.join(outdir, "%s-fusion-candidates.txt" % basename), 'w') as f:
        for read in pysam.Samfile(filename):
            total += 1
            if read.is_unmapped:
                unmapped += 1
                continue
            mapped += 1
            # only handle CIGAR entries with length two (i.e. xMyS)
            if len(read.cigar) != 2:
                length_not_2 += 1
                continue

            e1, e2 = read.cigar
            strand = ('forward', 'reverse')[read.is_reverse]
            if not read.is_reverse:
                # note: 0 is M, 4 is S
                if e1[0] == 0 and e2[0] == 4:
                    matched = read.seq[:e1[1]]
                    cut = read.seq[e1[1]:]
                    break_pos = read.pos + e1[1]
                    correct += 1
                    f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (read.qname, break_pos,
                                                          matched, cut, strand, read.mapq))
                    #fasta_file.write(">%s;;%s\n%s\n" % (read.qname, break_pos, cut))
                else:
                    not_xMyS += 1
            else:
                reverse += 1
                if e1[0] == 0 and e2[0] == 4:
                    matched = read.seq[:e1[1]]
                    cut = read.seq[e1[1]:]
                    break_pos = read.pos + e1[1]
                    correct += 1
                    f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (read.qname, break_pos,
                                                          matched, cut, strand, read.mapq))
                    #fasta_file.write(">%s;;%s\n%s\n" % (read.qname, break_pos, cut))
                else:
                    not_xMyS += 1

                
    stats =  {'unmapped':unmapped,
              'mapped':mapped,
              'CIGAR entry length not 2':length_not_2,
              'CIGAR format not xMyS':not_xMyS,
              'mapped to reverse strand':reverse,
              'total':total,
              'correct':correct,
              'ignored':total-correct}
    
    # sort keys first so output format is always the same
    keys = ['unmapped', 'mapped', 'CIGAR entry length not 2',
            'CIGAR format not xMyS', 'mapped to reverse strand',
            'correct', 'ignored', 'total']

    with open(os.path.join(statsdir, "%s-fusion-stats.txt" % basename), 'w') as f:
        for k in keys:
            f.write("%s\t%s\n" % (k, stats[k]))

    return stats

            
                
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="dir",
                      help="directory to output files (default: {arg}-out)",
                      default=None, metavar="FILE")
    (options, args) = parser.parse_args()

    if options.dir is None:
        options.dir = os.path.dirname(args[0])
        
    if not os.path.exists(options.dir):
        os.mkdir(options.dir)

    if len(args) != 1:
        parser.error("Incorrect number of arguments.")
    elif not os.path.exists(args[0]):
        parser.error("File '%s' does not exist." % args[0])

    stats = extract_fusion_candidates(args[0], options.dir)
