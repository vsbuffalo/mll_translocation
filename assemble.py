## assemble.py - take a FASTA file in which headers are mapped
## position, and assembly a consensus seq.

import sys
from optparse import OptionParser
import operator

class SequenceConsensus(object):
    """
    """

    def __init__(self, preallocate=10000):
        """
        """
        self.length = preallocate
        self.consensus = [{'A':0, 'T':0, 'C':0, 'G':0, 'N':0} for i in range(preallocate)]
        self.seq = ""
        self.positions = list()
        
    def add_seq(self, seq, pos):
        for i, base in enumerate(seq):
            try:
                self.consensus[pos+i][base] += 1
            except IndexError:
                self.consensus.append({'A':0, 'T':0, 'C':0, 'G':0, 'N':0})
                self.consensus[pos+i][base] += 1
                self.length += 1

    def pileup(self, width=100, output=False):
        max_count = max([max(d.values()) for d in self.consensus])
        collapse = int(max_count/width)
        for base_dict in self.consensus:
            tups = sorted(base_dict.items(), key=operator.itemgetter(1), reverse=True)
            if max([c for b, c in tups]) > 0:
                if output:
                    print tups[0][0] + ' ' + '|'*(tups[0][1]/collapse) + ''.join([t[0]*t[1] for t in tups[1:]])
                self.seq += tups[0][0]



if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-p", "--pileup", dest="pileup",
                      help="print pileup to screen",
                      default=False, action="store_true")
    parser.add_option("-w", "--width", dest="width",
                      help="screen width (for histogram)",
                      default=100)

    options, args = parser.parse_args()

    if len(args) != 1:
        parser.error("Incorrect number of arguments supplied.")

    consensus = SequenceConsensus()
    fasta_file = open(args[0])

    seqs = list()
    
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            position = int(line[1:])
            consensus.positions.append(position)
            seq = fasta_file.next().strip()
        seqs.append((seq, position))

    start = min([j for i, j in seqs])
    for seq, pos in seqs:
        consensus.add_seq(seq, pos-start)

    if options.pileup:
        consensus.pileup(int(options.width), output=True)
    else:
        consensus.pileup()
        sys.stdout.write(">pos.range:%s,%s\n" % (min(consensus.positions), max(consensus.positions)))
        sys.stdout.write(consensus.seq + '\n')
