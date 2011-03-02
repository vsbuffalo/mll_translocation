## assemble.py - take a FASTA file in which headers are mapped
## position, and assembly a consensus seq.

from optparse import OptionParser
import operator

class SequenceConsensus(object):
    """
    """

    def __init__(self, preallocate=100):
        """
        """
        self.length = preallocate
        self.consensus = [{'A':0, 'T':0, 'C':0, 'G':0, 'N':0} for i in range(preallocate)]
        self.seq = ""
        
    def add_seq(self, seq, pos):
        for i, base in enumerate(seq):
            try:
                self.consensus[pos+i][base] += 1
            except IndexError:
                self.consensus.append({'A':0, 'T':0, 'C':0, 'G':0, 'N':0})
                self.consensus[pos+i][base] += 1
                self.length += 1

    def pileup(self):
        for base_dict in self.consensus:
            tups = sorted(base_dict.items(), key=operator.itemgetter(1), reverse=True)
            self.seq += tups[0][0]
            print tups[0][0] + ' ' + '|'*tups[0][1] + ''.join([t[0]*t[1] for t in tups[1:]])


if __name__ == "__main__":
    parser = OptionParser()
    # parser.add_option("-l", "--len", dest="len",
    #                   help="max length of sequence (used for pre-allocation)",
    #                   default="mll_template/mll.fasta")

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
            seq = fasta_file.next().strip()
        seqs.append((seq, position))

    start = min([j for i, j in seqs])
    for seq, pos in seqs:
        consensus.add_seq(seq, pos-start)

    consensus.pileup()
