"""
"""
import sys
import os
import pysam
from optparse import OptionParser
import pdb

class MappedPairs(object):
    """
    A class for gathering, manipulating, and outputing paired-end
    reads, specifically those on different chromosomes.
    """

    def __init__(self, filename):
        self.mode = "r"
        self.filename = filename
        self.basename = filename.split('.')[0]
        
        if filename.split('.')[1] == "bam":
            self.mode = "rb"
        self.samfile = pysam.Samfile(filename, self.mode)
        self.mapped = 0
        self.total = 0
        self.unmapped = 0
        self.reverse_mapped = 0
        self.pairs = dict()
        
    def paired_add(self, read):
        """
        Add a read to a dictionary by searching if the sequence header
        name has been inserted in the dictionary already. This gathers
        the pairs together.
        
        Note that this uses "rname", which is depricated. Pysam-0.4
        requires "tid" to be used instead.
        """
        refname = self.samfile.getrname(read.rname)
        result = self.pairs.get(read.qname, False)

        if not result:
            self.pairs[read.qname] = dict()
        self.pairs[read.qname][refname] = read

                    
    def gather_pairs(self):
        """
        Iterate through all SAM file reads and add files to a
        dictionary with their sequence header name as key.
        """
        for read in self.samfile:
            self.total += 1
            if not read.is_unmapped:
                if not read.is_paired:
                    sys.exit()

                self.mapped += 1
                self.paired_add(read)              
            else:
                self.unmapped += 1
        

    def gather_unmapped_ends(self, outdir=None):
        """
        Gather and output reads that have one side that mapped and
        another end that isn't mapped. These are candidates for the
        unmapped read containing the breakpoint and part of the
        translocated chromosome.
        """
        if outdir is None:
            outdir = self.basename + "-out"
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        # First, build up a hash of all reads with the key being
        # qname. If a read exists and is mapped, (and the current read
        # during the iteration is likewise so), pop the existing and
        # don't add the current read.
        self.pairs_one_unmapped = dict()
        self.samfile = pysam.Samfile(self.filename, self.mode)
        for read in self.samfile:
            in_sample = self.pairs_one_unmapped.get(read.qname, False)
            if in_sample:
                assert(len(in_sample) == 1)
                # If both reads are mapped, pop this off.
                if in_sample[0].is_unmapped != read.is_unmapped:
                    self.pairs_one_unmapped[read.qname].append(read)
                else:
                    self.pairs_one_unmapped.pop(read.qname)
            else:
                self.pairs_one_unmapped[read.qname] = [read]

        # Now, output to file
        with open(os.path.join(outdir, self.basename + '-singles.txt'), 'w') as f:
            for qname, readset in self.pairs_one_unmapped.items():
                #pdb.set_trace()                                
                unmapped = [i for i, r in enumerate(readset) if r.is_unmapped][0]
                mapped = [i for i, r in enumerate(readset) if not r.is_unmapped][0]
                reads = [readset[unmapped], readset[mapped]]
                f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (qname,
                                                      reads[0].seq, reads[0].pos, 
                                                      self.samfile.getrname(reads[1].rname),
                                                      reads[1].seq, reads[1].pos))
        
    def find_odd_pairs(self):
        """
        Given all pairs in a file, find the ones that mapped to
        different places on the reference.
        """
        self.odd_pairs = dict()
        for qname in self.pairs:
            locations = self.pairs[qname].keys()
            if len(locations) != 2:
                continue
            if locations[0] != locations[1]:
                reads = [self.pairs[qname][k] for k in locations]
                self.odd_pairs[qname] = dict(zip(locations, reads))

    def output_pairs(self, basename=None, outdir=None):
        """
        Write a seperate file for each paired-end combination.

        First, join all locations into a single key, the for each of
        these keys, output a file with name {basename}_pairs.txt.
        """
        if basename is None:
            basename = self.basename
        if outdir is None:
            outdir = basename + "-out"
        
        self.grouped_odd_pairs = dict()
        for qname in self.odd_pairs:
            locations = self.odd_pairs[qname].keys()
            key = '-'.join(locations)

            result = self.grouped_odd_pairs.get(key, False)
            if not result:
                self.grouped_odd_pairs[key] = dict()
            self.grouped_odd_pairs[key][qname] = [self.odd_pairs[qname][k] for k in locations]

        # Write to file
        for key in self.grouped_odd_pairs:
            filename = basename + '-' + key + ".txt"
            if dir is not None:
                filename = os.path.join(outdir, filename)
                if not os.path.exists(outdir):
                    os.mkdir(outdir)
            with open(filename, 'w') as f:
                for readname, readset in self.grouped_odd_pairs[key].items():
                    strands = [("forward", "reverse")[int(r.is_reverse)] for r in readset]                                    
                    f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (readname,
                                                          readset[0].seq, readset[1].seq,
                                                          readset[0].pos, readset[1].pos,
                                                          strands[0], strands[1]))
                    

    def extract_forward_split_matches(self, outdir=None):
        """
        Given a SAM file, gather all forward matches that have a xMyS
        format, indicating that the match has been cut by BW
        bwasw. The position of these and the forward and ending
        sequences are sent to a file.
        """
        self.tossed = 0
        if outdir is None:
            outdir = self.basename + "-out"
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        with open(self.basename + "-split-unmapped.txt", 'w') as f:
            for read in self.samfile:
                if read.is_unmapped:
                    self.unmapped += 1
                    continue

                # only handle CIGAR entries with length two (i.e. xM yS)
                if len(read.cigar) != 2:
                    continue
                e1, e2 = read.cigar
                if not read.is_reverse:
                    if e1[0] == 0 and e2[0] == 4:
                        matched = read.seq[:e1[1]]
                        cut = read.seq[e1[1]:]
                        break_pos = read.pos + e1[1]
                        f.write("%s\t%s\t%s\t%s\tforward\n" % (read.qname, break_pos,
                                                               matched, cut))
                    else:
                        self.tossed += 1
                else:
                    self.reverse_mapped += 1
                    self.tossed += 1

        print "unmapped: %s\ntossed: %s\nreverse: %s\n" % (self.unmapped,
                                                           self.tossed,
                                                           self.reverse_mapped)


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="dir",
                      help="directory to output files (default: name of first arg)",
                      default=None)
    (options, args) = parser.parse_args()

    if len(args) < 2:
        parser.error("Not enough arguments.")
        if not os.path.exist(args[1]):
            parser.error("File '%s'does not exist." % args[1])
    elif len(args) > 2:
        parser.error("Too many arguments supplied; specify a single SAM or BAM file.")    

    if args[0] == 'findpairs':
        # Gather all paired ends, find the odd pairs (on different
        # chromosomes), and output pairs grouped by chromosome combination
        # to files.
        m = MappedPairs(args[1])
        m.gather_pairs()
        m.find_odd_pairs()
        m.output_pairs(outdir=options.dir)
        
        m.gather_unmapped_ends()

        # Write statistics about file
        basename = args[1].split('.')[0]
        with open(basename + "_stats.txt", 'w') as f:
            f.write("total: %d\nunmapped: %d\mapped: %d\nodd pairs: %d\n" 
                    % (m.total, m.unmapped, m.mapped, len(m.odd_pairs)))

    elif args[0] == 'findsplit':
       m = MappedPairs(args[1])
       m.extract_forward_split_matches(outdir=options.dir)
    else:
        parser.error("Unknown subcommand '%s'." % args[0])
