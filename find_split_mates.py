"""
find-split-mates.py

find-split-mates.py takes a SAM file from a paired-end mapping (only
tested with BWA), and returns:

 - files (in a directory) of all paired-end reads with mates that
   mapped to different locations/chromosomes.
 - a file of all reads that have one mate mapped and another mate unmapped.
 - statistics (to stdout) of the SAM file.
 
"""
from optparse import OptionParser
import pysam
import os
import sys


class PairedReads(object):
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
        self.reverse = 0
        self.unpaired = 0
        self.pairs = dict()

    @property
    def stats(self):
        """
        Return a dictionary of statistics.
        """
        return {'unmapped':self.unmapped,
                'mapped':self.mapped,
                'mapped to reverse strand':self.reverse,
                'total':self.total,
                'unpaired':self.unpaired}
        
    def add_pair(self, read):
        """
        Add a read to a dictionary by searching if the sequence header
        name has been inserted in the dictionary already. This gathers
        the pairs together.
        
        Note that this uses "rname", which is depricated. Pysam-0.4
        requires "tid" to be used instead.
        """
        refname = self.samfile.getrname(read.rname)
        result = self.pairs.get(read.qname, False)

        if read.is_reverse:
            self.reverse += 1

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
                    self.unpaired += 1

                self.mapped += 1
                self.add_pair(read)              
            else:
                self.unmapped += 1
        

    def gather_unmapped_ends(self, outdir):
        """
        Gather and output reads that have one side that mapped and
        another end that isn't mapped. These are candidates for the
        unmapped read containing the breakpoint and part of the
        translocated chromosome.
        """
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
                format_line = "%s\t" * 8 + "\n"
                f.write(format_line % (qname,
                                       reads[0].seq, reads[0].pos,
                                       reads[0].mapq,
                                       self.samfile.getrname(reads[1].rname),
                                       reads[1].seq, reads[1].pos,
                                       reads[1].mapq))
        
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

    def output_pairs(self, outdir):
        """
        Write a seperate file for each paired-end combination.

        First, join all locations into a single key, the for each of
        these keys, output a file with name {basename}_pairs.txt.
        """
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
            filename = self.basename + '-' + key + ".txt"
            if dir is not None:
                filename = os.path.join(outdir, filename)
                if not os.path.exists(outdir):
                    os.mkdir(outdir)
            with open(filename, 'w') as f:
                for readname, readset in self.grouped_odd_pairs[key].items():
                    strands = [("forward", "reverse")[int(r.is_reverse)] for r in readset]
                    format_line = "%s\t" * 9 + "\n"
                    f.write(format_line % (readname,
                                           readset[0].seq, readset[1].seq,
                                           readset[0].pos, readset[1].pos,
                                           strands[0], strands[1],
                                           readset[0].mapq, readset[1].mapq))

    
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="dir",
                      help="directory to output files (default: name of first arg)",
                      default=None)
    (options, args) = parser.parse_args()

    # Get basename (useful for barcoded data), make directory if doesn't exist
    basename = os.path.basename(args[0])
    basename = basename.split('.')[0]
    if options.dir is None:
        options.dir = basename + "-odds-and-singles"
    if not os.path.exists(options.dir):
        os.mkdir(options.dir)

    if len(args) != 1:
        parser.error("Incorrect number of arguments.")
    elif not os.path.exists(args[0]):
        parser.error("File '%s' does not exist." % args[0])

    # Gather and output all relevant files
    m = PairedReads(args[0])
    m.gather_pairs()
    m.find_odd_pairs()
    m.output_pairs(outdir=options.dir)
    m.gather_unmapped_ends(outdir=options.dir)
    
    
    # sort keys first so output format is always the same
    stats = m.stats
    keys = stats.keys()
    keys.sort()

    for k in keys:
        sys.stdout.write("%s\t%s\n" % (k, stats[k]))
