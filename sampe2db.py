"""
sam_db.py - Load SAM database results into a SQLite database.

"""
from optparse import OptionParser
import pysam
import sqlite3
import pdb
import os
from string import maketrans
complement_table = maketrans("ATCGatcg", "TAGCtagc")

def reverse_complement(seq):
    """
    Find the reverse complement of a sequence.
    """
    return seq.translate(complement_table)[::-1]


def create_table(conn, tablename, query):
    """
    Create a table if it doesn't exist.
    """
    check_query = """select name from sqlite_master where type='table' and name=?"""
    c = conn.cursor()
    c.execute(check_query, (tablename,))
    result = c.fetchone()
    if result is not None:
        result = result[0]
    if result != tablename:
        print "creating table '%s'" % tablename
        c.execute(query)
        conn.commit()

class PairedReads(object):
    """
    A class for gathering and loading paired-end reads into a
    database.
    """
    def __init__(self, filename, dbname, barcode=None):
        self.mode = "r"
        self.filename = filename
        if barcode is None:
            self.barcode = filename.split('.')[0]
        if filename.split('.')[1] == "bam":
            self.mode = "rb"
        self.samfile = pysam.Samfile(filename, self.mode)
        self.mapped = 0
        self.total = 0
        self.unmapped = 0
        self.reverse = 0
        self.unpaired = 0
        self.pairs = dict()
        self.dbname = dbname
        self.conn = sqlite3.connect(self.dbname)
        self.split_mates_table = 'split_mates'
        self.unmapped_mates_table = 'unmapped_mates'

        split_table = """create table %s
(barcode text,
name text,
chr_1 text, chr_2 text,
seq_1 text, seq_2 text,
pos_1 integer, pos_2 integer,
strand_1 text, strand_2 text,
mqual_1 integer, mqual_2 integer)""" % self.split_mates_table

        unmapped_table = """create table %s
(barcode text,
name text,
unmapped_seq text, unmapped_mqual integer, unmapped_strand text,
mapped_chr text, mapped_seq text,
mapped_pos integer, mapped_mqual integer, mapped_strand text)""" % self.unmapped_mates_table
        
        create_table(self.conn, self.split_mates_table, split_table)
        create_table(self.conn, self.unmapped_mates_table, unmapped_table)

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

        # Now, find odd pairs (mapped to different locations)
        self.find_odd_pairs()

        # Now, gather and push to table
        self.grouped_odd_pairs = dict()
        for qname in self.odd_pairs:
            locations = self.odd_pairs[qname].keys()
            locations.sort()
            # make chr11 first location for grouping
            if 'chr11' in locations:
                altchr = [k for k in locations if k != 'chr11'][0]
                locations = ['chr11', altchr]

            key = '-'.join(locations)

            result = self.grouped_odd_pairs.get(key, False)
            if not result:
                self.grouped_odd_pairs[key] = dict()
            self.grouped_odd_pairs[key][qname] = [self.odd_pairs[qname][k] for k in locations]

        # Write to database
        for key in self.grouped_odd_pairs:
            name_1, name_2 = key.split('-')
            for readname, readset in self.grouped_odd_pairs[key].items():
                strands = [("forward", "reverse")[int(r.is_reverse)] for r in readset]
                values = (self.barcode,
                          readname,
                          name_1, name_2,
                          readset[0].seq, readset[1].seq,
                          readset[0].pos, readset[1].pos,
                          strands[0], strands[1],
                          readset[0].mapq, readset[1].mapq)
                c = self.conn.cursor()
                c.execute("""insert into %s values
                (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""" % self.split_mates_table, values)
            self.conn.commit()
        
    def find_odd_pairs(self):
        """
        Given all pairs in a file, find the ones that mapped to
        different places on the reference.
        """
        self.odd_pairs = dict()
        for qname in self.pairs:
            locations = self.pairs[qname].keys()
            locations.sort()

            if len(locations) != 2:
                continue
            
            if locations[0] != locations[1]:
                reads = [self.pairs[qname][k] for k in locations]
                self.odd_pairs[qname] = dict(zip(locations, reads))        

    def gather_unmapped_ends(self):
        """
        Gather and output reads that have one side that mapped and
        another end that isn't mapped. These are candidates for the
        unmapped read containing the breakpoint and part of the
        translocated chromosome.

        Note: unmapped entries' sequences can be in the SAM file as
        reverse complements, so we need to check the bitflag to see if
        they are, and possibly reverse complement the sequence *back*
        to what it is in the original FASTA/FASTAQ file.
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

        # Now push into table, reverse complementing unmapped reads
        for qname, readset in self.pairs_one_unmapped.items():
            unmapped = [i for i, r in enumerate(readset) if r.is_unmapped][0]
            mapped = [i for i, r in enumerate(readset) if not r.is_unmapped][0]
            reads = [readset[unmapped], readset[mapped]]
            strands = [("forward", "reverse")[int(r.is_reverse)] for r in reads]
            
            # reverse complement reads *back* if they are unmapped and reverse
            if (reads[0].is_reverse):
                unmapped_seq = reverse_complement(reads[0].seq)
            else:
                unmapped_seq = reads[0].seq

            values = (self.barcode, qname,
                      unmapped_seq, reads[0].mapq, strands[0],
                      self.samfile.getrname(reads[1].rname),
                      reads[1].seq, reads[1].pos, reads[1].mapq, strands[1])
            c = self.conn.cursor()
            c.execute("""insert into %s values
            (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""" % self.unmapped_mates_table, values)
        self.conn.commit()



if __name__ == "__main__":
    parser = OptionParser()
    options, args = parser.parse_args()

    if len(args) < 1:
        parser.error("SAM file argument required.")

    if not os.path.exists(args[0]):
        parser.error("File in argument does not exist.")

    dbname = os.path.basename(args[0]).split('.')[0]
    m = PairedReads(args[0], dbname)
    m.gather_pairs()
    m.gather_unmapped_ends()
