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
        self.pairs_table = 'pairs'

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

        pairs_table = """create table %s
(barcode text,
name text,
chr_1 text, chr_2 text,
seq_1 text, seq_2 text,
pos_1 integer, pos_2 integer,
strand_1 text, strand_2 text,
mqual_1 integer, mqual_2 integer,
proper_1 int, proper_2 int)
""" % self.pairs_table

        create_table(self.conn, self.split_mates_table, split_table)
        create_table(self.conn, self.unmapped_mates_table, unmapped_table)
        create_table(self.conn, self.pairs_table, pairs_table)

    def add_pair(self, read):
        """
        Add a read to a dictionary by searching if the sequence header
        name has been inserted in the dictionary already. This groups
        the pairs together.
        """
        result = self.pairs.get(read.qname, False)

        if read.is_reverse:
            self.reverse += 1

        if not result:
            self.pairs[read.qname] = list()
        self.pairs[read.qname].append(read)


    def process_pairs(self, first='chr11'):
        """
        For each pair, add to approperiate table.

        There is some redundant code here, so this is a candidate to
        be refactored.
        """
        c = self.conn.cursor()
        for pair in self.pairs:
            readset = self.pairs[pair]

            # if len(readset) < 2:
            #     self.unpaired += 1
            #     continue
            assert(len(readset) == 2)

            all_mapped = all([not r.is_unmapped for r in readset])
            mapping_positions = [r.rname for r in readset]
            
            if len(set(mapping_positions)) == 2 and all_mapped:
                # reads are split mates; all mapped and on different chromosomes

                # Check that there are no unmapped reads this far,
                # according to reference mapping location (if it's -1)
                assert(-1 not in [r.rname for r in readset])

                # Make dictionary of reference mapping locations and reads
                keys = [self.samfile.getrname(r.rname) for r in readset]
                values = [r for r in readset]
                readset = dict(zip(keys, values))

                if first in readset.keys():
                    # make this the first element in reads
                    reads = [readset[first], [v for k, v in readset.items() if k != first][0]]
                    rnames = [first, [k for k in readset.keys() if k != first][0]]
                else:
                    reads = readset.values()
                    rnames = readset.keys()

                strands = [("forward", "reverse")[int(r.is_reverse)] for r in reads]
                values = (self.barcode,
                          reads[0].qname,
                          rnames[0], rnames[1],
                          reads[0].seq, reads[1].seq,
                          reads[0].pos, reads[1].pos,
                          strands[0], strands[1],
                          reads[0].mapq, reads[1].mapq)
                c.execute("""insert into %s values
                (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""" % self.split_mates_table, values)

            elif len(set(mapping_positions)) == 1 and all_mapped:
                # reads are pairs mapped to same chromosome.

                # Again, check that there are no unmapped reads this far:
                assert(-1 not in [r.rname for r in readset])

                strands = [("forward", "reverse")[int(r.is_reverse)] for r in readset]
                values = (self.barcode,
                          readset[0].qname,
                          self.samfile.getrname(readset[0].rname),
                          self.samfile.getrname(readset[1].rname),
                          readset[0].seq, readset[1].seq,
                          readset[0].pos, readset[1].pos,
                          strands[0], strands[1],
                          readset[0].mapq, readset[1].mapq,
                          int(readset[0].is_proper_pair), int(readset[1].is_proper_pair))
                c.execute("""insert into %s values
                (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""" % self.pairs_table, values)
                
            elif not all_mapped:
                # mate with other mate unmapped
                unmapped = [i for i, r in enumerate(readset) if r.is_unmapped][0]
                mapped = [i for i, r in enumerate(readset) if not r.is_unmapped]
                if len(mapped) == 0:
                    self.unmapped += 1
                    continue
                else:
                    mapped = mapped[0]
                
                reads = [readset[unmapped], readset[mapped]]
                strands = [("forward", "reverse")[int(r.is_reverse)] for r in reads]
            
                # reverse complement reads *back* if they are unmapped and reverse
                if (reads[0].is_reverse):
                    unmapped_seq = reverse_complement(reads[0].seq)
                else:
                    unmapped_seq = reads[0].seq

                values = (self.barcode, reads[0].qname,
                          unmapped_seq, reads[0].mapq, strands[0],
                          self.samfile.getrname(reads[1].rname),
                          reads[1].seq, reads[1].pos, reads[1].mapq, strands[1])
                c.execute("""insert into %s values
                (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""" % self.unmapped_mates_table, values)
            else:
                print "warning: unhandled case"
        self.conn.commit()

                 
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
            else:
                self.unmapped += 1
            self.add_pair(read)              



if __name__ == "__main__":
    parser = OptionParser()
    options, args = parser.parse_args()

    if len(args) < 1:
        parser.error("SAM file argument required.")

    if not os.path.exists(args[0]):
        parser.error("File in argument does not exist.")

    dbname = os.path.basename(args[0]).split('.')[0] + '.db'
    m = PairedReads(args[0], dbname)
    m.gather_pairs()
    m.process_pairs()
