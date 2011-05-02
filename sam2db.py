"""
sam2sql.py - load aligned reads from a SAM file in a SQLite database.
"""

import os
import pysam
import sqlite3
from optparse import OptionParser
import pdb
from string import maketrans
complement_table = maketrans("ATCGatcg", "TAGCtagc")

def reverse_complement(seq):
    """
    Find the reverse complement of a sequence.
    """
    return seq.translate(complement_table)[::-1]

schema = """
(id integer primary key,
name text,
cigar text,
qname text,
is_duplicate integer,
is_paired integer,
is_proper_pair integer,
is_qcfail integer,
readnum integer,
is_reverse integer,
is_secondary integer,
is_unmapped integer,
isize integer,
mapq integer,
mate_is_reverse integer,
mate_is_unmapped integer,
mrnm text,
pos integer,
qend integer,
qlen integer,
qqual integer,
qstart integer,
qual text,
query text,
rlen text,
tid text,
seq text
)
"""


class AlignedReads(object):
    """
    """
    
    def __init__(self, filename, dbname, name=None):
        """
        
        """
        self.mode = 'r'
        self.filename = filename
        self.dbname = dbname
        if name is None:
            self.name = filename.split('.')[0]
        else:
            self.name = name

        # Detect BAM file
        if filename.split('.')[1] == "bam":
            self.mode = 'rb'

        self.samfile = pysam.Samfile(filename, self.mode)
        self.conn = sqlite3.connect(self.dbname)
        self.total = 0
        self.table = "aligned_reads"

        # Make database/table if it doesn't
        self.create_table(self.table, schema)        

    @staticmethod
    def make_cigar(cigarlist):
        """
        Turn pysam's CIGAR list back to into a string.
        """
        if cigarlist is None:
            return ""
        cigar = "MIDNSHP=X"
        return ''.join([cigar[op] + str(len) for op, len in cigarlist])

    def create_table(self, tablename, schema, verbose=True):
        """
        Check for existence of table; if they don't exist, create.
        """
        check_query = """select name from sqlite_master where type='table' and name=?"""
        c = self.conn.cursor()
        c.execute(check_query, (tablename,))
        result = c.fetchone()
        if result is not None:
            result = result[0]
        if result != tablename:
            if verbose:
                print "creating table '%s'" % tablename
            c.execute("create table %s %s" % (tablename, schema))
            self.conn.commit()

    @staticmethod
    def check_value(value, target_type):
        """
        Check a value from a pysam AlignedReads attribute. Converts
        None values to target type.
        """
        if target_type is str:
            if value is None:
                return ""
            return str(value)
        if target_type is int:
            if value is None:
                return -1 ## pysam attributes shouldn't clash with this
            return int(value)
        if target_type is bool:
            if value is None:
                return -1
            if value not in [0, 1]:
                raise Exception, ("Database takes integer bool type (0, 1); "
                                  "got integer outside this set.")
            return int(value)
        raise Exception, "Could not coerce type to insert into database"

    def getrname(self, tid):
        """
        Get a reference sequence ID from a tid number; this wraps
        pysam's getrname, but with handling of the case when tid is
        -1.
        """
        if tid == -1:
            return ""
        return self.samfile.getrname(tid)

    def read_to_values(self, read):
        """
        Take a read and format into values to insert into table.

        Note: unmapped entries' sequences can be in the SAM file as
        reverse complements, so we need to check the bitflag to see if
        they are, and possibly reverse complement the sequence *back*
        to what it is in the original FASTA/FASTAQ file.
        """
        readnum = 1 if read.is_read1 else 2

        # reverse complement reads *back* if they are unmapped and reverse
        if (read.is_reverse):
            seq = reverse_complement(read.seq)
        else:
            seq = read.seq

        atr_values = [
            (self.name, str),
            (AlignedReads.make_cigar(read.cigar), str),
            (read.qname, str),
            (read.is_duplicate, bool),
            (read.is_paired, bool),
            (read.is_proper_pair, bool),
            (read.is_qcfail, bool),
            (readnum, int),
            (read.is_reverse, bool),
            (read.is_secondary, bool),
            (read.is_unmapped, bool),
            (read.isize, int),
            (read.mapq, int),
            (read.mate_is_reverse, bool),
            (read.mate_is_unmapped, bool),
            (read.mrnm, str),
            (read.pos, int),
            (read.qend, int),
            (read.qlen, int),
            (read.qqual, str),
            (read.qstart, int),
            (read.qual, str),
            (read.query, str),
            (read.rlen, int),
            (self.getrname(read.tid), str),
            (seq, str)]
        values = [AlignedReads.check_value(v, t) for v, t in atr_values]
        return values

    def insert_all(self):
        """
        Iterate over entire SAM file, pushing each entry into
        database.
        """
        c = self.conn.cursor()
        for read in self.samfile:
            self.total += 1
            values = self.read_to_values(read)
            query = ("insert into %s values (null, %s)" %
                     (self.table, ', '.join(['?' for _ in range(len(values))])))
            c.execute(query, values)
        self.conn.commit()


if __name__ == "__main__":
    parser = OptionParser()
    options, args = parser.parse_args()

    if not os.path.exists(args[0]):
        parser.error("Specified SAM file does not exist.")

    a = AlignedReads(args[0], "samdb", os.path.basename(args[0]).split('.')[0])
    a.insert_all()

    # Create index on qname
    c = a.conn.cursor()
    c.execute("create index qname_index on aligned_reads (qname)")
    a.conn.commit()
