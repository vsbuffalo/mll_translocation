"""
divide-singles.py

divide-singles.py takes the output from a *-singles.txt file (from
find-odd-mates.py), and seperates by each mapped mate's
location/chromosom into FASTA files
"""

from optparse import OptionParser
import os
import sys

def build_location_dict(filename, outdir):
    """
    Given a *-singles.txt, build a dictionary of all unmapped
    sequences by the mapped mate's position, and write each location
    to a seperate FASTA file, to be aligned to the MLL template using
    BWA bwasw.
    """
    location_dict = dict()

    for line in open(filename):
        (qname, unmapped_seq, unmapped_pos,
         mapped_location, mapped_seq, mapped_pos) = line.split('\t')

        result = location_dict.get(mapped_location, False)

        if not result:
            location_dict[mapped_location] = dict()
        assert(location_dict[mapped_location].get(qname, False) is False)
        location_dict[mapped_location][qname] = unmapped_seq

    # Write FASTA files
    basename = os.path.basename(filename).split('-')[0]
    for location in location_dict:
        filename = os.path.join(outdir, basename + '-' + location + '-unmapped.fasta')
        with open(filename, 'w') as f:
            for header, seq in location_dict[location].items():
                f.write(">%s\n%s\n" % (header, seq))
    
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="dir",
                      help="directory to output files (default: name of first arg)",
                      default=None)
    (options, args) = parser.parse_args()

    basename = os.path.basename(args[0])
    basename = basename.split('.')[0]
    if options.dir is None:
        options.dir = basename + "-unmapped"
    if not os.path.exists(options.dir):
        os.mkdir(options.dir)

    if len(args) != 1:
        parser.error("Incorrect number of arguments.")
    elif not os.path.exists(args[0]):
        parser.error("File '%s' does not exist." % args[0])
    
    build_location_dict(args[0], options.dir)
