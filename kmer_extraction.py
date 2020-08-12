# /usr/bin/env python
"""
kmer extraction: a script to crate clr kmer profles from 4mers
usage: kmer_extraction.py infile outfile ksize

"""

import os
import sys
import csv

from Bio import SeqIO
from Bio.Seq import Seq
import skbio.stats.composition

import vica

def _write_kmers_as_csv(infile, outfile, ksize, kmers, ):
    """Calculate centered log ratio transformed transformed kmer compositions for sequences in a fasta
       file.

   Takes a multi-sequence fasta file and a list of kmers and calculates the
   centered log-ratio transformed kmer composition for each sequence,
   writing a CSV file with the data.

   Args:
       infile (str): a Fasta file
       outfile (str): a path to a CSV output file
       ksize (int): the kmer size, 4-8
       kmers (list): A list of the kmers to count

   Returns:
       None
    """
    # The length id non-redundant kmer vectors for each k
    try:
        with open(infile, 'r') as f1:
            with open(outfile, 'w', buffering=16777216) as csvfile:
                mywriter = csv.writer(csvfile, lineterminator='\n')
                header = ["id"]
                header.extend(kmers)
                recnum = 0
                for record in SeqIO.parse(f1, 'fasta'):
                    rl = [record.id]
                    kmer_frequency = vica.khmer_features.get_composition(ksize=int(ksize),
                                                     seq=str(record.seq).upper(),
                                                     kmers=kmers)
                    kmer_z = skbio.stats.composition.multiplicative_replacement(kmer_frequency)
                    kmer_ilr = skbio.stats.composition.clr(kmer_z)
                    rl.extend(kmer_ilr)
                    mywriter.writerow(rl)
                    recnum += 1
                print("Wrote {} kmer records to {}.".format(recnum, outfile))
    except Exception as err:
        print("Could not write kmer profiles to file")
        raise err


def main(infile, outfile, ksize):
    kmers = vica.khmer_features.itrate_kmer(k=ksize, rc=True)
    _write_kmers_as_csv(infile=infile, outfile=outfile, ksize=ksize, kmers=kmers)

if __name__ == '__main__':
    main(infile=sys.argv[1], outfile=sys.argv[2], ksize=sys.argv[3])
