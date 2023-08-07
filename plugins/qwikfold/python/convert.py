import os 
import sys
from Bio import SeqIO

filename=sys.argv[1]

SeqIO.convert(filename, "stockholm", '/tmp/alignment.aln', "clustal")

