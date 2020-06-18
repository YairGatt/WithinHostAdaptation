"""
This script add genes that had no mutations to the summary.genes file so it is representative of the whole proteome
"""

import os,sys
import gzip as gz
from Bio import SeqIO
from subprocess import call

genes_summary_file = sys.argv[1]
proteome = sys.argv[2]

outfile = genes_summary_file.split(".txt")[0] + "_additional.txt"
#outdir = os.path.dirname(genes_summary_file)
outfile_final = genes_summary_file.split(".txt")[0] + "-all.txt"

#get all proteins
with gz.open(proteome,"rb") as fl:
    proteins = list(SeqIO.parse(fl,"fasta"))
protein_names = [i.id for i in proteins]

#get all proteins in genes summary file
with open(genes_summary_file) as fl:
    content = fl.readlines()
proteins_present = [i.split()[2] for i in content if not i.startswith("#")]
length = len(content[1].split())

proteins_to_add = []
for i in protein_names:
    if i not in proteins_present:
        proteins_to_add.append(i)

lines = []
for protein in proteins_to_add:
    #add 4 first elemetns: X and X and genename and geneID, protein as transcript name and nothing as BioType
    protein_line = "X\tX\t%s\t" % protein
    for i in xrange(length-4):
        protein_line += "\t0"
    protein_line += "\n"
    lines.append(protein_line)

with open(outfile,"w") as outfl:
    for line in lines:
        outfl.write(line)

call("cat %s %s > %s" % (genes_summary_file, outfile, outfile_final), shell=True)
