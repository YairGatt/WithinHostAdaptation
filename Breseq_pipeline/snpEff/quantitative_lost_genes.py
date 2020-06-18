"""
This script assess which genes were actually lost in progeny strain out of available genes in progenitor strain based on snpEff results
"""
import os,sys
#input
lost_genes_file = sys.argv[1]
background_file = sys.argv[2]
output_lost_genes_file = sys.argv[3]
#scan files to vectors
with open(lost_genes_file) as fl: lost_genes = [i.strip() for i in fl.readlines()]
with open(background_file) as fl: background = [i.strip() for i in fl.readlines()]
#count for each gene how many times it had been lsot out of the background
genes = []
#for gene in lost_genes:
for gene in background:
    gene_background = background.count(gene)
    gene_lost = float(lost_genes.count(gene))
    try: genes.append((gene,gene_background, gene_lost/gene_background))
    except ZeroDivisionError: continue
#output format
#write to output file
with open(output_lost_genes_file,"w") as outfl:
    for gene in genes:
        line = "%s\t%s\t%s\n" % gene
        outfl.write(line)
