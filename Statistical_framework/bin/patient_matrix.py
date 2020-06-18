#!/usr/bin/env python

"""
This script creates a patient matrix from lists of genes who were lost in every patient and runs a statistical test.
"""

#general imports
import os, sys
import argparse
import gzip as gz
import csv
#biopthon
from Bio import SeqIO
#specialized functions
import statistics_for_patient_matrix
import plots_for_patient_matrix

def process_command_line(argv):
    """
    Return an args list
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Process input.', add_help=False)

    #define options here

    parser.add_argument(
        '-l', '--lost_genes_files', nargs='+',
        help='Lists of genes that were lost from each patient.')

    parser.add_argument(
        '-o', '--outfile', default="./temp.txt",
        help='Output file.')

    parser.add_argument(
        '-p', '--reference_proteome',
        help='Reference proteome to inspect NCBI genes from')

    parser.add_argument(
        '-i', '--millstone_input_files', nargs='+',
        help='Millstone files detailing strains and timepoints.')

    parser.add_argument(
        '-a', '--alpha', type=float, default=0.05,
        help='Alpha value for statistical test.')

    parser.add_argument(# customized description; put --help last
        '-h', '--help', action='help',
        help='Show this help message and exit.')

    settings = parser.parse_args(argv)

    return settings

def write_to_files(outfile, patients, gene_vector_dict, p_value_dict, padjust_dict, alpha):
    """
    Write results to csv files
    """
    #write all genes vectors
    with open (outfile,"wb") as fl:
        #initialize writer
        writer = csv.writer(fl, delimiter='\t')
        #write header
        header = ['Genes/Patients']+patients
        writer.writerow(header)
        #write rows
        for gene in gene_vector_dict:
            row = [gene]+gene_vector_dict[gene]
            writer.writerow(row)
    #write genes to files
    with open(".".join(outfile.split(".")[:-1]) + "_all_genes.txt","w") as fl_all:
        with open(".".join(outfile.split(".")[:-1]) + "_significant_genes_corrected.txt","w") as fl:
            with open(".".join(outfile.split(".")[:-1]) + "_significant_genes.txt","w") as fl_uncorrected:
                #iterate through genes
                for gene in gene_vector_dict:
                    #define number of patients gene was lost in
                    number = sum(gene_vector_dict[gene])
                    #write all genes to all file
                    fl_all.write("%s\t%s\n" % (gene,number))
                    #if gene above threshold write to significant file
                    if padjust_dict[gene] <= alpha*2: fl.write("%s\t%s\n" % (gene,number))
                    #if gene above threshold write to significant file
                    if p_value_dict[gene] <= alpha: fl_uncorrected.write("%s\t%s\n" % (gene,number))

def gene_loss_for_patients(genes_files):
    #define patient vector
    patients = []
    patients_lost_genes = {}
    #iterate through patient files
    for patient_file in genes_files:
        #decide if patient has multiple isolates
        assembly_path = os.path.dirname(patient_file)+"/assemblies"
        if len(os.listdir(assembly_path)) < 2: continue
        #define patient name
        patient = os.path.basename(os.path.dirname(patient_file))
        patients.append(patient)
        #define genes lost in each patient
        with open(patient_file) as fl: lines = [i.strip() for i in fl.readlines()]
        #see which genes were lost for patinet
        lost_genes = []
        lost_percentage = []
        for line in lines:
            split_line = line.split("\t")
            name = split_line[0].replace(" ", "_").replace(",","_")
            percentage = split_line[2]
            #append
            lost_genes.append(name)
            lost_percentage.append(float(percentage))
        #add to lost genes dictionary
        patients_lost_genes[patient] = (lost_genes,lost_percentage)
    #return results
    return patients,patients_lost_genes

def define_gene_vector_dict(genes, patients, patients_lost_genes):
    #for each gene, define patients where they were lost by creating a gene vector with the order of the patients and Y for a patient where it was lost and N for a patient where it was not lost
    gene_vector_dict = {}
    for gene in genes:
        gene_vector = []
        for patient in patients:
            if gene in patients_lost_genes[patient][0]:
                #append fraction
                index = patients_lost_genes[patient][0].index(gene)
                gene_vector.append(patients_lost_genes[patient][1][index])
            #if the gene was not lost append 0
            else: gene_vector.append(0)
        #add to dict
        gene_vector_dict[gene] = gene_vector
    #return results
    return gene_vector_dict

def get_probabilities(patients_lost_genes, genes):
    #define loss probability for each patient
    probabilities = []
    for patient in patients_lost_genes:
        #calculate probability    for each patient
        probability = float(sum(patients_lost_genes[patient][1]))/len(genes)
        if probability > 1: raise Exception("Loss probability above 1 for patient %s, make sure your lost genes derive from your genome" % patient)
        probabilities.append(probability)
    #return results
    return probabilities

def get_genes(proteome_file):
    #parse cds file to define gene list
    genes = []
    #open file
    if proteome_file[-3:] == ".gz":
        with gz.open(proteome_file) as fl:
            for gene in SeqIO.parse(fl,"fasta"): genes.append(gene.id)
    else:
        with open(proteome_file) as fl:
            for gene in SeqIO.parse(fl,"fasta"): genes.append(gene.id)
    #adjust names
    genes = [i.replace(" ","_").replace(",","_") for i in genes]
    #return
    return genes

def main(argv=None):
    #process command line
    settings = process_command_line(argv)
    #remove repeating lost genes files
    genes_files = list(set(settings.lost_genes_files))
    #parse cds file to define gene list
    genes = get_genes(settings.reference_proteome)
    #define patient vector
    patients, patients_lost_genes = gene_loss_for_patients(genes_files)
    #for each gene, define patients where they were lost by creating a gene vector with the order of the patients and Y for a patient where it was lost and N for a patient where it was not lost
    gene_vector_dict = define_gene_vector_dict(genes, patients, patients_lost_genes)
    #define loss probability for each patient
    probabilities = get_probabilities(patients_lost_genes, genes)
    #get p values by simulation
    p_value_dict = statistics_for_patient_matrix.p_values_by_simulation(probabilities,gene_vector_dict)
    #correct for multiple hypotheses
    padjust_dict = statistics_for_patient_matrix.correct_multiple_hypotheses(p_value_dict, settings.alpha)
    #write rows to file
    write_to_files(settings.outfile, patients, gene_vector_dict, p_value_dict, padjust_dict, settings.alpha)

if __name__ == "__main__":
        exit(main())
