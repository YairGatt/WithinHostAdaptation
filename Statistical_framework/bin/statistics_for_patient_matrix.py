#!/usr/bin/env python

"""
Statistics for patient_matrix
"""
#general imports
import sys, os
#statistical tools
from statsmodels.stats.multitest import multipletests
from scipy import special, stats
import random

def p_values_by_simulation(probabilities, gene_vector_dict, NUM_SIMULATIONS=1000000):
    """
    Get p-values by simulations NUM_SIMULATIONS genes and seeing in how many strains they are lost by chance alone
    """
    #intiialize background_distribution
    background_distribution = []
    #carry out simulations
    for i in range(NUM_SIMULATIONS):
        #initialize number of strains gene was lost in
        count = 0
        #iterate thrrough strains
        for prob in probabilities:
            #random
            rand = random.uniform(0,1)
            #was teh gene lost in the strain?
            if rand < prob: count += 1
        #add to distribution
        background_distribution.append(count)
    #initialize p values
    p_value_dict = {}
    #iterate through actual genes
    for gene in gene_vector_dict:
        #see how man strains the gene was lost in
        number = sum(gene_vector_dict[gene])
        #get p-value
        if number >= 1: p_value_dict[gene] = float(sum([1 for i in background_distribution if i >= int(number)])) / NUM_SIMULATIONS
        elif number >= 0 and number < 1: p_value_dict[gene] = 1.0 #if the number is between 0 and 1 we give it a 1.0 v-balue
    #return results
    return p_value_dict

def correct_multiple_hypotheses(p_value_dict, alpha):
    """
    Correct p-values for multiple hypotheses using the Benjamini-Hochberg method
    """
    #initialize adjusted p values dict
    padjust_dict = {}
    #run multiple hypothesis correction in python using the Benjamini-Hochberg method
    truths, padjust, correctedp1, correctedp2 = multipletests(p_value_dict.values(), alpha*2, method = 'fdr_bh')
    #define p_adjust dict
    for index,gene in enumerate(p_value_dict.keys()):
        padjust_dict[gene] = padjust[index]
    #return
    return padjust_dict

def main(argv=None):
    raise Exception("This should be called from statistical_framework.py")

if __name__ == "__main__":
        exit(main())
