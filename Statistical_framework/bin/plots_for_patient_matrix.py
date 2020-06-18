#!/usr/bin/env python

"""
Plots for patient_matrix
"""

import sys, os
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from itertools import combinations
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def plotting(gene_vector_dict, outfile, number_of_patients, threshold=0, log=False):
        """
        Plot figures with distribution of number of patients genes were lost in
        """
        #values for histogram
        for_hist = []
        for gene in gene_vector_dict:
                number = sum(gene_vector_dict[gene])
                for_hist.append(number)
        #plot histogram with distribution
        dictionary = plt.figure(figsize=(20,20))
        plt.title("Distribution of the number of genes lost in specific numbers of patient, among %s patients" % numbents, fontsize=25)
        plt.xlabel("Number of patients the gene was lost in", fontsize=20)
        if log: plt.ylabel("Number of genes (log)",fontsize=20)
        else: plt.ylabel("Number of genes",fontsize=20)
        plt.tick_params(labelsize=20)
        #leave only data over thresholds
        for_hist_above_threshold = [i for i in for_hist if i>=threshold]
        data = np.array(for_hist_above_threshold)
        try: d = np.diff(np.unique(data)).min()
        except ValueError: d = 0.5
        left_of_first_bin = data.min() - float(d)/2
        right_of_last_bin = data.max() + float(d)/2
        if log: plt.hist(data, np.arange(left_of_first_bin, right_of_last_bin + d, d), log=True)
        else: plt.hist(data, np.arange(left_of_first_bin, right_of_last_bin + d, d))
        directory = os.path.dirname(outfile)
        if directory == "": directory = "."
        if log: log_indicator = "_log"
        else: log_indicator = ""
        #add fitting curve
        if threshold > 0 and not log:
                xdata = []
                ydata = []
                for i in xrange(threshold, max(for_hist_above_threshold)+1):
                        xdata.append(i)
                        ydata.append(for_hist_above_threshold.count(i))
                xdata = np.array(xdata)
                ydata = np.array(ydata)

                popt, pcov = curve_fit(lambda x,gamma: (xdata+1)**-gamma*len(for_hist), xdata, ydata)
                plt.plot(xdata, (xdata+1)**-popt*len(for_hist), 'r-', label='fit: gamma=%5.3f' % tuple(popt))
        #save
        plt.savefig("%s_lost_distribution_above_%s_patients%s.png" % (".".join(outfile.split(".")[:-1]), threshold, lor))

def compute_jaccard_index(set_1, set_2):
        set_1 = set(set_1)
        set_2 = set(set_2)
        n = len(set_1.intersection(set_2))
        try: return n / float(len(set_1) + len(set_2) - n)
        except ZeroDivisionError: return 1

def clustergram(patients_lost_genes):
        #cluster patients before continuing
        raw_linkage = [1-compute_jaccard_index(i[0],i[1]) for i in list(combinations([i[0] for i in patients_lost_gen)],2))]
        #raw_linkage = [1/compute_jaccard_index(i[0],i[1])for i in list(combinations(patient_dict.values(),2))]
        Z = linkage(raw_linkage, 'average')
        c,coph = cophenet(Z,raw_linkage)
        # calculate full dendrogram
        dictionary2 = plt.figure(figsize=(25, 22))
        plt.title('Hierarchical Clustering Dendrogram', fontsize=20)
        plt.xlabel('Patient', fontsize=20)
        plt.ylabel('Distance (1-Jaccard index)', fontsize=20)
        dendrogram(Z, leaf_rotation=90., leaf_font_size=16., labels=patients_lost_genes.keys())
        #dendrogram(Z, truncate_mode='lastp', p=12, leaf_rotation=90., leaf_font_size=12., show_contracted=True)
        cluster_outfile = "_cluster.".join("test_stuff.pdf".rsplit(".",1))
        plt.savefig(cluster_outfile)

def raw_prob(n, probabilities):
        raw_probabilities = []
        value = 1
        for num in xrange(n-1):
                value = faster_iterative_multiplication(probabilities,num+1)
                raw_probabilities.append((num+1, value))
                #append final probabilities
        raw_probabilities.append((n,reduce(lambda x, y: x*y, probabilities)))
        #return
        return raw_probabilities

def main():
        print "This should be called from patient_matrix.py"

if __name__ == "__main__":
                exit(main())
