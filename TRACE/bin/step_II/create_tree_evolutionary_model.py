#!/usr/bin/env python

"""
This script recieves the assemblies of multiple bacterial samples that were isolated from the same patient, and were determined to be derived from the same strain using previous scripts.
It then uses the breseq results/kSNP-produced SNP matrix for those samples to determine the most likely evolutionary model.
This is done by:
1. Assessing the number of evolutionary events required between any two isolates of different time points (IE how many lines of differences from reference strain are not identical)
2. Creating all possible trees and ranking them according to the previous assessment. The tree with the fewest events (highest parsimony) is selected as the most likely evolutionary model,
3. The assessment of how likely each branch is - in how many alternative likely trees did this branch not occur?

Creation of the trees:
The algorithm depends on the time points of the samples to create an evolutionary path, proxy samples are added for time points with multiple samples
"""

import sys, os
import argparse
import math
import numpy as np
from evolutionary_utilities import mkdir, samples_to_objects, convert_to_dict, plot_histogram
from breseq_utilities import rank_pairs_breseq, breseq_filter
from kSNP_utilities import rank_pairs_matrix
from tree_utilities import create_all_trees
from best_tree_utilities import create_best_tree
import warnings

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
		'-w', '--workdir', default="./temp/",
		help='Workdir where results will be written.')

	parser.add_argument(
		'-s', '--strain', default="./strain_1.txt",
		help='Strain file with times of the samples in relevant strain..')

	parser.add_argument(
		'-m', '--matrix',
		help='Matrix file created by kSNP in fasta format.')

	parser.add_argument(
		'-p', '--phylip',
		help='Phylip fasta file created from kSNP matrix.')

	parser.add_argument(
		'-b', '--breseq_files', nargs="+",
		help='Workdir where results will be written.')

	parser.add_argument(
		'-c', '--converted', default="./converted_1.txt",
		help='Converted file matching the accessions of the different samples in the strain.')

	parser.add_argument(
		'-t', '--threshold', type=float, default=0.4,
		help='Converted file matching the accessions of the different samples in the strain.')

	parser.add_argument(# customized description; put --help last
		'-h', '--help', action='help',
		help='Show this help message and exit.')

	settings = parser.parse_args(argv)

	return settings

def all_pairs(trees):
	"""
	This function recieves a list of all samples in the strain. It then uses the associated time points to create a list of all unique progenitor-progeny pairs.
	"""
	#initialize
	pairs = []
	count = 0
	#iterate through trees
	for tree in trees:
		#iterate through nodes
		for item in tree.nodes:
			#get progenitro
			if type(item.progenitor) == str: progenitor = item.progenitor
			else: progenitor = item.progenitor.name
			#get progeny
			if type(item.progeny) == str: progeny = item.progeny
			else: progeny = item.progeny.name
			#add to pairs
			pairs.append((progenitor,progeny))
	#get unique pairs
	pairs = list(set(pairs))
	#return pairs
	return pairs

def rank_all_trees(trees, pairs_dictionary, percent):
	"""
	This function ranks all possible trees according to the distances between the pairs that make them up, and determines the top tree (minimum evolutionary events)
	"""
	#intialize
	ranked_trees = []
	#iterate through trees
	for tree in trees:
		#initialize rank as 0
		tree.rank = 0
		#iterate through nodes
		for node in tree.nodes:
			#define progenitor and progeny in dict format
			try: progenitor = node.progenitor.name
			except AttributeError:
				if node.progenitor == "fiction": progenitor = "fiction"
				else: progenitor = "%s|%s" % (min(node.progenitor.split("|")),max(node.progenitor.split("|")))
			#progeny
			try: progeny = node.progeny.name
			except AttributeError: progeny = "%s|%s" % (min(node.progeny.split("|")),max(node.progeny.split("|")))
			#find value
			node_value = pairs_dictionary[progenitor][progeny]
			#add to tree rank
			tree.rank += node_value
		#add value to list of values
		ranked_trees.append(tree.rank)
	#find tree with fewest evolutionary events
	minimum_value = min(ranked_trees)
	#get very best trees
	top_trees = [i for i in trees if i.rank == minimum_value]
	#sort trees
	sorted_trees = []
	#iterate
	for n,tree in enumerate(top_trees):
		#initialize sorted
		sorted_tree = []
		#iterate through nodes
		for node in tree.nodes:
			#progenitor
			try: progenitor = node.progenitor.name
			except AttributeError:
				if node.progenitor == "fiction": progenitor = "fiction"
				else: progenitor = "%s|%s" % (min(node.progenitor.split("|")),max(node.progenitor.split("|")))
			#progeny
			try: progeny = node.progeny.name
			except AttributeError: progeny = "%s|%s" % (min(node.progeny.split("|")),max(node.progeny.split("|")))
			#add to sorted
			sorted_tree.append((progenitor,progeny))
		#sort
		sorted_tree.sort()
		#add to sorted trees if not  found yet
		if sorted_tree not in sorted_trees: sorted_trees.append(sorted_tree)
		else: top_trees[n] = None #otherwise remove tree
	#remove removed trees
	while None in top_trees:
		top_trees.remove(None)
	#see if there is actually more than one top tree
	if len(top_trees) > 1: warnings.warn("More than one top tree! Inspect further")
	#get top percentage of trees
	top_percentage = int(math.ceil(float(len(trees))*(float(percent)/100)))
	top_percentage_trees = [tree for rank,tree in sorted(zip(ranked_trees,trees))][:top_percentage]
	if percent == 100: top_percentage_trees = top_trees
	#return top tree  and top percent% trees
	return top_trees, top_percentage_trees

def support_for_nodes(top_trees, top_percentage_trees, threshold):
	"""
	Calculate the percentage of the top percent% trees that include the node - a proxy for the support of the node
	This function also chooses the tree with the best average support if a few top trees are available, if there are ties the first will be chosen arbitrarily
	"""
	#initialize vectors
	support = []
	supported_pairs = []
	#iterate through nodes
	for top_tree in top_trees:
		support_tree = []
		supported_pairs_tree = []
		for node1 in top_tree.nodes:
			#initialize supporting trees
			supporting_trees = 0
			#iterate through top 1% trees
			for tree in top_percentage_trees:
				for node2 in tree.nodes:
					if node1.progenitor == node2.progenitor and node1.progeny == node2.progeny:
						supporting_trees += 1
						break
			#normalize supporting trees by number of trees in top percentage
			normalized_support = float(supporting_trees)/len(top_percentage_trees)
			support_tree.append(normalized_support)
			#add pair to suppoerted pairs
			if normalized_support >= threshold: supported_pairs_tree.append((node1.progenitor,node1.progeny))
		#add to vector of vectors
		support.append(support_tree)
		supported_pairs.append(supported_pairs_tree)
	#call between trees
	measure = [np.mean(i) for i in support]
	top_index = measure.index(max(measure))
	#return
	return support[top_index], supported_pairs[top_index]

def write_to_file(supported_pairs, workdir, breseq=False):
	"""
	Write results to a comparisons converted file
	"""
	#determine name
	if breseq: outname = "breseq_comparisons_converted.txt"
	else: outname = "comparisons_converted.txt"
	#open file
	with open(workdir + "/" + outname,"w") as outfl:
		#add all pairs
		for pair in supported_pairs:
			#define progenitor
			progenitor = pair[0]
			if type(progenitor) != str: progenitor = progenitor.name
			#progeny
			progeny = pair[1]
			if type(progeny) != str: progeny = progeny.name
			#create line
			line = "%s\t%s\n" % (progenitor, progeny)
			#write
			outfl.write(line)

def main(argv=None):
	#process command line
	settings = process_command_line(argv)
	#create workdir
	mkdir(settings.workdir)
	#define all samples as sample objects
	Samples = samples_to_objects(settings.strain, settings.converted)
	#filter out samples with no breseq files if breseq mode
	if settings.breseq_files: Samples = breseq_filter(Samples, settings.breseq_files)
	#create all trees
	trees = create_all_trees(Samples)
	#inspect all pairs
	pairs = all_pairs(trees)
	#if breseq data available use it to assess differences between all pairs
	if settings.breseq_files: pairs_values = rank_pairs_breseq(pairs, settings.breseq_files)
	#otherwise use a kSNP matrix
	elif settings.matrix: pairs_values = rank_pairs_matrix(pairs, settings.matrix, settings.strain, settings.converted)
	#if the kSNP matrix is only available in phylip fasta format convert it and use it
	elif settings.phylip: pairs_values = rank_pairs_matrix(pairs, settings.phylip, settings.strain, settings.converted, mode="phylip")
	else: raise Exception("Script must use either breseq files or a kSNP matrix in fasta form or phylip form")
	#convert to dict
	pairs_dictionary, reverse_pairs_dictionary = convert_to_dict(pairs,pairs_values)
	#if trees were drawn at random, add very best theoretical trees
	if len(Samples) >= 10:
		best_tree = create_best_tree(reverse_pairs_dictionary, settings.strain, settings.converted)
		#add to trees
		trees.append(best_tree)
	#rank all trees and return trees with best value and a top percentage of the trees
	#percent determines what percentage of the top trees to return
	if len(trees) >= 10000: percent = 1
	elif len(trees) >= 100: percent = 10
	else: percent = 100
	#get top trees
	top_trees, top_percentage_trees = rank_all_trees(trees, pairs_dictionary, percent)
	#get strength of evidence for each pair and pairs with support above threshold
	support,supported_pairs = support_for_nodes(top_trees, top_percentage_trees, settings.threshold)
	#write to file
	write_to_file(supported_pairs, settings.workdir, settings.breseq_files)

if __name__ == "__main__":
		exit(main())
