#!/usr/bin/env python

"""
Utilities for creating trees for the create_tree_evolutionary_module script
"""

import sys, os
import random
import itertools
from single_tree_utilities import create_tree
from evolutionary_utilities import Sample

def draw_random_trees(list_of_lists,Samples,index,n=100000):
	"""
	Draw n random trees
	"""
	random_trees = []
	for i in xrange(n):
		#create raw_tree
		raw_tree = []
		for alist in list_of_lists:
			origin = random.choice(alist)	
			raw_tree.append(origin)
		#create tree
		print [i.name for i in raw_tree], [i.name for i in Samples]
		final_trees = create_tree(raw_tree, Samples, index)
		if final_trees: random_trees += final_trees
	#return
	return random_trees

def create_all_trees(Samples):
	"""
	This function creates all possible trees from the samples included in the strain
	"""
	trees = []
	#add theoretical sample for any cases where multiple samples exist from same timepoint
	#iterate to define options for progenitors for each sample
	times = [sample.time for sample in Samples]
	for sample1 in Samples:
		for sample0 in Samples:
			if sample0.time <= sample1.time and sample0.name != sample1.name: sample1.options.append(sample0)
		#add fictional node as root if there are two time0 samples, otherwise time0 sample is the root
		if sample1.time == min(times) and times.count(min(times)) > 1: sample1.options.append("fiction")
	#build combined list of all possible origins
	list_of_lists = []
	for i in Samples: list_of_lists.append(i.options)
	#all possible combinations
	if [] in list_of_lists: #time0 sample, root
		index = list_of_lists.index([])
		list_of_lists.pop(index)
	else:
		index = None
	#if there are two []s
	if [] in list_of_lists: raise Exception("Two samples with no origin, investigate further") #this should be solved due to the fictional origin
	#create raw trees with all combinations of options for the different samples
	if len(Samples) < 10:
		#raw tree - each element in the list is the origin of the respective element in samples, this defines a tree
		trees_raw = list(itertools.product(*list_of_lists))
		#go through the raw trees and change their format to pairs of each sample and its origin
		for tree in trees_raw:
			#create trees
			final_trees = create_tree(tree, Samples, index)
			#add to final list
                        #if final_tree:
                            #for node in final_tree.nodes:
                                #print node.progenitor, node.progeny
			if final_trees: trees += final_trees
                        #if final_tree: trees += [final_tree]

	else: #if there are too many trees, use random trees instead of all trees
		trees = draw_random_trees(list_of_lists, Samples, index, 1000000)
	#return
	return trees

if __name__ == "__main__":
		exit(main())
