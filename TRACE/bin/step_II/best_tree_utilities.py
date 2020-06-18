#!/usr/bin/env python

"""
Best tree utilities
"""

import sys, os
from single_tree_utilities import get_internodes, Tree, Node
from evolutionary_utilities import Sample, strain_converted, separate, plot_histogram
from combinatorics_utilities import sublisting
import warnings
from tree_utilities import create_all_trees

def create_best_tree(reverse_pairs_dictionary, strain, converted):
	"""
	Calculate theoretical best tree from pairs
	"""
	#final tree
	final_tree = Tree()
	#call best origin for each timepoint
	#define samples
	points = reverse_pairs_dictionary.keys()
	#get timepoints
	timepoints = {}
	for i in points:
		#get time id
		try: time_id = strain_converted(i, strain, converted)
		except ValueError: time_id = strain_converted(i.split("|")[0], strain, converted)
		#get timepoint
		try: timepoints[separate(time_id,"number")].append(i)
		except KeyError: timepoints[separate(time_id,"number")] = [i]
	#initiaalize
	finished_timepoints = []
	#iterate through samples and get best for easy timepoints
	for timepoint in timepoints:
		#simple timepoint with one sample
		if len(timepoints[timepoint]) == 1:
			progeny = timepoints[timepoint][0]
			#determine values of progenitor options
			value, node = call_value(reverse_pairs_dictionary, progeny)
			#add to final true
			final_tree.nodes.append(node)
		#complicated timepoint with multiple samples
		else:
			#considere all subtrees for timepoint and find the one with the best value
			#get_all_subgroups
			subdivide = [i for i in timepoints[timepoint] if "|" not in i]
			subdivisions = sublisting(subdivide)
			#get the best value of the subgroup and the leader of each subgroup
			values = []
			subdivisions_nodes = []
			#iterate through subdivisions
			for subdivision in subdivisions:
				nodes, value = values_for_subgroup(reverse_pairs_dictionary,subdivision)
				#add to lists
				subdivisions_nodes.append(nodes)
				values.append(value)
			#final value is value+leader best value
			#determine best subgroup division
			value = min(values)
			best_subdivision = subdivisions_nodes[values.index(min(values))]
			#add to final tree
			final_tree.nodes += best_subdivision
		#make sure not to do this timepoint again
		finished_timepoints.append(timepoint)
	#return
	return final_tree

def call_value(pairs_dictionary,progeny):
	#determine values of progenitor options
	values = pairs_dictionary[progeny].values()
	#if multilple options are tied for best there might be a problem
	if values.count(min(values)) > 1: warnings.warn("Multiple best prognitors for %s" % progeny)
	#find the best option
	for progenitor in pairs_dictionary[progeny]:
		if pairs_dictionary[progeny][progenitor] == min(values): node = Node(progenitor,progeny)
	#return
	return min(values), node

def values_for_subgroup(reverse_pairs_dictionary,subdivision):
	"""
	Get the best value of the subgroup
	"""
	#created limited reverse_pairs_dictionary
	limited_reverse_pairs_dictionary = {}
	for i in reverse_pairs_dictionary:
		limited_reverse_pairs_dictionary[i] = {}
		for j in reverse_pairs_dictionary[i]:
			if "|" not in j: limited_reverse_pairs_dictionary[i][j] = reverse_pairs_dictionary[i][j]
	#initialize
	value = 0
	subtree = []
	complicated = []
	#iterate through sublists in the division
	for sublist in subdivision:
		#if simple
		if len(sublist) == 1:
			progeny = sublist[0]
			local, node = call_value(limited_reverse_pairs_dictionary, progeny)
			#add value to overall value
			value += local
			#add node to tree
			subtree.append(node)
		#if complicated
		else:
			#must create all subtrees and determine their values, along with the best value for their leader
			Samples = [Sample(i,"0a") for i in sublist]
			subsubtrees = create_all_trees(Samples)
			#get value of subtrees
			values = []
			for subsubtree in subsubtrees:
				subvalue = 0
				for node in subsubtree.nodes:
					#progenitor
					if type(node.progenitor) == str:
						progenitor = node.progenitor
						if "|" in progenitor: progenitor = "%s|%s" % (min(progenitor.split("|")),max(progenitor.split("|")))
					else: progenitor = node.progenitor.name
					#progeny
					if type(node.progeny) == str:
						progeny = node.progeny
						if "|" in progeny: progeny = "%s|%s" % (min(progeny.split("|")),max(progeny.split("|")))
					else: progeny = node.progeny.name
					#check values
					if progeny not in reverse_pairs_dictionary: print "progeny has no progenitor %s" % progeny
					if progenitor not in reverse_pairs_dictionary[progeny] and progenitor != "fiction": print "progenitor %s not possible for progeny %s" % (progenitor, progeny)
					if progenitor == "fiction": local = 10000000
					else: local = reverse_pairs_dictionary[progeny][progenitor]
					#add to value of subsubtree
					subvalue += local
				#get root
				root = subsubtree.get_root_simple()
				if len(root) > 1:
					warnings.warn("Skipping potential tree for multiple roots")
					continue
				#get value of root edges
				progeny = root[0]
				local, node = call_value(limited_reverse_pairs_dictionary, progeny)
				#add value to overall value
				subvalue += local #final value of subsubtree
				subvalue -= 10000000
				#add to final list
				values.append(subvalue)
				subsubtree.nodes.append(node)
			#decide on the very best subsubtree
			value += min(values)
			best_subsubtree = subsubtrees[values.index(min(values))]
			subtree += best_subsubtree.nodes
	#return
	return subtree, value

def get_trees(pair, sublist, count):
	"""
	Get all subtrees of a timepoint based on pair and sublist
	"""
	if count == -1:
		return [[pair]], list(pair)
	else:
		new_trees = []
		node = sublist[count]
		curr_trees, mems = get_trees(pair,sublist,count-1)
		#skip if node already present
		if node in mems: return curr_trees, mems #skip
		#iterate through trees
		for tree in curr_trees:
			#add all new possible pairs to all trees
			for i in mems:
				new_node = (i,node)
				new_trees.append(tree + [new_node])
				
		return new_trees, mems + [node]

if __name__ == "__main__":
		exit(main())