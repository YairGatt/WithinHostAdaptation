#!/usr/bin/env python

"""
Create tree evolutionary model breseq module, includes all breseq related functions
"""

import sys, os

def parse_breseq_files(breseq_files):
	#parse breseq files to dict
	breseq_dict = {}
	#create dictionary of files
	for gdfile in breseq_files:
		filename = os.path.basename(gdfile).replace(".gd","")
		breseq_dict[filename] = {}
		#open file
		with open(gdfile) as fl:
			content = fl.readlines()
		#iterate through lines and add non-evidence lines
		for line in content:
			#get first word in line
			starter = line.split()[0]
			if len(starter) != 3: continue #evidence lines or whatever
			#form relevant dict
			position = "_".join(line.split()[3:])
			#add to dict
			if position not in breseq_dict[filename]: breseq_dict[filename][position] = None
			else: raise Exception("Position %s appears twice in breseq file %s" % (position,gdfile))				
	#return
	return breseq_dict


def breseq_filter(Samples, breseq_files):
	"""
	Only keep samples with breseq files
	"""
	#initialize results vector
	final_samples = []
	#get the names of all breseq files
	breseq_names = []
	for gdfile in breseq_files:
		filename = os.path.basename(gdfile).replace(".gd","")
		breseq_names.append(filename)
	for sample in Samples:
		if sample.name in breseq_names: final_samples.append(sample)
	#return
	return final_samples

def rank_pairs_breseq(pairs, breseq_files):
	"""
	This function ranks all the pairs according to the number of changes between them. It returns a list in the same order of pairs with the number of evolutionary events in the pair
	"""
	#initialize
	pairs_values = []
	#create breseq dict
	breseq_dict = parse_breseq_files(breseq_files)
	#iterate through pairs and peform comparisons
	for pair in pairs:
		if pair[0] == "fiction":
			pairs_values.append(10000000) #add very large value, it's the only option anyway, but if there's anything else, prefer it
			continue
		#define the progenitor dict
		try: progenitor = breseq_dict[pair[0]]
		except KeyError:
			breseq_dict[pair[0]] = form_internode(pair[0],breseq_dict)
			progenitor = breseq_dict[pair[0]]
		#define progebny dict
		try: progeny = breseq_dict[pair[1]]
		except KeyError:
			breseq_dict[pair[1]] = form_internode(pair[1],breseq_dict)
			progeny = breseq_dict[pair[1]]
		#diff is the difference between the two breseq files
		diff = compare_breseqs(progenitor, progeny)
		#print "difference:", progenitor, progeny, diff
		#append to list
		pairs_values.append(diff)
	#return
	return pairs_values

def compare_breseqs(breseq_dict1,breseq_dict2):
	#see how many lines defer between two breseq files
	#intiialize
	diff = 0
	#get diffs from breseq1 not in breseq2
	for i in breseq_dict1:
		if i not in breseq_dict2: diff += 1
	#get diffs from breseq2 not in breseq1
	for i in breseq_dict2:
		if i not in breseq_dict1: diff += 1
	#return
	return diff

def form_internode(pair, breseq_dict):
	"""
	Create an internode breseq
	"""
	#crete dict
	pair_dict = {}
	#define samples
	sample1 = pair.split("|")[0]
	sample2 = pair.split("|")[1]
	#intersect dicts
	for i in breseq_dict[sample1]:
		#add intersection to internode dict
		if i in breseq_dict[sample2]: pair_dict[i] = breseq_dict[sample1][i]
	#return
	return pair_dict

def compare_breseq_lines(line1,line2):
	#see if two breseq lines point at same evolutionary event
	if "\t".join(line1.split("\t")[3:]) == "\t".join(line2.split("\t")[3:]):
		return True
	else:
		return False

if __name__ == "__main__":
		exit(main())
