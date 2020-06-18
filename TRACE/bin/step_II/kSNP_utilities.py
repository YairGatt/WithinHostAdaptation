#!/usr/bin/env python

"""
Create tree evolutionary model kSNP module, includes all kSNP related functions
"""

import sys, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from evolutionary_utilities import convert_strain, separate

def parse_phylip(phylip_file):
	"""
	When kSNP matrix not available a phylip-converted matrix can be used. Convert back to fasta like format since SeqIO has trouble handling phylip
	"""
	#initialize
	matrix_dict = {}
	entries = []
	#open file
	with open(phylip_file) as fl:
		content = fl.readlines()
		lines = content[1:]
	#iterate through samples
	for line in lines:
		sample = line.split()[0]
		seq = line.split()[1]
		record = SeqRecord(Seq(seq),id=sample)
		if sample != "reference": entries.append(record)
		else: reference = record
	#create dictionary of files
	for i in entries:
		matrix_dict[i.id] = i.seq
	#return
	return matrix_dict, reference

def parse_matrix(matrix, strain, converted):
	#parse matrix to dict
	matrix_dict = {}
	#open with SeqIO
	with open(matrix,"r") as fl:
		raw_entries = list(SeqIO.parse(matrix,"fasta"))
		entries = [i for i in raw_entries if i.id != "reference"]
		reference = [i for i in raw_entries if i.id == "reference"][0]
	#create dictionary of files
	for i in entries:
		id = i.id
		#convert to sample name
		num = separate(id,"number")
		if num == None: continue #if not format
		#get filename
		try: filename = convert_strain(id, strain, converted)
		except ValueError: continue
		#add sequence
		matrix_dict[filename] = i.seq
	#return
	return matrix_dict, reference

def rank_pairs_matrix(pairs, matrix, strain, converted, mode="matrix"):
	"""
	This function ranks all the pairs according to the number of changes between them. It returns a list in the same order of pairs with the number of evolutionary events in the pair
	"""
	#initialize
	pairs_values = []
	#create matrix dict
	if mode == "matrix": matrix_dict, reference = parse_matrix(matrix, strain, converted)
	elif mode == "phylip": matrix_dict, reference = parse_phylip(matrix)
	#iterate through pairs and peform comparisons
	for pair in pairs:
		if pair[0] == "fiction":
			pairs_values.append(10000000)
			continue
		#define the progenitor dict
		try: progenitor = matrix_dict[pair[0]]
		except KeyError:
			matrix_dict[pair[0]] = form_internode(pair[0],matrix_dict,reference)
			progenitor = matrix_dict[pair[0]]
		#define the progeny dict
		try: progeny = matrix_dict[pair[1]]
		except KeyError:
			matrix_dict[pair[1]] = form_internode(pair[1],matrix_dict,reference)
			progeny = matrix_dict[pair[1]]
		#diff is the difference between the two entries
		diff = compare_entries(progenitor, progeny)
		#append to list
		pairs_values.append(diff)
	#return
	return pairs_values

def form_internode(pair, matrix_dict, reference):
	"""
	Create an internode entry
	"""
	pair_dict = ""
	#define samples
	sample1 = pair.split("|")[0]
	sample2 = pair.split("|")[1]
	#intersect dicts
	for i in xrange(len(matrix_dict[sample1])):
		#add intersections
		if matrix_dict[sample2][i] == matrix_dict[sample1][i]: pair_dict += matrix_dict[sample1][i]
		else: pair_dict += reference[i]
	#return
	return pair_dict

def compare_entries(matrix_dict1,  matrix_dict2):
	#see how many letters defer between two entries
	diff = 0
	for i in xrange(len(matrix_dict1)):
		if matrix_dict2[i] != matrix_dict1[i]: diff += 1
	#return
	return diff

if __name__ == "__main__":
		exit(main())
