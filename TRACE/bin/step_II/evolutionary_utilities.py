#!/usr/bin/env python

"""
Utility functions for create_tree_evolutionary_model.py
"""

import sys, os
import errno
import argparse
import re
import itertools
import numpy as np
import matplotlib.pyplot as plt

class Sample(object):
	"""
	A sample isolated from a patient
	"""
	def __init__(self, name, time):
		self.name = name
		self.time = separate(time, "number")
		self.options = []
		
	def __str__(self):
		return self.name

def mkdir(path):
	#create directory and don't crash if it already exists
	try:
		os.mkdir(path)
	except OSError as exc:
		if exc.errno != errno.EEXIST:
			raise exc
			pass

def separate(string,to_return):
    #separate number and letters in string
    match = re.match(r"([0-9]+)([a-z]+)", string, re.I)
    #print string
    #print match
    if match:
        if to_return == "number":
            #print string
            return int(match.groups()[0])
        if to_return == "letters":
            return match.groups()[1]
    else:
        return None

def unpack_strains(strain):
	#this function reads list of samples in a file to a list
	#open file
	with open(strain) as fl: samples = [i.strip() for i in fl.readlines()]
	#return
	return samples

def convert_strain(id, strain, converted):
	"""
	Convert time id to sample name
	"""
	#define samples
	samples = unpack_strains(converted)
	samples_times = unpack_strains(strain)
	#check matching index
	converted_id = samples[samples_times.index(id)]
	#return
	return converted_id

def strain_converted(id, strain, converted):
	"""
	Convert sample name to time id
	"""
	#define samples
	samples = unpack_strains(converted)
	samples_times = unpack_strains(strain)
	#check matching index
	converted_id = samples_times[samples.index(id)]
	#return
	return converted_id

def samples_to_objects(strain, converted):
	"""
	Create sample objects from the samples in the experiment
	"""
	#initialize
	Samples = []
	#define samples
	samples = unpack_strains(converted)
	samples_times = unpack_strains(strain)
	#match and define objects
	for n, accession in enumerate(samples):
		Samples.append(Sample(name=accession, time=samples_times[n]))
	#return
	return Samples

def convert_to_dict(pairs,pairs_values):
	"""
	Convert the two lists of the pairs and their vaules to proper dictionary form for easy retrieval of data
	"""
	#initialize dict
	pairs_dictionary = {}
	reverse_pairs_dictionary = {}
	#iterate through all pairs
	for n, pair in enumerate(pairs):
		#define progenitor with same pattern always
		progenitor = pair[0]
		if "|" in progenitor: progenitor = "%s|%s" % (min(progenitor.split("|")),max(progenitor.split("|")))
		#define progeny with same patter always
		progeny = pair[1]
		if "|" in progeny: progeny = "%s|%s" % (min(progeny.split("|")),max(progeny.split("|")))
		#add to dict
		if progenitor in pairs_dictionary:
			if progeny in pairs_dictionary[progenitor]:
				if pairs_dictionary[progenitor][progeny] != pairs_values[n]: raise Exception("Identical pairs %s %s" % (progenitor, progeny))
			else: pairs_dictionary[progenitor][progeny] = pairs_values[n]
		else:
			pairs_dictionary[progenitor] = {}
			pairs_dictionary[progenitor][progeny] = pairs_values[n]
		#add to reverse dict
		if progeny in reverse_pairs_dictionary:
			if progenitor in reverse_pairs_dictionary[progeny]:
				if reverse_pairs_dictionary[progeny][progenitor] != pairs_values[n]: raise Exception("Identical pairs with differing values %s %s" % (progenitor, progeny))
			else: reverse_pairs_dictionary[progeny][progenitor] = pairs_values[n]
		else:
			reverse_pairs_dictionary[progeny] = {}
			reverse_pairs_dictionary[progeny][progenitor] = pairs_values[n]
	#return
	return pairs_dictionary, reverse_pairs_dictionary

def plot_histogram(ranked_trees):
	n, bins, patches = plt.hist(x=ranked_trees, bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
	plt.grid(axis='y', alpha=0.75)
	plt.xlabel('Value')
	plt.ylabel('Frequency')
	plt.title('Number of evolutionary events in trees')
	maxfreq = n.max()
	# Set a clean upper y-axis limit.
	plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
	plt.show()

if __name__ == "__main__":
		exit()
