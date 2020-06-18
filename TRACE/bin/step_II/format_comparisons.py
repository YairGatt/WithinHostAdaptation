#!/usr/bin/env python

"""
Format comparisons to just get all the comparisons we want to do. Callable both from command line and python
"""

import sys, os
import argparse
import itertools

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
		'-c', '--comparisons_file', default="./comparisons_converted.txt",
		help='File with comparisons in full format.')
	
	parser.add_argument(# customized description; put --help last
		'-h', '--help', action='help',
		help='Show this help message and exit.')
	
	settings = parser.parse_args(argv)
	
	return settings


def absolute_comparisons(comparisons):
	"""
	Format all the comparisons changing them from full format to each node having alll subnodes
	"""
	#read lines
	with open(comparisons,"r") as fl:
		lines = fl.readlines()
	#initialize dict
	children_dict = {}
	#split to comparisons
	base_comparisons = [i.strip().split() for i in lines]
	#create basic children dict
	for comparison in base_comparisons:
		if "|" in comparison[0]:
			try:
				children_dict[comparison[0]].append(comparison[1])
			except KeyError:
				children_dict[comparison[0]] = [comparison[1]]
	#add children of children for each node until reaching end
	while any(list(itertools.chain.from_iterable([["|" in j for j in i] for i in children_dict.values()]))):
		for node in children_dict:
			for child in children_dict[node]:
				if child in children_dict:
					children_dict[node].remove(child)
					children_dict[node] += children_dict[child]
	#create absoulte comparisons
	absolute_comparisons = ["|".join(i) for i in children_dict.values()]
	#return
	return absolute_comparisons

def format_comparisons(comparisons):
	"""
	Format all the comparisons changing them from full format to actual comparisons
	"""
	#read lines
	with open(comparisons,"r") as fl:
		lines = fl.readlines()
	#initialize dict
	replacement_dict = {}
	#split to comparisons
	base_comparisons = [i.strip().split() for i in lines]
	#create basic replacement dict 
	for comparison in base_comparisons:
		if "|" in comparison[1]:
			replacement_dict[comparison[1]] = comparison[0]
	#replace until condition is satisfied
	while any(list(itertools.chain.from_iterable([["|" in j for j in i] for i in base_comparisons]))):
		for comparison in base_comparisons:
			if "|" in comparison[0]:
				try:
					comparison[0] = replacement_dict[comparison[0]]
				except KeyError: #nothing leading up to internode
					comparison[0] = "None"
			
			if "|" in comparison[1]:
				try:
					comparison[1] = replacement_dict[comparison[1]]
				except KeyError: #nothing leading up to internode
					comparison[1] = "None"
	
	#remove 0s
	base_comparisons = [i for i in base_comparisons if "None" not in i and "fiction" not in i and i[0] != i[1]]
	#return
	#print base_comparisons
	return base_comparisons

def write_to_file(outfile, formatted_comparisons):
	"""
	Write results to file
	"""
	with open(outfile,"w") as outfl:
		for comparison in formatted_comparisons:
			line = "%s\t%s\n" % (comparison[0],comparison[1])
			outfl.write(line)

def main(argv=None):
	#process command line
	settings = process_command_line(argv)
	#define outfile
	current_dir = os.path.dirname(settings.comparisons_file)
	if current_dir: outfile = current_dir + "/" + "formatted_comparisons_converted.txt"
	else: outfile = "formatted_comparisons_converted.txt"
	#format comparisons
	formatted_comparisons = format_comparisons(settings.comparisons_file)
	#write to file
	write_to_file(outfile, formatted_comparisons)

if __name__ == "__main__":
		exit(main())
