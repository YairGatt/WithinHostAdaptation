#!/usr/bin/env python

"""
Single tree utilities
"""

import sys, os
from evolutionary_utilities import Sample 

class Tree(object):
	"""
	An optional tree made of pairs (nodes)
	"""
	def __init__(self):
		self.nodes = []
		self.rank = None
		self.dict = {}
		
	def included(self,a,b):
		for i in self.nodes:
			if i.progenitor == a and i.progeny == b:
				return True
		return False
	
	def to_dict(self):
		#reset self.dict
		self.dict = {}
		#define self.dict from nodes
		for i in self.nodes:
			try: self.dict[i.progeny.name] = i.progenitor.name
			except AttributeError: self.dict[i.progeny.name] = i.progenitor
	
	def draw_progenitor(self,i):
		try: return self.dict[i]
		except KeyError: return None
			
	def loops(self):
		"""
		Check for loops
		"""
		for i in self.dict:
			a = i
			i_list = [i]
			#keep going back recursively and see if you end up with anythin already on the path
			while a != None:
				#step back
				a = self.draw_progenitor(a)
				#check for loop
				if a in i_list: return True
				#add to nodes already passed
				i_list.append(a)
		#return no loops found
		return False
	
	def double_progenitor(self):
		"""
		Check if there is a progeny with two progenitors for some reason
		"""
		#reset self.dict
		self.dict = {}
		#define self.dict from nodes
		for i in self.nodes:
			#get progeny
			if type(i.progeny) == str: progeny = i.progeny
			else: progeny = i.progeny.name
			#get progenitor
			if type(i.progenitor) == str: progenitor = i.progenitor
			else: progenitor = i.progenitor.name
			#check if progeny already has another progenitor
			if progeny in self.dict and progenitor != self.dict[progeny]: return True
			#add to dict
			self.dict[progeny] = progenitor
		#no double progenitor foudn
		return False
	
	def get_root(self):
		self.root = []
		for i in self.nodes:
			flag = 0
			for n in self.nodes:
				if n.progeny == i.progenitor and n.progenitor != "fiction": flag = 1
			if flag == 0: self.root.append(i.progenitor)
		
		return list(set(self.root))
	
	def get_root_simple(self):
		self.root = []
		for i in self.nodes:
			if i.progenitor == "fiction": self.root.append(i.progeny)
		
		return list(set(self.root))
	
	def sort_nodes(self):
		for n,i in enumerate(self.nodes):
			flag = 0
			if type(i.progeny) == str:
				progeny = "%s|%s" % (min(i.progeny.split("|")),max(i.progeny.split("|")))
				flag = 1
			else:
				progeny = i.progeny
			if type(i.progenitor) == str and i.progenitor != "fiction":
				progenitor = "%s|%s" % (min(i.progenitor.split("|")),max(i.progenitor.split("|")))
				flag = 1
			else:
				progenitor = i.progenitor
			#add regular nodes
			if flag: self.nodes[n] = Node(progenitor,progeny)
			
class Node(object):
	"""
	A node in a potential phylogenetic tree
	"""
	def __init__(self, progenitor,progeny):
		self.progenitor = progenitor
		self.progeny = progeny
							
def create_tree(raw_tree, Samples, index):
	"""
	Create a single tree from raw_tree
	"""
	#initialize index of sample
	count = 0
	if count == index: count += 1 #index to be skipped
	#initialize final tree
	final_tree = Tree()
	#add each sample to final tree in proper format
	for origin in raw_tree:
		#add node
		final_tree.nodes.append(Node(origin, Samples[count]))
		#add to index
		count += 1
		if count == index: count += 1 #index to be skipped
	#don't append tree if has loops
	final_tree.to_dict()
	if final_tree.loops(): return None
	#if pairs of samples from same time point exist, change the format to include and internode
	final_tree = get_internodes(final_tree)
        if final_tree.double_progenitor(): return None
	#sort nodes
	final_tree.sort_nodes()
	#return
	return final_tree

def get_internodes(final_tree):
	"""
	See if there are nodes with samples from the same time point and create an internodes to replace them
	"""
	#initialize nodes
	iterable = list(final_tree.nodes)
	#dict internodes
	internode_dict = {}
	#round one, simply add internodes
	for n, pair in enumerate(iterable):
		try:
			if pair.progenitor.time == pair.progeny.time:
				#replace node with special annotation
				internode = "%s|%s" % (pair.progenitor.name, pair.progeny.name)
				#add to dict
				try:
					internode_dict[pair.progenitor.name].append(internode)
				except KeyError:
					internode_dict[pair.progenitor.name] = [internode]
				#Create proxy nodes
				new_link1 = Node(internode, pair.progenitor)
				new_link2 = Node(internode, pair.progeny)
				#Add nodes
				final_tree.nodes.append(new_link1)
				final_tree.nodes.append(new_link2)
				final_tree.nodes[n] = None
		except AttributeError: pass
	#remove Nones
	while None in final_tree.nodes: final_tree.nodes.remove(None)
	#round two, look for nodes that lead up to the internodes (their progeny is the progenitor strain of the internode), and replace accordingly
        for n, node in enumerate(final_tree.nodes):
            if node.progeny.name in internode_dict: final_tree.nodes[n] = Node(node.progenitor, Sample(internode_dict[node.progeny.name], str(node.progeny.time)+"a"))
	#return
	return final_tree
	
def get_trees_internodes(final_tree, internode_dict):
	#replace non-internode nodes with internode progeny and remove from dict
	for n, pair in enumerate(final_tree.nodes):
		if type(pair.progenitor) != str: #non_internode
			if pair.progeny.name in internode_dict:
				#choose random internode as progeny
				internode = random.choice(internode_dict[pair.progeny.name])
				#remove from dict
				internode_dict[pair.progeny.name].remove(internode)
				#replace
				final_tree.nodes[n] = Node(pair.progenitor, internode)
	#replace internode nodes with internode progeny and remove from dict
	for n, pair in enumerate(final_tree.nodes):
		if type(pair.progenitor) == str: #internode
			if pair.progeny.name in internode_dict:
				#choose random internode as progeny
				non_self_hit = [i for i in internode_dict[pair.progeny.name] if i != pair.progenitor]
				if not non_self_hit: continue
				internode = random.choice(non_self_hit)
				#remove from dict
				internode_dict[pair.progeny.name].remove(internode)
				#replace
				final_tree.nodes[n] = Node(pair.progenitor, internode)
				
	return final_tree



if __name__ == "__main__":
		exit(main())
