#!/usr/bin/env python

"""
Utilities for combinatorics
"""

import sys, os

def sublisting(alist):
	"""
	Getting all sublists of alist
	"""
	#define number of items
	n = len(alist)
	#get all divisions to subgroups as numerical vectors
	final_list = finalizing(n)
	#convert to items
	final_list_converted = []
	#iterate through all subdivisions
	for sublist in final_list:
		#initialize converted list
		converted_list = [[] for i in xrange(max(sublist)+1)]
		#iterate through subdivision
		for n,i in enumerate(sublist):
			#add appropriate item
			converted_list[i].append(alist[n])
		#add converted sublist to list of converted sublists
		final_list_converted.append(sorted(converted_list))
	#remove lists that repeat multiple times
	finalest_list = []
	for i in final_list_converted:
		if i not in finalest_list: finalest_list.append(i)
	#return
	return finalest_list

def finalizing(n):
	"""
	Get all subdivisions of all lengths in numerical vectors for the number of items in n
	"""
	final_list = []
	#iterate through all possible numbers of subgroups
	for number in xrange(n):
		number_of_subgroups = number + 1
		count = 0
		char_list = ""
		#define char list
		for i in xrange(number_of_subgroups):
			char_list += str(count)
			count+= 1
		#get all subdivisions for this number of subgroups
		final_list += print_sequences(char_list,n)
	#return
	return final_list

def print_sequences(char_list, n):
	# gets all possible combinations of characters in char_list in length n
	#initialize list of sublists
	list_of_lists = []
	#add all sublists
	for seq in sub_sequences(char_list, n):
		list_of_lists.append([int(i) for i in seq])
	#remove lists with incorrect number of groups
	final_list_of_lists = []
	for list in list_of_lists:
		#flag means the list is not already in the final list
		flag = 0
		for char in char_list:
			#if it does not include one of the characters, which signify subgroups, it does not have the correct number of groups
			if int(char) not in list: flag = 1
		#append if still relevant
		if flag == 0: final_list_of_lists.append(list)
	#return
	return final_list_of_lists

def sub_sequences(char_list, n):
	# supports print_sequences by creating the different combinations
	if n == 0: return []
	if n == 1: return char_list
	else:
		#recursively subsequences
		sub_lst = sub_sequences(char_list, n - 1)
		#return
		return [arg + char for arg in sub_lst for char in char_list]

if __name__ == "__main__":
		exit(main())