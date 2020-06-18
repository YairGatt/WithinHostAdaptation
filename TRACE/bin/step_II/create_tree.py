"""
This script creates a core SNP matrix using kSNP for TRACE II in case the matrix created in TRACE I was deleted (space saving, accidental, etc.)
"""

import os,sys
import errno
from subprocess import call
import shutil
import re
from os import listdir
from os.path import isfile, join

def mkdir(path):
	#create directory and don't crash if it already exists
	try: os.mkdir(path)
	except OSError as exc:
		if exc.errno != errno.EEXIST: raise exc
		pass

def separate(string,to_return):
	#separate number and letters in string
	match = re.match(r"([0-9]+)([a-z]+)", string, re.I)
	if match:
		if to_return == "number": return int(match.groups()[0])
		if to_return == "letters": return match.groups()[1]
	else: return False

#input directory full path
input_directory = sys.argv[1]
#input repository with kSNP format of organism strains from NCBI
repository = sys.argv[2]
#temp dir
temp_dir = sys.argv[3]
#converting file between strains and isolate names
convert_file = sys.argv[4]
#outdir
workdir = sys.argv[5]
#kSNP installation
kSNP_dir = sys.argv[6]
#strains mode, create tree only for isolate in strains
strains = sys.argv[7]
#create outdir
mkdir(workdir)
if strains:
	with open(strains) as fl: current_strains = [i.strip() for i in fl.readlines()]
#create directories
mkdir(workdir)
#initialize kSNP input file
input_file = workdir + "/kSNP_input.txt"
#initialize time points
times = []
#write input file
with open(input_file,"w") as input:
	#add timepoints
	with open(convert_file,"r") as convert:
		#get lines
		convert_lines = convert.readlines()[1:]
		#iterate through lines
		for line in convert_lines:
			#get time point
			split_line = line.split()
			the_time = split_line[-1]
			#convert to int or remove letter signifying scale (days, weeks, etc.)
			try: int(the_time)
			except ValueError: the_time = separate(the_time,"number")
			#add the number of appearances of the time point so far as a cahracter to distinguish isolates from the same time point
			value = str(the_time)+chr(ord('a')+times.count(the_time))
			#add to time points
			times.append(the_time)
			#get name of isolate
			name = split_line[1].replace("_1.fastq.gz","")
			#get assembly directory in input_directory
			sample_directory = input_directory + "/" + name + "/"
			#if there is no directory something is wrong
			if not os.path.isdir(sample_directory): continue
			#get files
			onlyfiles = [f for f in listdir(sample_directory) if isfile(join(sample_directory,f))]
			#check for existance of .fna file, if it does assembly_type = True else = False
			assembly_type = False
			for f in onlyfiles:
				if f.endswith(".fna"):
					assembly_type = True
					assembly = f
			#go
			if not assembly_type: #the name is scaffolds.fasta
				#see if relevant clone
				if strains:
					if value not in current_strains: continue
				#get file
				file_name = sample_directory + "scaffolds.fasta"
				#in case the size of fastq was not consistent between pair mates
				matched_file_name = file_name.replace("/scaffolds.fasta","_matched/scaffolds.fasta")
				#write to kSNP input file
				if os.path.isfile(file_name):
					input.write(file_name + "\t" + value + "\n")
				elif os.path.isfile(matched_file_name):
					input.write(matched_file_name + "\t" + value + "\n")
				else:
					raise Exception("Assembly does not exist %s" % name)
			else: #the name is something.fna
				#copy to temp to fix name if it has dots and such which kSNP can't handle
				old_name = sample_directory + "/" + assembly
				new_name = temp_dir + "/" + name.replace(".","-") + "_genomic.fna"
				#copt file
				with open(old_name) as infile:
					with open(new_name,"w") as outfile:
						print new_name
						for line in infile:
							outfile.write(line)
				#check if clone
				if strains:
					if value not in current_strains: continue
				#write to kSNP input
				input.write(new_name + "\t" + value + "\n")
	#this parameter decides how many of the repository organisms will be used, -1 means all. Here we just want the snp matrix between our guys so no need to add strains from NCBI if not neccesary
	balance = max(4 - len(current_strains), 1)
	#add a few strains from ncbi for balance
	with open(repository,"r") as rep:
		rep_lines = rep.readlines()
		for line in rep_lines:
				if balance != 0:
					input.write(line)
					balance -= 1
#run kSNP
#call(kSNP_dir + "/kSNP3 -in " + input_file + " -outdir " + workdir + " -k 13 -ML -NJ", shell=True)
call(kSNP_dir + "/kSNP3 -in " + input_file + " -outdir " + workdir + " -k 13 -core", shell=True)
#print kSNP_dir + "/kSNP3 -in " + input_file + " -outdir " + workdir + " -k 13 -core"
#delete temp files
try:
   shutil.rmtree("/home/users/yair/temp")
   pass
except OSError:
   pass
