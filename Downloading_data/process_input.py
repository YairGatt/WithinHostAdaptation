"""
This script accepts tab-delimited-CSV input from Pathoadaptation DB and outputs the fastq files for millstone, along with the target file to be uploaded
"""

import sys,os
import determine_best_genome, fastq_determine_best_genome
import csv
import errno
from Bio import Entrez
import argparse
import glob
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
	
	# define options here:
	parser.add_argument(
		'-D' , '--database', default='/home/users/yair/Documents/PhD_projects/project_B/data/final_database.csv',
		help='Database.')
	
	parser.add_argument(
		'-o', '--organism', default="Pseudomonas_aeruginosa",
		help='Organism to be analyzed.')

        parser.add_argument(
               '--EDIRECT', default="/home/users/yair/Software/edirect",
               help='edirect tools installation.')

        parser.add_argument(
               '--QUAST', default="/home/users/yair/Software/quast-4.5",
               help='Quast installation.')

        parser.add_argument(
               '--mail', default="yair.gatt@mail.huji.ac.il",
               help='Email for edirect.')

        parser.add_argument(
               '--sratools', default="/home/users/yair/Software/sratoolkit/bin",
               help='SRAtools installation.')
	
	parser.add_argument(
		'-t' , '--type', default='SRA',
		help='SRA/Run or Assembly or Local.')

        parser.add_argument(
                '-i' , '--installation', default='/home/users/yair/git/WithinHostAdaptation',
                help='Pipeline installation.')

	parser.add_argument(
		'--skip', action="store_true", default=False,
		help='Skip lines in db corresponding to folders already created.')
	
	parser.add_argument(
		'-w', '--workdir', default="./temp.txt",
		help='Workdir where results will be written.')
	
	parser.add_argument(# customized description; put --help last
		'-h', '--help', action='help',
		help='Show this help message and exit.')
	
	settings = parser.parse_args(argv)
	
	return settings

def process_csv(database, organism, type):
	"""
	This function reads the pathoadaptation DB file and extracts lines with the relevant organism and accession type. It then creates a list of tuples for each line, with (SRA,timepoint) and returns
	a list of tuple of the index of each line and the line list.
	"""
	#initialize sample list and index dict
	samples = []
	experiment_dict = {}
	organism = organism.replace("_"," ")
	#parse csv
	with open(database, "rb") as fl:
		reader = csv.DictReader(fl, delimiter="\t")
		for row in reader:
                        if organism in row["Bacteria"]:
        			#initialize sample list of tuples
				sample = []
				#extract accessions sorted
				accessions = sorted([(int(i.split("_")[0].split("Timepoint")[1]),row[i]) for i in row.keys() if "_accession" in i and row[i] != ""])
				#extract timepoints sorted
				timepoints = sorted([(int(i.split("_")[0].split("Timepoint")[1]),row[i]) for i in row.keys() if "_time" in i and row[i] != ""])
				#define index as experiment_patient
				experiment = row["Experiment"]
				try: experiment_dict[experiment] += 1
				except KeyError: experiment_dict[experiment] = 1
				index = experiment + "_" + str(experiment_dict[experiment])
				#go through the accessions and create tuples
				for i in accessions:
					accession = i[1]
					#find relevant timepoint
					position = i[0]
					for j in timepoints:
						if j[0] == position: timepoint = j[1]
                                        #list of relevant accessions in case the cell includes several accessions
	                        	if type.lower() == "assembly" and "assembly" in row["Accession_type"].lower():
						accession_list = [n for n in accession.split("/") if n.startswith("GC")]
                                        elif (type.lower() == "sra" or type.lower() == "run") and ("sra" in row["Accession_type"].lower() or "run" in row["Accession_type"].lower()):
				                accession_list = [i for i in accession.split("/") if not i.startswith("GC")]
		                        elif type.lower() == "local" and "local" in row["Accession_type"].lower():
						accession_list = [accession]
					else: continue

					#add accession to sample
					sample.append((accession_list, timepoint))
				#add to samples	
				samples.append((index, sample))
	#return
	return samples

def mkdir(path, skip):
	#create directory and don't crash if it already exists
	try:
		os.mkdir(path) 
		return path
	except OSError as exc:
		if not skip: return path
		if exc.errno != errno.EEXIST: raise exc

def create_accessions_files(workdir, samples, skip):
	"""
	This function accepts the samples list of sample lists from the process_csv function, and creates an accession file for each sample list, it then returns a list with all files. They can later be downloaded
	"""
	#initialize dict
	accessions_files_list = []
        #iterate through samples
	for line in samples:
                #create patient dir
		index = line[0]
		split_index = index.split("_")
		path = workdir + "/trial_%s_patient_%s" % (split_index[0], split_index[1])
		#create directory
		path = mkdir(path, skip)
		if path == None: continue
		#gather accessions
		accessions = []
		timepoints = []
		for sample in line[1]:
			accessions += sample[0]
			timepoints += [sample[1]]*len(sample[0]) #remember there may be multiple accessions for each sample	
		#write file
		if len(accessions) != len(timepoints): raise Exception("Something is wrong, accessions do not match timepoints")
		with open(path + "/accessions_file.txt", "w") as accessions_file:
			for accession in accessions: accessions_file.write(accession + "\n")
		with open(path + "/timepoints_file.txt", "w") as timepoints_file:
			for timepoint in timepoints: timepoints_file.write(timepoint + "\n")
		#add to list
		accessions_files_list.append(path + "/accessions_file.txt")
	#return
	return accessions_files_list

def write_convert_files(workdir, organism, accession_type="sra"):
	"""
	This function write every relevant line in samples to a convert file format so they can be used for the pipeline
	"""
	#extract relevant samples
	with open(workdir+"/accessions_file.txt","r") as accessions_file:
		accessions = accessions_file.readlines()
		accessions = [i.strip() for i in accessions]
	with open(workdir+"/timepoints_file.txt","r") as timepoints_file:
		timepoints = timepoints_file.readlines()
		timepoints = [i.strip() for i in timepoints]
	#split samples to timepoints dict
	timepoint_dict = {}
	for index, timepoint in enumerate(timepoints):
		if accession_type.lower() == "run" or accession_type.lower() == "sra": final_accessions = [accessions[index]]
		elif accession_type.lower() == "assembly":
			accessions_dir = os.listdir(workdir + "/assemblies/")
			accession = accessions[index]
                        #get final
			final_accessions = []
			for i in accessions_dir:
				if accession in i: final_accessions.append(i)
		elif accession_type.lower() == "local": final_accessions = [accessions[index]]
		else: raise Exception("type must be SRA/Run/Assembly")
		#add to dict
		try: timepoint_dict[timepoint] += final_accessions
		except KeyError: timepoint_dict[timepoint] = final_accessions
	#write to file by groups
	with open(workdir + "/convert_file.txt","w") as convert_file:
		#header
		count_list = []
		for a_timepoint in timepoint_dict.values():
			for a_sample in a_timepoint: count_list.append(len(glob.glob(workdir + "/fastqs/" + a_sample + "*")))
		try: example_count = max(count_list)
		except ValueError: example_count = 0
                #how many files are there		
		if example_count == 1 or accession_type.lower() == "assembly":
			convert_file.write("Sample_Name\tRead_1_Filename\tTimepoint")
		elif example_count == 2:
			convert_file.write("Sample_Name\tRead_1_Filename\tRead_2_Filename\tTimepoint")
		elif example_count == 0:
			convert_file.write("Sample_Name\tRead_1_Filename\tTimepoint")
		elif example_count > 2:
			raise Exception("File count for accession %s must be 1 or 2, investigate the disperency in %s" % (timepoint_dict.values()[0][0], workdir+"/fastqs/"))
		#iterate timepoints
		for timepoint in timepoint_dict.keys():
			#iterate through accessions
			for n, accession in enumerate(timepoint_dict[timepoint]):
				#see how many file there are
				#list files in workdir/fastqs containing accession
				count = len(glob.glob(workdir + "/fastqs/" + accession + "*"))
				if workdir[-1] == "/": truncated_workdir = workdir[:-1]
				else: truncated_workdir = workdir
				sample_name = os.path.basename(truncated_workdir) + "_" + "_" + timepoint + "_" + str(n+1)
				if accession_type.lower() == "run" or accession_type.lower() == "sra" or accession_type.lower() == "local":
					file1 = accession+"_1.fastq.gz"
					file2 = accession+"_2.fastq.gz"
				elif accession_type.lower() == "assembly": file1 = accession
				else: raise Exception("type must be SRA/Run/Assembly")
				#if 0, the file probably failed QUAST, user should confirm
				if count == 0 and (accession_type.lower() == "run" or accession_type.lower() == "sra" or accession_type.lower() == "local"):
					warnings.warn("No fastq files for accession %s, likely failed QUAST filtering. make sure it downloaded properly" % accession)
				#if 1, add to convert file as single end
				elif count == 1 or accession_type.lower() == "assembly":
					if example_count == 1:
						convert_file.write("\n" + sample_name + "\t" + file1 + "\t" + timepoint)
					elif example_count == 2:
						convert_file.write("\n" + sample_name + "\t" + file1 + "\t" + "\t" + timepoint)
					else:
						convert_file.write("\n" + sample_name + "\t" + file1 + "\t" + timepoint)
				#if 2, add to millstone as paired end
				elif count == 2:
					if example_count == 2:
						convert_file.write("\n" + sample_name + "\t" + file1 + "\t" + file2 + "\t" + timepoint)
					elif example_count == 1:
						convert_file.write("\n" + sample_name + "\t" + file1 + "\t" + timepoint)
				#if more, odd but unlikely problem
				elif count > 2:
					raise Exception("File count for accession %s must be 1 or 2, investigate the disperency in %s" % (accession, workdir+"/fastqs/"))
				else:
					raise Exception("File count is negative or not an interger")
		
def main(argv=None):
	#Process command line
	settings = process_command_line(argv)
	#Process CSV to output a list of the lines, where each element is a list of tuples (experiment, SRA,timepoint)
	samples = process_csv(settings.database, settings.organism, settings.type)
	#create directories - create workdir, inside it a directory for each line/patient, inside it a directory for each timepoint
	path = mkdir(settings.workdir, settings.skip)	
	#form accessions list with the accession files
	accessions_files_list = create_accessions_files(settings.workdir, samples, settings.skip)
	#run (fastq_)determine_best_genome
	if settings.type.lower() == "run" or settings.type.lower() == "sra":
		for accessions_file in accessions_files_list:
			arguements = argparse.Namespace(installation=settings.installation, workdir=os.path.dirname(accessions_file) ,organism=settings.organism,sras_file=accessions_file,representative=settings.representative)
			fastq_determine_best_genome.main(arguements)
	elif settings.type.lower() == "assembly":
		for accessions_file in accessions_files_list:
			arguements = argparse.Namespace(installation=settings.installation, QUAST=settings.QUAST, EDIRECT=settings.EDIRECT, mail=settings.mail, workdir=os.path.dirname(accessions_file),organism=settings.organism,assembly_file=accessions_file)
			determine_best_genome.main(arguements)
        #create convert file for each patient
	write_convert_files(os.path.dirname(accessions_file), settings.organism, settings.type)	
	
if __name__ == "__main__":
	exit(main())
