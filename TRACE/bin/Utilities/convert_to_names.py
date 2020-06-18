"""
This script converts the results of divide_to_strains.py to proper isoalte names
"""

import os,sys
import re

def separate(string,to_return):
    #separate number and letters in string
    match = re.match(r"([0-9]+)([a-z]+)", string, re.I)
    #print string
    if match:
        if to_return == "number":
            #print string
            return int(match.groups()[0])
        if to_return == "letters":
            return match.groups()[1]
    else:
        return int(string)

#input
strain = sys.argv[1]
convert_file = sys.argv[2]
#converting dict
convert_dict = {}
#open convert file
with open(convert_file,"r") as convert:
    #get lines
    lines = convert.readlines()[1:]
    #initialize list of time points
    times = []
    #iterate through lines
    for line in lines:
        #get time point
        split_line = line.split()
        the_time = split_line[-1]
        #add number of appearances of time point so far as character to diffrentiate isolates from the same time point
        value = str(separate(the_time,"number")) + chr(ord('a')+times.count(the_time))
        #if too much add uppercase
        if ord('a') + times.count(the_time) > 122: value = str(separate(the_time,"number")) + chr(ord('a') + times.count(the_time) - 58)
        #add to list of time points
        times.append(the_time)
        #get proper name of isolate
        name = split_line[1].replace("_1.fastq.gz","")
        #add to converting dict
        convert_dict[value] = name
#read strain infromation
with open(strain,"r") as fl:
    #write to converted file
    with open(strain.split("strain_")[0] + "converted_" + strain.split("strain_")[1],"w") as outfl:
        #iterate through lines
        for line in fl:
            #convert
            new_line = convert_dict[line.strip()]
            #write
            outfl.write(new_line+"\n")
