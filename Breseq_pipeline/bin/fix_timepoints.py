import os,sys

#python /home/usesr/yair/PhD_projects/project_B/bin/gene_loss_pipeline/breseq_pipeline/fix_timepoints.py "${workdir}/temp_files/diff.gd" $timepoint1 $timepoint2
"""
fix timepiotn format to match diff file 
"""

diff_file = sys.argv[1]
timepoint1 = sys.argv[2]
timepoint2 = sys.argv[3]

with open(diff_file) as fl:
    line = fl.readline()
    line = fl.readline()
    
split_line = line.split()
frequencies = [i.split("=")[0].split("frequency_")[1] for i in split_line if "frequency_" in i]

for i in frequencies:
    if i in timepoint1:
        time1 = i
    if i in timepoint2:
        time2 = i

print time1,time2

