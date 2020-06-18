import sys,os
import csv

def get_lost(direction):
    #get direction
    direction_name = direction.split("_")[-2]
    #if no comparison remove direction
    if direction_name == "snpEff": direction_name = ""
    else: direction_name = "_" + direction_name
    name = "snpEff%s_summary.genes.txt" % direction_name
    full_name = "snpEff%s_summary.genes-all.txt" % direction_name
    base = direction.replace(name,"").replace(full_name,"")
    output = base + "lost_genes%s.txt" % direction_name
    output_high = base + "high_lost_genes%s.txt" % direction_name
    #parse_file
    lost = []
    lost_high = []
    not_lost = []
    not_lost_high = []
    with open(direction) as fl:
        content = fl.readlines()[1:]
        reader = csv.DictReader(content,delimiter="\t")
        for row in reader:
            gene = row['TranscriptId']
            high = 0
            moderate = 0
            if 'variants_impact_HIGH' in row:
                high = int(row['variants_impact_HIGH'])
                if high > 0:
                    lost.append(gene)
                    lost_high.append(gene)
                    continue
            not_lost_high.append(gene)
            if 'variants_impact_MODERATE' in row:
                moderate = int(row['variants_impact_MODERATE'])
                if moderate > 0:
                    lost.append(gene)
                    continue
            not_lost.append(gene)
    #sort_list
    lost.sort()
    not_lost.sort()
    lost_high.sort()
    not_lost_high.sort()
    #write_output
    with open(output,"w") as outfl:
        for i in lost: outfl.write(i + "\n")
    with open(output_high,"w") as outfl:
        for i in lost_high: outfl.write(i + "\n")
    #if no direction, write non lost as well
    if not direction_name:
        #get output name
        non_output = base + "non_lost_genes.txt"
        non_output_high = base + "high_non_lost_genes.txt"
        #write
        with open(non_output,"w") as outfl:
            for i in not_lost: outfl.write(i + "\n")
        with open(non_output_high,"w") as outfl:
            for i in not_lost_high: outfl.write(i + "\n")

#input
#straight
straight = sys.argv[1]
get_lost(straight)
#reverse
if len(sys.argv) > 2:
    reverse = sys.argv[2]
    get_lost(reverse)
