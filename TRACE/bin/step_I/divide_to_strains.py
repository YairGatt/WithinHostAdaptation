"""
This tree takes a maximum likelihood tree with isolates from a strain and other NCBI complete genomes and divides the strains isolates to clones
"""
import sys,os
from Bio import Phylo
import argparse

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
        '-t', '--tree', nargs='+',
        help='List of trees to be processed.')

    parser.add_argument(
        '-f', '--format', default=["newick"], nargs='+',
        help='Format of tree files, can be given as one format or a list of length matching tree_files for each tree.')

    parser.add_argument(
        '-w', '--workdir', default="./temp/",
        help='Workdir where results will be written.')

    parser.add_argument(# customized description; put --help last
        '-h', '--help', action='help',
        help='Show this help message and exit.')

    settings = parser.parse_args(argv)

    return settings

def strains_from_tree(tree_file, tree_format):
    """
    Parse tree to get strains
    """
    #build tree
    tree = Phylo.read(tree_file,tree_format)
    #root tree to reference
    tree.root_with_outgroup("reference")
    #initialize strains dict
    strains_dict = {}
    #define terminals as initial group whose parents need to be annotated
    terminals = tree.get_terminals()
    group = terminals
    #run recursively through tree until reaching root - also works for unrooted trees
    while len(group) > 0:
        #prepare next group
        new_group = []
        #iterate through current group and find parents to annotate
        for terminal in group:
             #check if parent not root, otherwise skip
             try:
                step_back = tree.get_path(terminal)[-2]
             except IndexError:
                 continue
            #continue if the parent is annotated
             if step_back.name: continue
             #define all children
             children = step_back.clades
             cut_flag = 0
             #check if all children are annotated, otherwise don't annotate yet
             if all([i.name for i in children]):
                #find all the names of the children
                names = [i.name for i in children]
                #if "cut" in children, cut!
                if "cut" in names: cut_flag = 1
                #count GCFs (NCBI strains)
                strains = []
                GCFs = []
                #iterate through children annotations and see how many NCBI strains are there
                for i in names:
                    if i.startswith("GCF") or i == "reference": GCFs.append(i)
                    else: strains += i.split("#")
                #determine step_up name
                if len(GCFs) >= 1:
                    step_back.name = "#".join(GCFs)
                    #add strains_dict if needed 
                    if len(strains) >= 1:
                        if step_back.name in strains_dict: strains_dict[step_back.name] += strains
                        else: strains_dict[step_back.name] = strains
                #if no GCFs, remember the strains in the node name and continue
                elif len(GCFs) ==0: step_back.name = "#".join(strains)
                #add parent node for next iteration
                new_group.append(step_back)
        #define group for next iteration
        group = new_group

    #after all groups are done, add annotation for base
    if all([i.name for i in tree.root.clades]):
        #find all the names of the children
        names = [i.name for i in tree.root.clades]
        #count GCFs (NCBI strains)
        strains = []
        GCFs = []
        #iterate through children names
        for i in names:
            if i.startswith("GCF") or i == "reference": GCFs.append(i)
            else: strains += i.split("#")
        #if 1 GCF, add to the strains dict
        if len(strains) >= 1 and len(GCFs) >= 1:
            the_name = "#".join(GCFs)
            #add to strains
            if the_name in strains_dict: strains_dict[the_name] += strains
            else: strains_dict[the_name] = strains
    #remove duplicates
    for i in strains_dict: strains_dict[i] = set(strains_dict[i])
    #return
    return strains_dict

def write_strains(strains, workdir):
    """
    Write each strain group to file
    """
    #determine first strain in each strain
    the_strains = ["#".join(i.split("#")[:3]) for i in strains]
    #create set to see if any strain appears twice
    #set_strains = set(the_strains)
    #initialize clone number
    count = 0
    #iterate through clones
    for i in strains:
        #current clone number
        count += 1
        #filename
        name = "strain_%s" % count
        #write
        with open(workdir+"/%s.txt" % name,"w") as outfl:    
            for l in strains[i]:
                outfl.write(l+"\n")

def main(argv=None):
    #Parse command line
    if argv == None: settings = process_command_line(argv)
    #initialize input and output
    #trees
    tree_files = settings.tree
    #tree formats
    tree_formats = settings.format
    #if there is only one format make a list to match all tree files
    if len(tree_formats) == 1: tree_formats = tree_formats*len(tree_files)
    #get comparisons from trees
    trees_comparisons = []
    #iterate through trees
    for tree_file, tree_format in zip(tree_files,tree_formats):
        #define strains
        strains = strains_from_tree(tree_file, tree_format)
        #write to file
        write_strains(strains, settings.workdir)

if __name__ == "__main__":
    exit(main())
