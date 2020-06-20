"""
This script recieves a file with assemblies and the name of the organism, searches NCBI for all completed genomes for that organism, downloads them and their gff files, and uses the QUAST program to determine how good of a match each genome is for the assemblies.
It then uses the cosen metric (default is Genome fraction) to rank the genomes by their mean value for the metric for all assemblies, and chooses the best one.
"""

import sys,os
from subprocess import call
import warnings
import glob
import optparse

def process_command_line(argv):
    """
    Return an args list
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]
        
    # initialize the parser object:
    parser = optparse.OptionParser(
             formatter=optparse.TitledHelpFormatter(width=78),
             add_help_option=None)

    # define options here:

    parser.add_option(
        '-w', '--workdir',
        help='Workdir where temporary and final files will be saved.')

    parser.add_option(
        '-i', '--sras_file',
        help='File with a list of sras for which a reference genome is to be determined.')

    parser.add_option(
        '--installation',
        help='Pipeline installation.')

    parser.add_option(
           '--EDIRECT',
           help='edirect tools installation.')

    parser.add_option(
           '--QUAST',
           help='Quast installation.')

    parser.add_option(
           '--mail',
           help='Email for edirect.')

    parser.add_option(
           '--sratools',
           help='SRAtools installation.')

    parser.add_option(
        '-o', '--organism',
        help='Organism to be searched for on NCBI Assembly.')

    parser.add_option(
        '--dont_delete', action="store_true",
        help='Do not delete temporary files after running.')
    
    parser.add_option(
        '--representative', action="store_true",
        help='Download only representative genomes.')
    
    parser.add_option(
        '-s', '--script', default="/home/users/yair/Documents/PhD_projects/project_B/bin/downloading_database/fastq_determine_best_genome.sh",
        help='Path of determine_best_genome.sh script')
    
    parser.add_option(
        '--assembler_installation', default="/home/users/yair/Software/SPAdes-3.10.1-Linux/bin/",
        help='Path of assembler installation')
            
    parser.add_option(    # customized description; put --help last
        '-h', '--help', action='help',
        help='Show this help message and exit.')

    settings, args = parser.parse_args(argv)

    return settings, args

def process_call_from_script(argv):
    command_line_settings, args = process_command_line(["-s","Downloading_data/fastq_determine_best_genome.sh"])
    for i in dir(command_line_settings):
        try: getattr(argv, i)
        except AttributeError: setattr(argv, i, getattr(command_line_settings, i))
    return argv

def main(argv=None):
    #Parse command line
    if argv == None: settings, args = process_command_line(argv)
    else: settings = process_call_from_script(argv)
    #determine deletion
    try:
        if settings.dont_delete: no_deletion = "Y"
        else: no_deletion = "N"
    except AttributeError: no_deletion = "N"
    #run script
    call([settings.installation + settings.script, settings.assembler_installation, settings.EDIRECT, settings.QUAST, settings.mail, settings.sratools, settings.organism, settings.sras_file, settings.workdir, no_deletion])
    return 1

if __name__ == "__main__":
        exit(main())
