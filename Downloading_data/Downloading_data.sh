#! /bin/tcsh
#Download data from Supplementary Table 1 in Gatt and Margalit, 2020, for specified organisms.
#Run process_input.py with type sra, and determined paths.
#constants
set pipeline_installation = /home/users/yair/git/WithinHostAdaptation/
set EDIRECT_PATH = /home/users/yair/Software/edirect/
set QUAST_PATH = /home/users/yair/Software/quast-4.5/
set SRATOOLS_PATH = /home/users/yair/Software/sratoolkit/bin/
set SPADES_PATH = /home/users/yair/Software/SPAdes-3.10.1-Linux/bin/
set mail = yair.gatt@mail.huji.ac.il
#input
set outdir = $1
set organism = $2
#determine database location
set database = $pipeline_installation/Downloading_data/data/Table_S1.txt
#run pipeline
python $pipeline_installation/Downloading_data/process_input.py -o $organism -D $database --EDIRECT $EDIRECT_PATH --QUAST $QUAST_PATH --mail $mail --sratools $SRATOOLS_PATH --assembler_installation $SPADES_PATH -t SRA -i $pipeline_installation -w $outdir
