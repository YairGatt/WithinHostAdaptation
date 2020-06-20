#Downloading data

-Downloading_data.sh <-- You are here
-Breseq_preTRACE.sh
-Create_repository.sh
-TRACE_step_I.sh
-TRACE_step_II.sh
-Breseq_postTRACE.sh
-Statistical_framework.sh

Requirements:
-Bio 1.68
-SRA Toolkit 2.9.1
-EDirect 9.50
-Quast 4.5
-SPAdes 3.10.1

Description:
Download fastq data for chosen organism from Supplementary Table S1 from Gatt and Margalit, 2020, assemble to assemblies,
quality filter, and process it to relevant formats.

How to run:
Run Downloading_data.sh:
-Edit "pipeline_installation" variable to location of git on your computer.
-Edit EDIRECT_PATH variable to location of edirect on your computer.
-Edit QUAST_PATH variable to location of quast on your computer.
-Edit SRATOOLS_PATH to location of SRA Tools bin directory on your computer.
-Edit SPADES_PATH to location of SPAdes bin directory on your computer.
-Edit mail to your NCBI mail (for use in authentication when using EDirect)

-Run:
--------------------------------------------------------------------------------
Downloading_data.sh $outdir $organism
--------------------------------------------------------------------------------
-outdir: Output directory for results, it is recommended that the user create a subdirectory for each organism in his main
directory, and use that as the output directory.
-organism: Scientific name of organism whose data is to be processed.

Sample data:
Table S1 is provided in WithinHostAdaptation/Downloading_data/data/Table_S1.txt
We give of Enterococcus faecalis (which only has data from a single patient) as an example run:
Run line:
Downloading_data.sh outdir/ Enterococcus_faecalis
