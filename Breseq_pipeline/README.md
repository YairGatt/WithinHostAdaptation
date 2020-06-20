#Breseq pipeline

-Downloading_data.sh
-Breseq_preTRACE.sh <-- You are here
-Create_repository.sh
-TRACE_step_I.sh
-TRACE_step_II.sh
-Breseq_postTRACE.sh <-- You are here
-Statistical_framework.sh

Requirements:
-breseq 0.32.1
-gdtools 0.32.1
-snpEff 4.1
-Bio 1.68
-SRA Toolkit 2.9.1 (only for downloading sample data)

Description:
Run Breseq for all isolates in each strain to detect differences from reference genome, and then run again after
determining progenitor-progeny pairs to determine differences between Breseq results of different isolates, and
assess which occured during the infection.

How to run:
Run Breseq_preTRACE.sh for each strain:
-Edit BRESEQ_DIR variable to location of breseq on your computer.
-Edit SNPEFF_PATH variable to location of snpEff on your computer.
-Edit "pipeline_installation" variable to location of git on your computer.
-Run:
--------------------------------------------------------------------------------
Breseq_preTRACE.sh $patient $genome_dir $snpEff_name
--------------------------------------------------------------------------------
-patient: Directory created by process_input.py for patient. Should have subdirectory names fastqs with sequencing results in fastq format.
-genome_dir: Directory with .gbk format reference genome and .faa.gz format proteome (can be downloaded from NCBI)
-snpEff_name: Name of reference genome in snpEff.

After running TRACE proceed:

Run Breseq_postTRACE.sh for each strain:
-Edit BRESEQ_DIR variable to location of breseq on your computer.
-Edit SNPEFF_PATH variable to location of snpEff on your computer.
-Edit "pipeline_installation" variable to location of git on your computer.
-Run:
--------------------------------------------------------------------------------
Breseq_postTRACE.sh $patient $genome_dir $snpEff_name $tree_dir
--------------------------------------------------------------------------------
-patient: Directory created by process_input.py for patient. Should have the results of Breseq_preTRACE.sh
in a subdirectory named breseq.
-genome_dir: Directory with .gbk format reference genome and .faa.gz format proteome (can be downloaded from NCBI)
-snpEff_name: Name of reference genome in snpEff.
-tree_dir: TRACE_step_II.sh results

Sample data:
Sample data is provided in WithinHostAdaptation/Sample_data/Staphylococcus_aureus/trial_122_patient_4
Sample fastqs are too heavy for github and should first be downloaded using WithinHostAdaptation/download_fastqs.sh
How to run download_fastqs.sh:
-Edit "pipeline_installation" variable to location of git on your computer.
-Edit SRATOOLS_PATH variable to location of sratoolkit on your computer.
--------------------
download_fastqs.sh
--------------------
Run line:
Breseq_preTRACE.sh WithinHostAdaptation/Sample_data/Staphylococcus_aureus/trial_122_patient_4/ WithinHostAdaptation/Sample_data/genomes/Staphylococcus_aureus Staphylococcus_aureus_NCTC_8325_uid57795
Breseq_postTRACE.sh WithinHostAdaptation/Sample_data/Staphylococcus_aureus/trial_122_patient_4/ WithinHostAdaptation/Sample_data/genomes/Staphylococcus_aureus Staphylococcus_aureus_NCTC_8325_uid57795 WithinHostAdaptation/Sample_data/Staphylococcus_aureus/trial_122_patient_4/trees/
