#Statistical framework

-process_input.py
-Breseq_preTRACE.sh
-TRACE_step_I.sh
-TRACE_step_II.sh
-Breseq_postTRACE.sh
-Statistical_framework.sh <-- You are here

Requirements:
-Bio 1.68
-statsmodels 0.8.0
-scipy 1.2.2

Description:
Run Statistical framework for strains from different species to assess genes that were underwent adaptive changes
and loss of function during the infection.

Each organism must have a subdirectory in the main directory, and each strain must have a subdirectory inside the
organism subdirectory in the format of trial_X_patient_X.
All strains must have lists of the genes that were changed and lost during the infection in those
strains, as created by Breseq_postTRACE.sh.

How to run:
Run Statistical_framework.sh:
-Edit "pipeline_installation" variable to location of git on your computer.
-Run:
--------------------------------------------------------------------------------
Statistical_framework.sh $patient_list $reference_proteomes $outdir $temp_dir
--------------------------------------------------------------------------------
-patient list: Either file with full paths of strains, or the super-directory (with subdirectory for each organism etc.)
-reference proteomes: Either file with full paths of proteomes in .faa.gz format, or a super-directory, with a subdirectory for each organism, with a proteome of format .faa.gz inside the subdirectory.
-outdir: Output directory for results
-temp_dir: Directory for temporary files, defult is /tmp

Sample data:
For the users convenience we have attached Sample data under WithinHostAdaptation/Sample_data/Staphylococcus_aureus
Run line:
Statistical_framework.sh WithinHostAdaptation/Sample_data/Staphylococcus_aureus WithinHostAdaptation/Sample_data/genomes/ outdir/
