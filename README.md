# WithinHostAdaptation
Repository related to main pipeline in Within Host Adaptation paper (Gatt and Margalit, Mol Biol Evol, in review)

This Repository is comprised of several main src directories:

-Download_data: Encomprises the pipeline to download the data from Table_S1 from the paper and parse it in the user's machine.
-Breseq_pipeline: This directory has two main scripts. Breseq_preTRACE.sh runs breseq for all sequencing libraries of isolates from a single
patient and runs snpEff to assess the genes that are found in each time point. After running the TRACE pipeline for the patient, the user
can proceed to run Breseq_postTRACE.sh, which assesses all changes occuring during the infection.
-TRACE: This directory includes the TRACE program to assess the phylogenetic relationship between the isolates from a patient.
Create_repository downloades complete genome for an organism of choice from NCBI assembly to be used for dividing isolates to clones later.
TRACE_step_I.sh is the first step and can be run with kSNP alone, it divides the isolates to clones. TRACE_step_II.sh can be used immediately after TRACE_step_I.sh, if
the user first ran Breseq_preTRACE.sh the breseq variation data can be used in addition to the kSNP matrix to assess the phylogenetic
relationship between the strains.
-Statistical_framework: This directory includes a single script running the statistical framework described in the paper for all strains
and organisms.

Each src directory has a README.md of its own, describing its use!

Note that sample fastq files are too large for github, and therefore need to be downloaded using the ./download_fastqs.sh script
