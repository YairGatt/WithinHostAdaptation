#Statistical framework

-process_input.py
-Breseq_preTRACE.sh
-Create_repository.sh <-- You are here
-TRACE_step_I.sh <-- You are here
-TRACE_step_II.sh <-- You are here
-Breseq_postTRACE.sh
-Statistical_framework.sh

Requirements:
-kSNP 3.1
-Bio 1.68
-numpy 1.16.4
-matplotlib 2.0.2

Description:
Run TRACE for isolates of strains to assess the phylogenetic relationship between them.
Isolates should be in a single directory, with a subdirectory for each isolate, including either a
scaffolds.fasta file, or an assembly ending in .fna.

If needed, first create a repository: list of complete genomes from NCBI to be used for dividing isolates to clones.

How to run:
Run Create_repository.sh for organism:
-Edit kSNP_INSTALLATION variable to location of kSNP on your computer.
-Run:
--------------------------------------------------------------------------------
Create_repository.sh $organism $repository $reference $skip $REPOSITORY_DIR
--------------------------------------------------------------------------------
-organism: Scientific name of organism for repository.
-repository: File with list of genome paths
-reference: reference genome for rooting trees
-skip: If "Y", skip downloading and rewrite repository file only
-REPOSITORY_DIR: Directory to download genomes to

Run TRACE_step_I:
-Edit "pipeline_installation" variable to location of git on your computer.
-Edit kSNP_INSTALLATION variable to location of kSNP on your computer.
-Run:
--------------------------------------------------------------------------------
TRACE_step_I.sh $input_directory $repository $outdir $convert_file $temp_dir
--------------------------------------------------------------------------------
-input_directory: Directory with subdirectory for each isolate, with a scaffolds.fasta file or .fna file
-repository: List of paths of complete genomes for the organism downloaded from NCBI assembly (by Create_repository.sh)
-outdir: Output directory for results
-convert_file: File converting assembly names to times, see Sample data for format (created by process_input.py)
-temp_dir: Directory for temporary files, defult is /tmp

Run TRACE_step_II:
-Edit "pipeline_installation" variable to location of git on your computer.
-Edit kSNP_INSTALLATION variable to location of kSNP on your computer.
-Run:
--------------------------------------------------------------------------------
TRACE_step_II.sh $strains_dir $assemblies_dir $convert_file $workdir $repository $temp_dir $breseq_dir $core
--------------------------------------------------------------------------------
-strains_dir: Directory with results of TRACE_step_I.sh
-assemblies_dir: Directory with subdirectory for each isolate, with a scaffolds.fasta file or .fna file
-convert_file: File converting assembly names to times, see Sample data for format (created by process_input.py)
-workdir: Output directory for results
-repository: List of paths of complete genomes for the organism downloaded from NCBI assembly (by Create_repository.$
-temp_dir: Directory for temporary files, defult is /tmp
-breseq_dir: Directory with results of Breseq_preTRACE.sh (optional!)
-core: put Y to use only core snps (--core in kSNP)

Sample data:
Sample data is provided in WithinHostAdaptation/TRACE/Sample_data
Run line:
TRACE_step_I.sh WithinHostAdaptation/TRACE/Sample_data/assemblies WithinHostAdaptation/TRACE/Sample_data/repository/Staphylococcus_aureus.txt outdir WithinHostAdaptation/TRACE/Sample_data/convert_file/Convert_file.txt
TRACE_step_II.sh outdir WithinHostAdaptation/TRACE/Sample_data/assemblies WithinHostAdaptation/TRACE/Sample_data/convert_file/Convert_file.txt outdir1 WithinHostAdaptation/TRACE/Sample_data/repository/Staphylococcus_aureus.txt /tmp/TRACE/ WithinHostAdaptation/TRACE/Sample_data/breseq/
