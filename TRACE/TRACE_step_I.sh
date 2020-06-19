#! /bin/tcsh

set pipeline_installation = /home/users/yair/git/WithinHostAdaptation
set kSNP_INSTALLATION = /home/users/yair/Software/kSNP3.1_Linux_package/

#set directory with isolates,  with a subdirectory for each assembly (the assembly itself can either be called scaffolds.fasta or anything.fna)
set input_directory = $1
if (! -d $input_directory) then
  echo "$input_directory does not exist"
endif
#set repository file
set repository = $2
if (! -e $repository) then
  echo "$repository is not a file"
endif
#set directory to write strains
set outdir = $3
mkdir -p $outdir
#set file with names of samples for conversion. format is specified in README
set convert_file = $4
if ("$convert_file" == "") then
  set convert_file == ${input_directory}/convert_file.txt
endif
set temp_dir = $5
if ("$temp_dir" == "") then
  set temp_dir = /tmp/TRACE/
  mkdir -p $temp_dir
endif
#create tree
#run kSNP
python $pipeline_installation/TRACE/bin/step_I/create_strain.py $input_directory $repository $temp_dir $convert_file $outdir $kSNP_INSTALLATION/kSNP3
#divide to strains
if (-s ${outdir}/tree.ML.tre) then
  python $pipeline_installation/TRACE/bin/step_I/divide_to_strains.py -t ${outdir}/tree.ML.tre -f newick -w $outdir
else if (-s ${outdir}/tree.core.tre) then
  python $pipeline_installation/TRACE/bin/step_I/divide_to_strains.py -t ${outdir}/tree.core.tre -f newick -w $outdir
else
  echo "Tree could not be constructed for $patient."
endif
#convert to names
foreach strain (${outdir}/strain_*.txt)
  python $pipeline_installation/TRACE/bin/Utilities/convert_to_names.py $strain $convert_file
end
