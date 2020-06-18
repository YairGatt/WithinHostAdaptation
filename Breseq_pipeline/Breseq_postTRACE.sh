#! /bin/tcsh

#This script runs on a patient on which Breseq_preTRACE.sh and TRACE have already been ran.
#It compares breseq output between progenitor-progeny pairs determined by TRACE.
#It recieves a patient direcotry, a gbk format genome and the snpEff name

#initialize
set patient = $1
set genome = $2
set tree_dir = $3
set snpEff_name = $4
#constants
set BRESEQ_DIR = /home/users/yair/Software/breseq-0.32.1-Linux-x86_64/bin/
set SNPEFF_PATH = /home/users/yair/Software/snpEff_original/
set pipeline_installation = /home/users/yair/git/WithinHostAdaptation
#check for breseq
if (! -d $patient/breseq || "$patient" == "/home/hosts/disk20/projects_data/project_B/data/Klebsiella_pneumoniae/trial_54_patient_1" || "$patient" == "/home/hosts/disk20/projects_data/project_B/data/Klebsiella_pneumoniae/trial_54_patient_10") then
  echo "$patient has no breseq data, please run Breseq_preTRACE first."
  exit 1
endif
#check for trees
if ("$tree_dir" == "") then
  set tree_dir = ${patient}/trees/
endif
if (! -e $tree_dir/formatted_comparisons_converted.txt) then
  echo "$patient has no tree"
  exit 1
else
  set a = `cat $tree_dir/formatted_comparisons_converted.txt | grep -vc "|"`
endif
if ("$a" == "0") then
  echo "$patient only has internodes in tree"
  exit 1
endif
#compare results for different genomes according to tree
#initialize results vectors
set comparisons = ""
#iterate through all comparisons from tree
foreach line ( "`cat $tree_dir/formatted_comparisons_converted.txt`")
  #take only lines that are a comparison between two different isolates
  set good_line = `echo $line | awk '$1 != $2 {print $0}'`
  set fictious = `echo $line | grep -c "fiction"`
  if ("$fictious" > 0) then
    continue
  endif
  if ("$good_line" != "") then
    set timepoint1 = ${good_line[1]}
    set timepoint2 = ${good_line[2]}
    #create results directory
    set workdir = "${patient}/breseq/comparative/${timepoint1}#${timepoint2}/"
    mkdir -p "$workdir"
    mkdir -p "${workdir}/temp_files"
    mkdir -p "${workdir}/results"
    #run BRESEQ compare
    ${BRESEQ_DIR}/gdtools COMPARE -r $genome -f GD -o "${workdir}/temp_files/diff.gd"  ${patient}/breseq/${timepoint1}/output/${timepoint1}.gd ${patient}/breseq/${timepoint2}/output/${timepoint2}.gd
    set times = `python $pipeline_installation/Breseq_pipeline/bin/fix_timepoints.py "${workdir}/temp_files/diff.gd" $timepoint1 $timepoint2`
    #keep only differential lines
    cat "${workdir}/temp_files/diff.gd" | grep -E "GENOME_DIFF|frequency_${times[1]}=0" > "${workdir}/temp_files/straight_diff.gd"
    cat "${workdir}/temp_files/diff.gd" | grep -E "GENOME_DIFF|frequency_${times[2]}=0" > "${workdir}/temp_files/reverse_diff.gd"
    #convert to VCF
    ${BRESEQ_DIR}/gdtools GD2VCF -r $genome -o "${workdir}/temp_files/straight_diff.vcf" "${workdir}/temp_files/straight_diff.gd"
    ${BRESEQ_DIR}/gdtools GD2VCF -r $genome -o "${workdir}/temp_files/reverse_diff.vcf" "${workdir}/temp_files/reverse_diff.gd"
    #running snpEff
    set organism_dir = `dirname $genome`
    set organism = `basename $organism_dir`
    #run it
    $pipeline_installation/Breseq_pipeline/snpEff/snpEff_pair.sh $SNPEFF_PATH $pipeline_installation $snpEff_name "${workdir}/temp_files/straight_diff.vcf" "${workdir}/temp_files/reverse_diff.vcf" "${workdir}/temp_files/"
    #add to list of results
    set comparisons = `echo "$comparisons" "${workdir}/"`
  endif
end
#check if there are comparisons
if ("$comparisons" == "") then
  echo "No good comparisons found for patient $patient"
  exit 1
endif
#get final lost_genes files for combination
set changed_files = ""
set high_lost_files = ""
foreach comparison ($comparisons)
  set changed_files = `echo $changed_files ${comparison}/changed_genes_full.txt`
  set high_lost_files = `echo $high_lost_files ${comparison}/high_lost_genes_full.txt`
end
#changed
python /home/users/yair/Documents/PhD_projects/project_B/bin/gene_loss_pipeline/breseq_pipeline/combine_comparisons_for_final_patient_matrix.py -l $changed_files -c ${patient}/trees/comparisons_converted.txt -o ${patient}/changed_genes.txt
#high lost
python /home/users/yair/Documents/PhD_projects/project_B/bin/gene_loss_pipeline/breseq_pipeline/combine_comparisons_for_final_patient_matrix.py -l $high_lost_files -c ${patient}/trees/comparisons_converted.txt -o ${patient}/high_lost_genes.txt
