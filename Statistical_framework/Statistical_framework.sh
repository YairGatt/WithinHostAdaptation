#! /bin/tcsh

#This script runs the statistical test for a list of patients downloaded by Download_DB.sh and processed by Gene_loss_pipeline.sh
#It recieves a list of the patients or a directory including a subdirectory for each species with subdirectories for each patient infected by that species
#Patients can also be given in the file as trials, for instance trial_15 will include all patients in the format of trial_15_patient_X
#Also requires list of reference proteomes or directory with subdirectories for each organism including the reference proteomes in the format of *_protein.faa (fasta format)

set pipeline_installation = /home/users/yair/git/WithinHostAdaptation

#input
set patient_list = $1
set reference_proteomes = $2
set outdir = $3
set temp_dir = $4
#determine if list is directory or file
if (-d "$patient_list") then
  set all_patients = `ls -d $patient_list/*/*`
else if (-e "$patient_list") then
  set all_patients = `cat $patient_list | grep -v "#"`
endif
#determine if proteome list is directory or file
if (-d "$reference_proteomes") then
  set proteomes = `ls -d $reference_proteomes/*/*_protein.faa.gz`
else if (-e "$reference_proteomes") then
  set proteomes = `cat $reference_proteomes | grep -v "#"`
endif
#temp dir
if ("$temp_dir" == "") then
  set temp_dir = /tmp/WithinHostAdaptation/
endif
mkdir -p $temp_dir
#determine temp files and prepare data for run
set organism_file = $temp_dir/organisms.txt
set changed_genes_file = $temp_dir//changed_genes.txt
set high_lost_genes_file = $temp_dir/high_lost_genes.txt
set millstone_files_file = $temp_dir/millstone_files.txt
#delete files if exist
rm -f $organism_file
rm -f $changed_genes_file
rm -f $high_lost_genes_file
rm -f $millstone_files_file
#iterate through different patients, they can also be given as trials
foreach trial ($all_patients)
  set trialname = `basename $trial`
  #get organism
  set organism = `dirname $trial`
  set the_organism = `basename $organism`
  #determine if single patient or trial
  set a = `echo $trialname | grep -c "patient"`
  if ("$a" == 0) then #trial
    #get all patients
    set data = `ls ${trial}_patient_*/ | grep -c "_genes.txt"`
    if ("$data" != 0) then
      #get changed and loss of function gene lists and add to temp file
      ls ${trial}_patient_*/changed_genes.txt >> $changed_genes_file
      ls ${trial}_patient_*/high_lost_genes.txt >> $high_lost_genes_file
      echo $the_organism >> $organism_file
      ls ${trial}_patient_*/millstone_file.txt >> $millstone_files_file
    endif
  else #single paitnet
    set data = `ls ${trial}/ | grep -c "_genes.txt"`
    if ("$data" != 0) then
      ls $trial/changed_genes.txt >> $changed_genes_file
      ls $trial/high_lost_genes.txt >> $high_lost_genes_file
      echo $the_organism >> $organism_file
      ls $trial/millstone_file.txt >> $millstone_files_file
    endif
  endif
end
#iterate through organisms, read results and run test
foreach organism (`cat $organism_file | sort -u`)
  set changed_genes_organism = `cat $changed_genes_file | grep $organism`
  set high_lost_genes_organism = `cat $high_lost_genes_file | grep $organism`
  set millstone_files_organism = `cat $millstone_files_file | grep $organism`
  #defien reference proteome
  set reference_proteome = `ls $proteomes | grep $organism`
  #run patient matrix
  python $pipeline_installation/Statistical_framework/bin/patient_matrix.py -l $changed_genes_organism -o $outdir/changed_${organism}.csv -p $reference_proteome -i $millstone_files_organism
  python $pipeline_installation/Statistical_framework/bin/patient_matrix.py -l $high_lost_genes_organism -o $outdir/high_${organism}.csv -p $reference_proteome -i $millstone_files_organism
end
