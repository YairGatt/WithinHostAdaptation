#! /bin/tcsh

#This script recieves a file with assemblies and the name of the organism, searches NCBI for all completed genomes for that organism, downloads them and their gff files, and uses the QUAST program to determine how good of a match each genome is for the assemblies.
#It then uses the cosen metric (default is Genome fraction) to rank the genomes by their mean value for the metric for all assemblies, and chooses the best one.

#set variables
set script_path = `dirname $0`
set SPADES_PATH = $1
set EDIRECT_PATH = $2
set QUAST_PATH = $3
set MAIL = $4
set SRATOOLS_PATH = $5
set organism = $6 #scientific name of the organism to be searched for in NCBI
set sras_file = $7 #file with a list of the paths of the assemblies to be assessed
set workdir = $8 #working directory where all final and temportay files will be outputted
set dont_delete = $9 # suppress deletion of temporary files at end of run
#set mail
${EDIRECT_PATH}/econtact -email $MAIL
#create directories
mkdir -p $workdir
mkdir -p ${workdir}/fastqs
mkdir -p ${workdir}/genomes/
mkdir -p ${workdir}/assemblies/
#searchable name of organism
set organism_white = `echo $organism | tr "_" " "`
#determine path for sra-tools
${SRATOOLS_PATH}/vdb-config --set repository/user/main/public/root=${workdir}/ #will downlaod sras to ${workdir}/sra and refseqs to ${workdir}/refseq
#downlaod SRAs from file, also downloads references
${SRATOOLS_PATH}/prefetch --option-file $sras_file >& ${workdir}/sra_tools_output.txt
#define the SRAs properly downloaded
cat ${workdir}/sra_tools_output.txt | grep -v "path not found" | grep -ve '^$' > ${workdir}/sra_tools_results.txt
#make list of SRAs from files
set sra_line = ""
foreach line (`cat $sras_file`)
  if (`cat ${workdir}/sra_tools_results.txt | grep -c "$line"` > 0) then
    set sra_line = `echo $sra_line $line`
  endif
end
#set sra path
set sra_path = $workdir/sra
#create assemblies and run determine_best_genome
set assembly_path = ${workdir}/assemblies/
foreach sra (${sra_path}/*.sra)
  #handle fastqs
  set fastq_path = ${workdir}/fastqs/
  set name = `echo $sra | rev | cut -d '/' -f 1 | rev | sed 's/.sra//'`
  #dump fastqs
  ${SRATOOLS_PATH}/fastq-dump $sra --split-files --gzip -O $fastq_path
  #handle assemblies
  #run SPAdes assembler
  set num_split = `echo ${fastq_path}/${name}*`
  if (${#num_split} == 1) then
    ${SPADES_PATH}/spades.py -s $num_split -o ${assembly_path}/${name}/
  else if (${#num_split} == 2) then
    ${SPADES_PATH}/spades.py -1 ${num_split[1]} -2 ${num_split[2]} -o ${assembly_path}/${name}/
  endif
end
echo "Fastq and Assembly files created"
#Run determine_best_genome for the assemblies
echo ${assembly_path}/*/scaffolds.fasta > ${workdir}/${organism}_assembly_list.txt
${script_path}/determine_best_genome.sh $EDIRECT_PATH $QUAST_PATH $MAIL $organism ${workdir}/${organism}_assembly_list.txt $workdir $dont_delete
set reference_genome = `echo ${workdir}/genome/*.fna.gz`
echo "Determine best genome finished running"
#delete temporary files if all is well
#redefine edirect path
${SRATOOLS_PATH}/vdb-config --restore-defaults
if (`echo $dont_delete | tr '[:upper:]' '[:lower:]'` != "y") then
  #delete sras
  rm -rf ${workdir}/sra/
  #delete refseqs
  rm -rf ${workdir}/refseq/
  #delete records
  rm -f ${workdir}/sra_tools_output.txt
  rm -f ${workdir}/sra_tools_results.txt
  rm -f ${workdir}/${organism}_assembly_list.txt
endif
