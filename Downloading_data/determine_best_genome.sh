#! /bin/tcsh

#This script recieves a file with assemblies and the name of the organism, searches NCBI for all completed genomes for that organism, downloads them and their gff files, and uses the QUAST program to determine how good of a match each genome is for the assemblies.
#It then uses the cosen metric (default is Genome fraction) to rank the genomes by their mean value for the metric for all assemblies, and chooses the best one.

#set variables
set script_path = `dirname $0`
set EDIRECT_PATH = $1
set QUAST_PATH = $2
set MAIL = $3
set organism = $4 #scientific name of the organism to be searched for in NCBI
set assemblies_file = $5 #file with a list of the paths of the assemblies to be assessed
set workdir = $6 #working directory where all final and temportay files will be outputted
set dont_delete = $7 # suppress deletion of temporary files at end of run
#set mail
${EDIRECT_PATH}/econtact -email $MAIL
#create directories
mkdir -p $workdir
mkdir -p ${workdir}/
mkdir -p ${workdir}/genome
mkdir -p ${workdir}/quast
mkdir -p ${workdir}/assemblies
mkdir -p ${workdir}/genomes
mkdir -p ${workdir}/quasts
#get organism name in searchable format
set organism_white = `echo $organism | tr "_" " "`
#flag representing no choice regarding reference genome (as only one is available)
set no_choice = 0
#make list of assemblies from file
foreach line (`cat $assemblies_file`)
  if (! -s $line) then
    ${EDIRECT_PATH}/esearch -db assembly -query "$line AND latest [filter]" | ${EDIRECT_PATH}/efetch -format docsum | ${EDIRECT_PATH}/xtract -pattern DocumentSummary -block FtpPath -match "@type:Genbank" -element FtpPath | awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz";filesuffix2="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;file2=asm"_"filesuffix2;print ftpdir"/"file;print ftpdir"/"file2}' >> ${workdir}/${organism}_assembly_list.txt
  endif
end
#download
perl $script_path/rsyn_files.pl ${workdir}/${organism}_assembly_list.txt ${workdir}/assemblies > /dev/null
#iterate through assemblies
foreach assembly (${workdir}/assemblies/*.fna.gz)
  set name = `basename $assembly | sed 's/_genomic.fna.gz//'`
  mkdir -p ${workdir}/assemblies/${name}
  #mv assembly
  mv $assembly ${workdir}/assemblies/${name}/
  #gunzip assembly
  gunzip ${workdir}/assemblies//${name}/${name}_genomic.fna.gz
  #mv gff
  set gff = `echo $assembly | sed 's/.fna.gz/.gff.gz/'`
  mv $gff ${workdir}/assemblies/${name}/
end
echo "Assemblies file parsed"
#download relevent genomes and gff files from ncbi
${EDIRECT_PATH}/esearch -db assembly -query "$organism_white [ORGN] AND latest [filter] AND "+'"complete genome"'+" [filter] AND representative [PROP]" | ${EDIRECT_PATH}/efetch -format docsum | ${EDIRECT_PATH}/xtract -pattern DocumentSummary -block FtpPath -match "@type:Genbank" -element FtpPath | awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz";filesuffix2="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;file2=asm"_"filesuffix2;print ftpdir"/"file;print ftpdir"/"file2}' > ${workdir}/${organism}_ftp_list.txt
#check results
if (`cat ${workdir}/${organism}_ftp_list.txt | wc -l` == 0) then
  echo "Error: no genomes found for ${organism_white} in NCBI Assembly"
  exit 110
else if (`cat ${workdir}/${organism}_ftp_list.txt | wc -l` < 3 && `cat ${workdir}/${organism}_ftp_list.txt | wc -l` > 0 ) then
  echo "One genome found on NCBI Assembly"
  set no_choice = 1
else
  echo "Relevant genomes found on NCBI Assembly"
endif
#rsync script from NCBI rsyn_files.pl
perl rsyn_files.pl ${workdir}/${organism}_ftp_list.txt ${workdir}/genomes > /dev/null
set genome_dir = ${workdir}/genomes/
#define assemblies
set assemblies = `echo ${workdir}/assemblies/*/*.fna ${workdir}/assemblies/*/scaffolds.fasta`
#select quast stat
set stat = "Genome fraction"
set stat_nowhite = `echo $stat | sed 's/ /_/'`
if ("$no_choice" == 0) then
  #run quast for each set of assemblies with each reference genome
  foreach line (${genome_dir}/*.fna.gz)
    set genome = `echo $line | sed 's/.fna.gz//'`
    set genome_name = `echo $line | rev | cut -d '/' -f 1 | rev`
    /home/users/yair/Software/quast-4.5/quast.py -o ${workdir}/quasts/${genome_name} -R $line -G ${genome}.gff.gz $assemblies --threads 4 --space-efficient --fast > /dev/null
  end
  set quast_dir = ${workdir}/quasts/
  #list stat for all genomes
  rm -f ${workdir}/quasts_${stat_nowhite}.txt
  foreach file (${quast_dir}/*/report.tsv)
    set genome = `echo $file | rev | cut -d '/' -f 2 | rev`
    grep "$stat" $file | awk -v genome=$genome -F "\t" '{for (i=2; i<=NF; i++) sum+=$i} END {print genome "\t" sum/(NF-1)}' >> ${workdir}/quasts_${stat_nowhite}.txt
  end
  #determine best genome
  set reference_genome = `cat ${workdir}/quasts_${stat_nowhite}.txt | sort -k 2 -nr | head -1 | cut -f 1`
else
  set reference_genome = `basename ${genome_dir}/*.fna.gz`
  set quast_dir = ${workdir}/quasts/
endif
#save results
set clean_genome = `echo $reference_genome | sed 's/.fna.gz//'`
#move to own directory
cp -f ${genome_dir}/${clean_genome}* ${workdir}/genome/
if (-d ${quast_dir}/${reference_genome}/) then
  cp -f ${quast_dir}/${reference_genome}/* ${workdir}/quast/ -r
endif
#clear bad quality genomes
${script_path}/quast_patient_combined.sh $QUAST_PATH ${workdir} ${workdir}/genome/*.fna.gz 50 N
#delete temporary files if all is well
if (`echo $dont_delete | tr '[:upper:]' '[:lower:]'` != "y") then
  #delete ftp list
  rm ${workdir}/${organism}_ftp_list.txt
  #delete list of genomes with ranks
  rm -f ${workdir}/quasts_${stat_nowhite}.txt
  #delete other genomes
  rm -rf ${workdir}/genomes
  #delete other quasts
  rm -rf ${workdir}/quasts
endif
