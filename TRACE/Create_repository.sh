#! /bin/tcsh

#This script creates a repository file of NCBI complete genomes that can be utilized by TRACE to divide assemblies from a strain to clones

#determine kSNP installation
set kSNP_INSTALLATION = /home/users/yair/Software/kSNP3.1_Linux_package/
#input
#organism of interest
set organism = $1
#output file for list of NCBI compelte genomes
set output = $2
#reference genome to be used for rooting future trees
set reference = $3
#determine whether to skip downlaoding all the genomes and just rewrite the repository file
set skip = $4
#determine repository dir
set REPOSITORY_DIR = $5
#go
if ("$skip" != "Y") then
  #download all genoems and assembly stats files from NCBI
  rsync --recursive --copy-links --times --verbose --exclude "*/*_from_*" --include "*/*_genomic.fna.gz" --include "*/*assembly_stats.txt" --exclude "*/*" --prune-empty-dirs rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/${organism}/latest_assembly_versions/ ${REPOSITORY_DIR}/${organism}
  #remove non complete genomes
  foreach stats ($REPOSITORY_DIR/${organism}/*/*_assembly_stats.txt)
    set genome = `dirname $stats`
    set level = `cat $stats | grep "Assembly level:" | cut -d ":" -f 2 | sed 's/ //'`
    if ("$level" !~ "Complete Genome*") then
      rm -rf $genome
    endif
  end
  #remove assembly stats files
  rm -f $REPOSITORY_DIR/${organism}/*/*_assembly_stats.txt
  #correct names if they have dots and such
  foreach organism_folder ($REPOSITORY_DIR/${organism}/*)
    set name = `basename $organism_folder | sed 's/\./-/g'`
    mv $organism_folder $REPOSITORY_DIR/${organism}/${name}
  end
  #gunzip all genomes for kSNP
  gunzip $REPOSITORY_DIR/${organism}/*/*.fna.gz
  #iterate through genomes and adjust name
  foreach organism_genome ($REPOSITORY_DIR/${organism}/*/*.fna)
    set dir = `dirname $organism_genome`
    set name = `basename $dir`
    mv $organism_genome ${dir}/${name}_genomic.fna
  end
endif
#create kSNP reference direcotry if doesn't exist
mkdir -p $kSNP_INSTALLATION/references/
#copy to kSNP reference directory
cp $reference $kSNP_INSTALLATION/references/${organism}_reference.fna.gz
gunzip $kSNP_INSTALLATION/references/${organism}_reference.fna.gz
#ad to repository on top
echo "$kSNP_INSTALLATION/references/${organism}_reference.fna	reference" > $output
#add all complete genomes to repository file
foreach the_organism ($REPOSITORY_DIR/${organism}/*/*.fna)
  set name = `basename $the_organism | sed 's/\.fna//' | tr '-' '_'`
  echo "$the_organism	$name" >> $output
end
