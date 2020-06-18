#! /bin/tcsh

#This script recieves the patient directory as downloaded by Download_data.sh and a directory with a gbk format genome,
#and a _protein.faa.gz format proteome
#It then runs breseq for all assemblies in that patient

#initialize
set patient = $1
set genome_dir = $2
set snpEff_name = $3

set genome = $genome_dir/*_genomic.gbk
set proteome = $genome_dir/*_protein.faa.gz

set BRESEQ_DIR = /home/users/yair/Software/breseq-0.32.1-Linux-x86_64/bin/
set SNPEFF_PATH = /home/users/yair/Software/snpEff_original/
set pipeline_installation = /home/users/yair/git/WithinHostAdaptation

#create directories
mkdir -p ${patient}/breseq
mkdir -p ${patient}/breseq/comparative
#run breseq for all fastq files
foreach fastq (${patient}/fastqs/*_1.fastq.gz)
  set fastq_base = `echo $fastq | sed 's/_1.fastq.gz//'`
  set output_dir = `basename $fastq_base`
  #create output directory for fastq breseq
  mkdir -p ${patient}/breseq/${output_dir}
  #run breseq
  ${BRESEQ_DIR}/breseq -j 8 --brief-html-output --no-javascript -r $genome -o ${patient}/breseq/${output_dir} ${fastq_base}*
  #change name and convert to vcf
  mv ${patient}/breseq/${output_dir}/output/output.gd ${patient}/breseq/${output_dir}/output/${output_dir}.gd
  ${BRESEQ_DIR}/gdtools GD2VCF -r $genome -o ${patient}/breseq/${output_dir}/output/${output_dir}.vcf ${patient}/breseq/${output_dir}/output/${output_dir}.gd
  #get some statistics
  ${BRESEQ_DIR}/gdtools COUNT -r $genome -o ${patient}/breseq/${output_dir}/output/${output_dir}_counts.csv ${patient}/breseq/${output_dir}/output/${output_dir}.gd
  #remove heavy temp files
  rm ${patient}/breseq/${output_dir}/data/reference.bam
  rm ${patient}/breseq/${output_dir}/08_mutation_identification/*.coverage.tab
  #determine mapping percentage
  set reads = `cat "${patient}/breseq/${output_dir}/output/${output_dir}_counts.csv" | tail -1 | cut -f 4 -d ","`
  set mapped_reads = `cat "${patient}/breseq/${output_dir}/output/${output_dir}_counts.csv" | tail -1 | cut -f 8 -d ","`
  set mapped_percentage = `awk "BEGIN {print 100*$mapped_reads/$reads}" | cut -d "." -f 1`
  if ("$mapped_percentage" < 80) then
    echo "$output_dir has less than 80% mapping percentage: ${mapped_percentage}%"
    continue
  endif
end
#run snpEff
$pipeline_installation/Breseq_pipeline/snpEff/snpEff_patient.sh $SNPEFF_PATH $pipeline_installation $patient $snpEff_name $proteome
