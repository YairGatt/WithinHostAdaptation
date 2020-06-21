#! /bin/tcsh

set QUAST_PATH = $1
set patient = $2
set reference = $3
set threshold = $4
set skip = $5

set a = `ls $patient | grep -c "assemblies"`

if ("$a" != "0") then
  mkdir -p ${patient}/below_threshold/
  mkdir -p ${patient}/below_threshold/fastqs
  mkdir -p ${patient}/below_threshold/assemblies
  set assemblies = ${patient}/assemblies/*
else
  set assemblies = $patient
endif

if ($skip != "Y") then
  foreach assembly ($assemblies)
    mkdir -p ${assembly}/quast
    set first_choice = `ls ${assembly} | grep -c "\.fna"` #first choice means a preassembled genome was downloaded
    #scaffolds
    if (! $first_choice) then
      set strain = `basename $assembly`
      set num_split = `echo ${patient}/fastqs/${strain}*`
      if (${#num_split} == 2) then
        $QUAST_PATH/quast.py -R $reference -o ${assembly}/quast $assembly/scaffolds.fasta --threads 4 --scaffolds --reads1 ${patient}/fastqs/${strain}_1.fastq.gz --reads2 ${patient}/fastqs/${strain}_2.fastq.gz > /dev/null
      else if (${#num_split} == 1) then
        $QUAST_PATH/quast.py -R $reference -o ${assembly}/quast $assembly/scaffolds.fasta --threads 4 --scaffolds --reads1 ${patient}/fastqs/${strain}_1.fastq.gz > /dev/null
      else if (${#num_split} == 0) then
        $QUAST_PATH/quast.py -R $reference -o ${assembly}/quast $assembly/scaffolds.fasta --threads 4 > /dev/null
      else
        echo "inapproriate number of fastq files for $assembly"
      endif
    else
      set first_choice = ${assembly}/*.fna
      $QUAST_PATH/quast.py -R $reference -o ${assembly}/quast $first_choice --threads 4 --space-efficient --fast > /dev/null
    endif
  end
endif

foreach assembly (${patient}/assemblies/*)
  set strain = `basename $assembly`
  if (-d ${assembly}/quast) then
    set genome_fraction = `cat ${assembly}/quast/report.tsv | grep "Genome fraction" | cut -f 2 | cut -d "." -f 1`
    echo $strain $genome_fraction
    if ($genome_fraction < $threshold) then
      mv $assembly ${patient}/below_threshold/assemblies
      mv ${patient}/fastqs/${strain}_* ${patient}/below_threshold/fastqs
    endif
  endif
end
