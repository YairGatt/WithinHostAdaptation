#! /bin/tcsh
#get which genes are found in each isolate
#input
set SNPEFF_PATH = $1
set pipeline_installation = $2
set patient = $3
set snpEff = $4
set proteome = $5
#iterate
foreach sample_gd ($patient/breseq/*/output/*.gd)
  #get sample vcf
  set sample_vcf = `echo $sample_gd | sed 's/\.gd/.vcf/'`
  #define directory
  set workdir = `dirname $sample_vcf`
  #define output base
  set outbase = `echo $sample_vcf | sed 's/.vcf//'`
  #run snpEff
  java -Xmx4g -jar ${SNPEFF_PATH}/snpEff.jar -no-downstream -no-upstream -s "${outbase}_snpEff_summary.csv" -csvStats $snpEff $sample_vcf > ${outbase}_snpEff.vcf
  #add all genes with no snps
  python $pipeline_installation/Breseq_pipeline/snpEff/fill_summary_up.py ${outbase}_snpEff_summary.genes.txt $proteome
  #extract lost genes
  python $pipeline_installation/Breseq_pipeline/snpEff/parse_snpEff_output.py ${outbase}_snpEff_summary.genes-all.txt
end

