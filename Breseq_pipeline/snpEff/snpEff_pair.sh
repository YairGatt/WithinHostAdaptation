#! /bin/tcsh

#input
set SNPEFF_PATH = $1
set pipeline_installation = $2
set genome = $3
set straight_vcf = $4
set reverse_vcf = $5
set outbase = $6
set upper = `dirname $outbase`

#run snpEff
java -Xmx4g -jar ${SNPEFF_PATH}/snpEff.jar -no-downstream -no-upstream -s "${outbase}snpEff_straight_summary.csv" -csvStats $genome $straight_vcf > "${outbase}snpEff_straight.vcf"
java -Xmx4g -jar ${SNPEFF_PATH}/snpEff.jar -no-downstream -no-upstream -s "${outbase}snpEff_reverse_summary.csv" -csvStats $genome $reverse_vcf > "${outbase}snpEff_reverse.vcf"
#extract lost genes
python $pipeline_installation/Breseq_pipeline/snpEff/parse_snpEff_output.py ${outbase}snpEff_straight_summary.genes.txt ${outbase}snpEff_reverse_summary.genes.txt
#keep genes lost only in straight and not in reverse
comm -23 ${outbase}lost_genes_straight.txt ${outbase}lost_genes_reverse.txt > ${upper}/lost_genes.txt
comm -23 ${outbase}high_lost_genes_straight.txt ${outbase}high_lost_genes_reverse.txt > ${upper}/high_lost_genes.txt
#move to quantitative
set sample1 = `dirname $outbase | tr '/' '\n' | awk '$1~/#/{print $0}' | cut -d "#" -f 1`
set comparison = `dirname $outbase`
set comparisons = `dirname $comparison`
set patient = `dirname $comparisons`
#to quantitative
python $pipeline_installation/Breseq_pipeline/snpEff/quantitative_lost_genes.py ${comparison}/lost_genes.txt $patient/$sample1/output/${sample1}_high_non_lost_genes.txt ${comparison}/changed_genes_full.txt
#quantitative
python $pipeline_installation/Breseq_pipeline/snpEff/quantitative_lost_genes.py ${comparison}/high_lost_genes.txt $patient/$sample1/output/${sample1}_high_non_lost_genes.txt ${comparison}/high_lost_genes_full.txt
