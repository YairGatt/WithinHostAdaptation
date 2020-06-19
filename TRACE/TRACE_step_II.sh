#! /bin/tcsh

#This is the second script of TRACE, it creates a tree for each strain and combined them to extract progenitor-progeny pairs

set pipeline_installation = /home/users/yair/git/WithinHostAdaptation
set kSNP_INSTALLATION = /home/users/yair/Software/kSNP3.1_Linux_package/

#directory with results of TRACE I
set strains_dir = $1
#assemblies
set assemblies_dir = $2
#millstone
set convert_file = $3
#output directory
set workdir = $4
#repository
set repository = $5
#temp dir
set temp_dir = $6
#directory with breseq results divided to directory for each isolate
set breseq_dir = $7
#core mode
set core = $8
#create temp dir
if ("$temp_dir" == "") then
  set temp_dir = /tmp/TRACE/
  mkdir -p $temp_dir
endif
#Run a tree for each strain
foreach strain (${strains_dir}/strain_*.txt)
  #create directory for results
  set strain_name = `basename $strain | sed 's/.txt//'`
  mkdir -p ${workdir}/${strain_name}
  #define possibility to continue
  set num_isolates = `cat $strain | wc -l`
  set num_timepoints = `cat $strain | rev | cut -c 2- | rev | sort -u | wc -l`
  #if there are less than two time points continue
  if ($num_timepoints < 2) then
    rm -f ${workdir}/${strain_name}/comparisons_converted.txt
    touch ${workdir}/${strain_name}/comparisons_converted.txt
    rm -f ${workdir}/${strain_name}/breseq_comparisons_converted.txt
    touch ${workdir}/${strain_name}/breseq_comparisons_converted.txt
    continue
  endif
  #define converted file
  set converted = `echo $strain | sed 's/strain_/converted_/'`
  #if there are only two isolates, there is only one possible tree, let's not waste our time
  if ($num_isolates == 2) then
    #determine first
    set first = `paste $strain $converted | sort -n | cut -f 2 | head -1`
    set second = `paste $strain $converted | sort -n | cut -f 2 | head -2 |tail -1`
    #write
    echo "$first\t$second" > ${workdir}/${strain_name}/comparisons_converted.txt
    echo "$first\t$second" > ${workdir}/${strain_name}/breseq_comparisons_converted.txt
    continue
  endif
  #if there are more than 7 samples in a single timepoint, skip
  set hs = `cat $strain | grep -c h`
  if ("$hs" > 0) then
    continue
  endif
  #create tree based on BRESEQ distances if available
  if ("$breseq_dir" != "") then
    if (-d $breseq_dir) then
      #determine isolates
      set strain_isolates = ""
      foreach isolate (`cat $converted`)
        if (-d $breseq_dir/$isolate) then
          set strain_isolates = `echo $strain_isolates $breseq_dir/${isolate}*/output/${isolate}*.gd`
        endif
      end
      #create tree
      python $pipeline_installation/TRACE/bin/step_II/create_tree_evolutionary_model.py -w ${workdir}/${strain_name}/ -s $strain -c $converted -b $strain_isolates -t 0.3
    endif
  endif
  #create tree based on kSNP matrix
  #create kSNP matrix, use existing one if already there
  #use matrix of all snps if all
  if ("$core" != "Y") then
    set matrix = ${workdir}/${strain_name}/SNPs_all_matrix.fasta
    set fasta_matrix = ${workdir}/${strain_name}/kSNP_tree_with_reference_corrected.tre
  #if only core, use only core
  else
    set matrix = ${workdir}/${strain_name}/core_SNPs_matrix.fasta
    set fasta_matrix = ${workdir}/${strain_name}/kSNP_tree_with_reference_corrected.tre
  endif
  #if the matrix doesn't exist
  if (! -e  $matrix) then
    #if no kSNP and no phylip are available, create matrix with kSNP
    python $pipeline_installation/TRACE/bin/step_II/create_tree.py $assemblies_dir $repository $temp_dir $convert_file ${workdir}/${strain_name} $kSNP_INSTALLATION/kSNP3/ $strain
  endif
  #create tree with matrix
  python $pipeline_installation/TRACE/bin/step_II/create_tree_evolutionary_model.py -w ${workdir}/${strain_name}/ -s $strain -c $converted -m $matrix -t 0.3
end
#combine all trees from all clones
#make sure there are any
set w = `ls $workdir | wc -l`
if ($w == 0) then
  continue
endif
#get all pairs from all clones
set a = `ls ${workdir}/* | grep -c comparisons_converted.txt`
set b = `ls ${workdir}/* | grep -c breseq_comparisons_converted.txt`
#combine
if ($a > 0) then
  cat ${workdir}/*/comparisons_converted.txt > ${workdir}/comparisons_converted.txt
endif
#comine breseq-based trees too if available
if ($b > 0) then
  cat ${workdir}/*/breseq_comparisons_converted.txt > ${workdir}/breseq_comparisons_converted.txt
endif
#format comparisons
python $pipeline_installation/TRACE/bin/step_II/format_comparisons.py -c $workdir/comparisons_converted.txt
