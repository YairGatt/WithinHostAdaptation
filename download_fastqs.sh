#! /bin/tcsh

set pipeline_installation = /home/users/yair/git/WithinHostAdaptation/
set SRATOOLS_PATH = /home/users/yair/Software/sratoolkit/bin/

$SRATOOLS_PATH/fastq-dump ERR212989 --split-files --gzip -O $pipeline_installation/Sample_data/Staphylococcus_aureus/trial_122_patient_4/fastqs/
$SRATOOLS_PATH/fastq-dump ERR715447 --split-files --gzip -O $pipeline_installation/Sample_data/Staphylococcus_aureus/trial_122_patient_4/fastqs/
