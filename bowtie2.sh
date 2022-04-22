#!/bin/bash
  #
#  load the module
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load tophat/2.1.2
module load bowtie2/2.2.9

bowtie2 -x mouse -1 4040-KH-4.4040-KH-4_0_filtered_R1_val_1.fq -2 4040-KH-4.4040-KH-4_0_filtered_R2_val_2.fq  -S aligned_genome_sequences.sam



