#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bowtie2/2.2.9

bowtie2-build -f GRCm39.genome.fa mouse
