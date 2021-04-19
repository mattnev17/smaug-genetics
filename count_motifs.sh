#!/bin/bash
module load samtools

for i in `seq 1 22`; do
	for j in `seq 0 4`; do
		perl data_mgmt/data_prep/count_motifs.pl $i $j
	done
done
