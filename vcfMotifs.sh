#!/bin/bash
module load bcftools
module load samtools
module load vcftools
perl data_mgmt/data_prep/vcf_to_summary.pl copy all
