#!/usr/bin/env bash

mkdir -p data_processed

FILES=(data/*/*.vcf)

for file in "${FILES[@]}"; do

fbname=$(echo $file | awk -F\. '{print $1}')

routine=$(echo $file | awk -F\/ '{print $2}')

output=$(basename $fbname)

bcftools query -f '%CHROM %POS %REF %ALT %INFO/DP [%GT]\n' $file > "data_processed/${output}_${routine}.csv"

done