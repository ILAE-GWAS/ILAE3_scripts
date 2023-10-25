#!/bin/bash

find * -prune -type d | while IFS= read -r d; do 
/mnt/als-epilepsy/tools/plink2/plink --bfile /mnt/als-epilepsy/data/postqc_plink/new_plink_post_QC/${d}/${d}-updated --extract /mnt/als-epilepsy/data/postqc_plink/new_plink_post_QC/homo/${d}.extract --exclude /mnt/als-epilepsy/data/postqc_plink/new_plink_post_QC/all_bad_snps_homogeneity_filter_ita_hetfix_rsids --a2-allele /mnt/als-epilepsy/data/postqc_plink/new_plink_post_QC/${d}/Force-Allele1-${d}-HRC.txt --make-bed --allow-no-sex --output-chr M --out final_post_qc_plink_files/${d}_temp
/mnt/als-epilepsy/tools/plink2/plink --bfile /mnt/als-epilepsy/data/postqc_plink/new_plink_post_QC/final_post_qc_plink_files/${d}_temp --bmerge /mnt/als-epilepsy/data/postqc_plink/new_plink_post_QC/${d}/${d}-updated-chr23 --a2-allele /mnt/als-epilepsy/data/postqc_plink/new_plink_post_QC/${d}/Force-Allele1-${d}-HRC.txt --make-bed --allow-no-sex --output-chr M --out final_post_qc_plink_files/${d}
rm final_post_qc_plink_files/${d}_temp*
/mnt/als-epilepsy/tools/plink2/plink --bfile /mnt/als-epilepsy/data/postqc_plink/new_plink_post_QC/final_post_qc_plink_files/$d --recode vcf bgz --a2-allele /mnt/als-epilepsy/data/postqc_plink/new_plink_post_QC/${d}/Force-Allele1-${d}-HRC.txt --output-chr M --allow-no-sex --out /mnt/als-epilepsy/data/postqc_plink/new_plink_post_QC/final_vcf_files/$d
done
