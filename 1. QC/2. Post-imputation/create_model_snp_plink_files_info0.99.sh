#Caucasians
for i in {1..22}
do
/mnt/als-epilepsy/tools/plink2/plink --bfile chr${i}_merged --keep-allele-order --maf 0.01 --geno 0.01 --hwe 0.000001 --extract range /mnt/als-epilepsy/data/imputed_2017/bgen_files/snpstats/info0.99/caucasian/chr${i}_info0.99 --make-bed --out modelsnps/caucasian/chromosomes/caucasian_chr${i}_modelsnps --keep /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.99/merged_datasets/stratified/caucasian.subjects --remove /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.99/merged_datasets/duplicates_and_outliers_to_exclude.txt
done

for f in modelsnps/caucasian/chromosomes/caucasian_chr*_modelsnps.fam;
do FILENAME=${f%%.*};
echo ${FILENAME};
done > modelsnps/caucasian/chromosomes/merge_list

/mnt/als-epilepsy/tools/plink2/plink --merge-list modelsnps/caucasian/chromosomes/merge_list --keep-allele-order --make-bed --out modelsnps/caucasian/caucasian_modelsnps_info0.99
cd modelsnps/caucasian/
/mnt/als-epilepsy/tools/ldak5.linux --bfile caucasian_modelsnps_info0.99 --cut-genes highld --genefile highld.txt --ignore-weights YES
/mnt/als-epilepsy/tools/ldak5.linux --thin thin --window-kb 1000 --window-prune .2 --bfile caucasian_modelsnps_info0.99 --exclude highld/genes.predictors.used 
nohup /mnt/als-epilepsy/tools/plink2/plink --bfile caucasian_modelsnps_info0.99 --extract thin.in --pca --out caucasian_modelsnps_info0.99_PCA &
cd /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.9/merged


#Asians
for i in {1..22}
do
/mnt/als-epilepsy/tools/plink2/plink --bfile chr${i}_merged --keep-allele-order --maf 0.01 --geno 0.01 --hwe 0.000001 --extract range /mnt/als-epilepsy/data/imputed_2017/bgen_files/snpstats/info0.99/asian/chr${i}_info0.99 --make-bed --out modelsnps/asian/chromosomes/asian_chr${i}_modelsnps --keep /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.99/merged_datasets/stratified/asian.subjects --remove /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.99/merged_datasets/duplicates_and_outliers_to_exclude.txt
done

for f in modelsnps/asian/chromosomes/asian_chr*_modelsnps.fam;
do FILENAME=${f%%.*};
echo ${FILENAME};
done > modelsnps/asian/chromosomes/merge_list

/mnt/als-epilepsy/tools/plink2/plink --merge-list modelsnps/asian/chromosomes/merge_list --keep-allele-order --make-bed --out modelsnps/asian/asian_modelsnps_info0.99
cd modelsnps/asian/
/mnt/als-epilepsy/tools/ldak5.linux --bfile asian_modelsnps_info0.99 --cut-genes highld --genefile highld.txt --ignore-weights YES
/mnt/als-epilepsy/tools/ldak5.linux --thin thin --window-kb 1000 --window-prune .2 --bfile asian_modelsnps_info0.99 --exclude highld/genes.predictors.used 
cd /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.9/merged



#African-Americans
for i in {1..22}
do
/mnt/als-epilepsy/tools/plink2/plink --bfile chr${i}_merged --keep-allele-order --maf 0.01 --geno 0.01 --hwe 0.000001 --extract range /mnt/als-epilepsy/data/imputed_2017/bgen_files/snpstats/info0.99/african/chr${i}_info0.99 --make-bed --out modelsnps/african/chromosomes/african_chr${i}_modelsnps --keep /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.99/merged_datasets/stratified/african-american.subjects --remove /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.99/merged_datasets/duplicates_and_outliers_to_exclude.txt
done

for f in modelsnps/african/chromosomes/african_chr*_modelsnps.fam;
do FILENAME=${f%%.*};
echo ${FILENAME};
done > modelsnps/african/chromosomes/merge_list

/mnt/als-epilepsy/tools/plink2/plink --merge-list modelsnps/african/chromosomes/merge_list --keep-allele-order --make-bed --out modelsnps/african/african_modelsnps_info0.99
cd modelsnps/african/
/mnt/als-epilepsy/tools/ldak5.linux --bfile african_modelsnps_info0.99 --cut-genes highld --genefile highld.txt --ignore-weights YES
/mnt/als-epilepsy/tools/ldak5.linux --thin thin --window-kb 1000 --window-prune .2 --bfile african_modelsnps_info0.99 --exclude highld/genes.predictors.used 
cd /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.9/merged