#!/bin/bash

#checks if dataset is specified
if [ $# -ne 1 ]; then
    echo 'Usage: input the name of the dataset to be converted, e.g. "./merge_chromosomes dataset_x"'
    exit 1
fi

#checks if dataset exists in this folder
if [ ! -d "$1" ]; then
  echo The dataset $1 does not exist in this folder
  exit
fi

cd $1

#Filtering
for i in {1..22} 
do
/mnt/als-epilepsy/tools/plink2/plink --bfile ${i} --maf 0.01 --hwe 0.000001 --make-bed --out ${i}_filtered
done
/mnt/als-epilepsy/tools/plink2/plink --bfile X --maf 0.01 --hwe 0.000001 --make-bed --out X_filtered

#Merging datasets
echo Merging dataset $1
/mnt/als-epilepsy/tools/plink2/plink --merge-list /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.99/mergelist_chromosomes --make-bed --out $1_merged --memory 5000 --allow-no-sex
exit
