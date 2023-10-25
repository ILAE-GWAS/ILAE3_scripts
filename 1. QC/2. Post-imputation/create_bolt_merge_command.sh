
for j in {1..22}
do
echo "/mnt/als-epilepsy/tools/bin/qctool2" > merged_datasets_per_chromosome/merge_commands/chr${j}.temp
for i in bonn_illumina brussels brz_mtle-hs dublin duke_dd duke_ee epgp_humancore epicure_sp1 epicure_sp5_ger epicure_sp5_ita epipgx_omni12 epipgx_omni24a epipgx_omni24 genepa hbc hk_11 hk_2 hk_3 hk_4 hk_epipgx kora livhist melhist nbs nothen nwe phili_a5 phili_c5 phili_co phili_ao pobi popgen tcd_irish_cons ucl uk_1958; 
do echo "-g $i/filtered/${i}_${j}_filtered.bgen -s $i/${i}.sample.stub.new"; done >> merged_datasets_per_chromosome/merge_commands/chr${j}.temp
echo "-og merged_datasets_per_chromosome/chr${j}_merged.bgen -os merged_datasets_per_chromosome/all_datasets_merged.sample" >> merged_datasets_per_chromosome/merge_commands/chr${j}.temp
sed ':a;N;$!ba;s/\n/ /g' merged_datasets_per_chromosome/merge_commands/chr${j}.temp > merged_datasets_per_chromosome/merge_commands/chr${j}.sh
rm merged_datasets_per_chromosome/merge_commands/chr${j}.temp
done
