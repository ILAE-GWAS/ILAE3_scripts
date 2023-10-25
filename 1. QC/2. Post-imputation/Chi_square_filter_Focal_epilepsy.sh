for i in {1..22}
do
#/mnt/als-epilepsy/tools/plink2/plink --bfile caucasian_chr${i} --ld-window-r2 0.3 --ld-window 250 --ld-window-kb 250 --r2 --out chi2_LD_couples/chr${i}
cd /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.9/merged/stratified/caucasian/chromosomes/chi2_LD_couples

#Creating P-value, MAF and Chi-square columns
awk 'NR==FNR{pats[$1]=$13; next} $3 in pats{print pats[$3]}' ../../results/focal_epilepsy_caucasian.tab chr${i}.ld > chi1${i}
awk 'NR==FNR{pats[$1]=$13; next} $6 in pats{print pats[$6]}' ../../results/focal_epilepsy_caucasian.tab chr${i}.ld > chi2${i}
awk 'NR==FNR{pats[$1]=$14; next} $3 in pats{print pats[$3]}' ../../results/focal_epilepsy_caucasian.tab chr${i}.ld > p-value1${i}
awk 'NR==FNR{pats[$1]=$14; next} $6 in pats{print pats[$6]}' ../../results/focal_epilepsy_caucasian.tab chr${i}.ld > p-value2${i}
awk 'NR==FNR{pats[$1]=$7; next} $3 in pats{print pats[$3]}' ../../results/focal_epilepsy_caucasian.tab chr${i}.ld > MAF1${i}
awk 'NR==FNR{pats[$1]=$7; next} $6 in pats{print pats[$6]}' ../../results/focal_epilepsy_caucasian.tab chr${i}.ld > MAF2${i}

#Adding column headers
cat header_chi1 chi1${i} > chi1_with_header${i}
cat header_chi2 chi2${i} > chi2_with_header${i}
cat header_p-value1 p-value1${i} > p-value1_with_header${i}
cat header_p-value1 p-value2${i} > p-value2_with_header${i}
cat header_MAF1 MAF1${i} > MAF1_with_header${i}
cat header_MAF2 MAF2${i} > MAF2_with_header${i}
paste chr${i}.ld chi1_with_header${i} chi2_with_header${i} p-value1_with_header${i} p-value2_with_header${i} MAF1_with_header${i} MAF2_with_header${i} > chr${i}_chi1_chi2.ld

#Adding Chi-square difference column and printing final file
awk '{print $8-$9}' chr${i}_chi1_chi2.ld > chr${i}_chi1_chi2_CHI_diff.ld_temp
cat header_chi-square_diff chr${i}_chi1_chi2_CHI_diff.ld_temp > chr${i}_chi1_chi2_CHI_diff.ld
paste chr${i}_chi1_chi2.ld chr${i}_chi1_chi2_CHI_diff.ld > focal/final/chr${i}_final
rm chr${i}_chi1_chi2_CHI_diff.ld_temp


#awk '{if ($8-$9>3.57 || $9-$8>3.57) print $3,$6}' chr${i}_chi1_chi2.ld > chr${i}_exclude_temp
#sed 's/ /\n/g' chr${i}_exclude_temp > chr${i}_exclude
#rm chr${i}_exclude_temp

cd /mnt/als-epilepsy/data/imputed_2017/new_plink_info0.9/merged/stratified/caucasian/chromosomes/
done
#cat chr*_exclude > all_chr_chisquare_exclusionlist

cd chi2_LD_couples/
#Concatenate chromosomes
cat focal/final/chr*_final > focal/final/all_final

#Making LD R-square bins of 0.1
awk '{if ($7>0.3 && $7<=0.4) print $0}' focal/final/all_final > focal/final/bins/bin_0.3_0.4
awk '{if ($7>0.4 && $7<=0.5) print $0}' focal/final/all_final > focal/final/bins/bin_0.4_0.5
awk '{if ($7>0.5 && $7<=0.6) print $0}' focal/final/all_final > focal/final/bins/bin_0.5_0.6
awk '{if ($7>0.6 && $7<=0.7) print $0}' focal/final/all_final > focal/final/bins/bin_0.6_0.7
awk '{if ($7>0.7 && $7<=0.8) print $0}' focal/final/all_final > focal/final/bins/bin_0.7_0.8
awk '{if ($7>0.8 && $7<=0.9) print $0}' focal/final/all_final > focal/final/bins/bin_0.8_0.9
awk '{if ($7>0.9 && $7<=1) print $0}' focal/final/all_final > focal/final/bins/bin_0.9_1


#Calculating SD per bin
awk '{delta = $14 - avg; avg += delta / NR; mean2 += delta * ($14 - avg); } END { print sqrt(mean2 / NR); }' focal/final/bins/bin_0.3_0.4 > focal/final/bins/bin_0.3_0.4_STD
awk '{delta = $14 - avg; avg += delta / NR; mean2 += delta * ($14 - avg); } END { print sqrt(mean2 / NR); }' focal/final/bins/bin_0.4_0.5 > focal/final/bins/bin_0.4_0.5_STD
awk '{delta = $14 - avg; avg += delta / NR; mean2 += delta * ($14 - avg); } END { print sqrt(mean2 / NR); }' focal/final/bins/bin_0.5_0.6 > focal/final/bins/bin_0.5_0.6_STD
awk '{delta = $14 - avg; avg += delta / NR; mean2 += delta * ($14 - avg); } END { print sqrt(mean2 / NR); }' focal/final/bins/bin_0.6_0.7 > focal/final/bins/bin_0.6_0.7_STD
awk '{delta = $14 - avg; avg += delta / NR; mean2 += delta * ($14 - avg); } END { print sqrt(mean2 / NR); }' focal/final/bins/bin_0.7_0.8 > focal/final/bins/bin_0.7_0.8_STD
awk '{delta = $14 - avg; avg += delta / NR; mean2 += delta * ($14 - avg); } END { print sqrt(mean2 / NR); }' focal/final/bins/bin_0.8_0.9 > focal/final/bins/bin_0.8_0.9_STD
awk '{delta = $14 - avg; avg += delta / NR; mean2 += delta * ($14 - avg); } END { print sqrt(mean2 / NR); }' focal/final/bins/bin_0.9_1 > focal/final/bins/bin_0.9_1_STD


#Removing SNPs more than 4STD from mean
#awk '{if ($8-$9>5.256 || $9-$8>5.256) print $3,$6}' final/bins/bin_0.5_0.6 > exclusions/exclude_bin_0.5_0.6_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.5_0.6_temp > exclusions/exclude_bin_0.5_0.6
#awk '{if ($8-$9>5.1 || $9-$8>5.1) print $3,$6}' final/bins/bin_0.6_0.7 > exclusions/exclude_bin_0.6_0.7_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.6_0.7_temp > exclusions/exclude_bin_0.6_0.7
#awk '{if ($8-$9>4.471 || $9-$8>5.256) print $3,$6}' final/bins/bin_0.7_0.8 > exclusions/exclude_bin_0.7_0.8_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.7_0.8_temp > exclusions/exclude_bin_0.7_0.8
#awk '{if ($8-$9>4.421 || $9-$8>4.421) print $3,$6}' final/bins/bin_0.8_0.9 > exclusions/exclude_bin_0.8_0.9_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.8_0.9_temp > exclusions/exclude_bin_0.8_0.9
#awk '{if ($8-$9>3.118 || $9-$8>3.118) print $3,$6}' final/bins/bin_0.9_1 > exclusions/exclude_bin_0.9_1_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.9_1_temp > exclusions/exclude_bin_0.9_1

#Removing SNPs more than 6STD from mean
#awk '{if ($8-$9>8.2 || $9-$8>8.2) print $3,$6}' final/bins/bin_0.3_0.4 > exclusions/exclude_bin_0.3_0.4_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.3_0.4_temp > exclusions/exclude_bin_0.3_0.4
#awk '{if ($8-$9>8 || $9-$8>8) print $3,$6}' final/bins/bin_0.4_0.5 > exclusions/exclude_bin_0.4_0.5_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.4_0.5_temp > exclusions/exclude_bin_0.4_0.5
#awk '{if ($8-$9>7.884 || $9-$8>7.884) print $3,$6}' final/bins/bin_0.5_0.6 > exclusions/exclude_bin_0.5_0.6_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.5_0.6_temp > exclusions/exclude_bin_0.5_0.6
#awk '{if ($8-$9>7.6253 || $9-$8>7.6253) print $3,$6}' final/bins/bin_0.6_0.7 > exclusions/exclude_bin_0.6_0.7_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.6_0.7_temp > exclusions/exclude_bin_0.6_0.7
#awk '{if ($8-$9>7.0777 || $9-$8>7.0777) print $3,$6}' final/bins/bin_0.7_0.8 > exclusions/exclude_bin_0.7_0.8_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.7_0.8_temp > exclusions/exclude_bin_0.7_0.8
#awk '{if ($8-$9>6.63252 || $9-$8>6.63252) print $3,$6}' final/bins/bin_0.8_0.9 > exclusions/exclude_bin_0.8_0.9_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.8_0.9_temp > exclusions/exclude_bin_0.8_0.9
#awk '{if ($8-$9>4.6775 || $9-$8>4.6775) print $3,$6}' final/bins/bin_0.9_1 > exclusions/exclude_bin_0.9_1_temp
#sed 's/ /\n/g' exclusions/exclude_bin_0.9_1_temp > exclusions/exclude_bin_0.9_1


#Removing SNPs more than 6STD from mean divided by R-square squared
awk '{if (($8-$9>((3*sqrt(($8+$9)/2))/($7*$7))) || ($9-$8>((3*sqrt(($8+$9)/2))/($7*$7)))) print $3,$6}' focal/final/bins/bin_0.3_0.4 > focal/exclusions/exclude_bin_0.3_0.4_temp
sed 's/ /\n/g' focal/exclusions/exclude_bin_0.3_0.4_temp > focal/exclusions/exclude_bin_0.3_0.4
awk '{if (($8-$9>((3*sqrt(($8+$9)/2))/($7*$7))) || ($9-$8>((3*sqrt(($8+$9)/2))/($7*$7)))) print $3,$6}' focal/final/bins/bin_0.4_0.5 > focal/exclusions/exclude_bin_0.4_0.5_temp
sed 's/ /\n/g' focal/exclusions/exclude_bin_0.4_0.5_temp > focal/exclusions/exclude_bin_0.4_0.5
awk '{if (($8-$9>((3*sqrt(($8+$9)/2))/($7*$7))) || ($9-$8>((3*sqrt(($8+$9)/2))/($7*$7)))) print $3,$6}' focal/final/bins/bin_0.5_0.6 > focal/exclusions/exclude_bin_0.5_0.6_temp
sed 's/ /\n/g' focal/exclusions/exclude_bin_0.5_0.6_temp > focal/exclusions/exclude_bin_0.5_0.6
awk '{if (($8-$9>((3*sqrt(($8+$9)/2))/($7*$7))) || ($9-$8>((3*sqrt(($8+$9)/2))/($7*$7)))) print $3,$6}' focal/final/bins/bin_0.6_0.7 > focal/exclusions/exclude_bin_0.6_0.7_temp
sed 's/ /\n/g' focal/exclusions/exclude_bin_0.6_0.7_temp > focal/exclusions/exclude_bin_0.6_0.7
awk '{if (($8-$9>((3*sqrt(($8+$9)/2))/($7*$7))) || ($9-$8>((3*sqrt(($8+$9)/2))/($7*$7)))) print $3,$6}' focal/final/bins/bin_0.7_0.8 > focal/exclusions/exclude_bin_0.7_0.8_temp
sed 's/ /\n/g' focal/exclusions/exclude_bin_0.7_0.8_temp > focal/exclusions/exclude_bin_0.7_0.8
awk '{if (($8-$9>((3*sqrt(($8+$9)/2))/($7*$7))) || ($9-$8>((3*sqrt(($8+$9)/2))/($7*$7)))) print $3,$6}' focal/final/bins/bin_0.8_0.9 > focal/exclusions/exclude_bin_0.8_0.9_temp
sed 's/ /\n/g' focal/exclusions/exclude_bin_0.8_0.9_temp > focal/exclusions/exclude_bin_0.8_0.9
awk '{if (($8-$9>((3*sqrt(($8+$9)/2))/($7*$7))) || ($9-$8>((3*sqrt(($8+$9)/2))/($7*$7)))) print $3,$6}' focal/final/bins/bin_0.9_1 > focal/exclusions/exclude_bin_0.9_1_temp
sed 's/ /\n/g' focal/exclusions/exclude_bin_0.9_1_temp > focal/exclusions/exclude_bin_0.9_1

rm focal/exclusions/exclude_bin_*_temp

cat focal/exclusions/exclude_bin_* > focal/focal_all_SNPS_to_exclude_4STD_dividedby_R2_bins


#awk '{if (f==1) { r[$1] } else if (! ($1 in r)) { print $0 } } ' f=1 all_SNPS_to_exclude_4STD_dividedby_R2_bins f=2 ../../results/all_epilepsy_caucasian_heterogeneityfilter.tab > all_epilepsy_LD_filtered_4STD_dividedbyR2