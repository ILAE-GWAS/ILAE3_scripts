# Meta-analysis weighted by standard error does not work well
# because different studies used very different transformations
# SCHEME   STDERR

# Not sure if genomic control is a good idea, given the large
# number of true associations in these three regions ...
# GENOMICCONTROL ON

# To help identify allele flips, it can be useful to track
# allele frequencies in the meta-analysis.
# AVERAGEFREQ ON
# MINMAXFREQ ON


# Meta-analysis weighted by standard error does not work well
# because different studies used very different transformations
# SCHEME   STDERR

# To help identify allele flips, it can be useful to track
# allele frequencies in the meta-analysis.
AVERAGEFREQ ON
# MINMAXFREQ ON
TRACKPOSITIONS ON
CHROMOSOMELABEL CHR
POSITIONLABEL BP
SCHEME SAMPLESIZE

MARKER   SNP
ALLELE   ALLELE1 ALLELE0
EFFECT   BETA
STDERR   SE
PVAL     P_BOLT_LMM_INF
WEIGHTLABEL	Effective_N
FREQLABEL	A1FREQ
PROCESS /mnt/ilae/work/meta_with_epi25/rsids/summary_stats/ILAE2/final_generalised_epilepsy_withX_new_phenos_het_LD_CR0.05_filtered

POSITIONLABEL POS
ADDFILTER imputationInfo > 0.3
ADDFILTER AF_Allele2 > 0.01
MARKER   SNPID
ALLELE   Allele2 Allele1
EFFECT   BETA
STDERR   SE
PVAL     p.value
WEIGHTLABEL	Effective_N
FREQLABEL	AF_Allele2
PROCESS /mnt/ilae/work/meta_with_epi25/rsids/summary_stats/Epi25/EUR/eur_rm-ibd_isGGE_cc_gwas_exl_ilae2.saige.rsid.vcf.dosage.auto_Effective_N
PROCESS /mnt/ilae/work/meta_with_epi25/rsids/summary_stats/Epi25/FIN/gender_fixed/eur_FIN-rm-ibd_isGGE_cc_gwas_exl_ilae2.saige.rsid.vcf.dosage.auto_Effective_N

#chrX
CHROMOSOMELABEL 23
PROCESS /mnt/ilae/work/meta_with_epi25/rsids/summary_stats/Epi25/chrX/eur_FIN-rm-ibd_isGGE_cc_gwas_exl_ilae2.saige.rsid.vcf.dosage.chrX2.txt
PROCESS /mnt/ilae/work/meta_with_epi25/rsids/summary_stats/Epi25/chrX/eur_rm-ibd_isGGE_cc_gwas_exl_ilae2.saige.rsid.vcf.dosage.chrX2.txt

OUTFILE ILAE3.7_Caucasian_generalised .tbl
#VERBOSE ON
MINWEIGHT 20000
ANALYZE HETEROGENEITY
