python3 ../bin/annotate_CpG.py -r hg19.RefSeq.union.bed -i test_01.bed6 -o OUT1

python3 ../bin/beta_profile.py -r hg19.RefSeq.union.bed -i test_02.bed6.gz -o OUT2

python3 ../bin/chrom_distribution.py -i test_03a.bed3.gz,test_03b.bed3.gz -n 450K,850K -s hg19.chrom.sizes -o chromDist

python3 ../bin/dmc_logit.py -i test_04_TwoGroup.tsv.gz -g test_04_TwoGroup.grp.csv.gz -o OUT_4

python3 ../bin/dmc_nonparametric.py -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp.csv.gz -o OUT_05
python3 ../bin/dmc_nonparametric.py -i test_06_ThreeGroup.tsv.gz -g test_06_ThreeGroup.grp.csv.gz -o OUT_06

python3 ../bin/dmc_ttest.py -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp.csv.gz -o OUT_05
python3 ../bin/dmc_ttest.py -i test_06_ThreeGroup.tsv.gz -g test_06_ThreeGroup.grp.csv.gz -o OUT_06

python3 ../bin/dmc_fisher.py -i a_data.txt -g a_grp.csv -o OUT_10

python3 ../bin/genomic_distribution_1.py -i test_03b.bed3.gz -r hg19.RefSeq.union.bed -o OUT_7
python3 ../bin/genomic_distribution_2.py -i test_03b.bed3.gz  -b  hg19_H3K4me3.bed4,hg19_CGI.bed4,hg19_H3K27ac_with_H3K4me1.bed4,hg19_H3K27me3.bed4 -o OUT_8
 
python3 ../bin/methyl_logo.py -i test_07_450K_CH.bed -r /database/hg19.fa -o OUT_9
