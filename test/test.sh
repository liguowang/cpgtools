python3 ../bin/annotate_CpG.py -r hg19.RefSeq.union.bed -i test_01.bed6 -o OUT1

python3 ../bin/beta_profile.py -r hg19.RefSeq.union.bed -i test_02.bed6.gz -o OUT2

python3 ../bin/chrom_distribution.py -i test_03a.bed3.gz,test_03b.bed3.gz -n 450K,850K -s hg19.chrom.sizes -o chromDist

python3 ../bin/dmc_nonparametric.py -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp.csv.gz -o OUT_05
python3 ../bin/dmc_nonparametric.py -i test_06_ThreeGroup.tsv.gz -g test_06_ThreeGroup.grp.csv.gz -o OUT_06

python3 ../bin/dmc_ttest.py -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp.csv.gz -o OUT_05
python3 ../bin/dmc_ttest.py -i test_06_ThreeGroup.tsv.gz -g test_06_ThreeGroup.grp.csv.gz -o OUT_06

python3 ../bin/dmc_fisher.py -i a_data.txt -g a_grp.csv -o OUT_10

python3 ../bin/genomic_distribution_1.py -i test_03b.bed3.gz -r hg19.RefSeq.union.bed -o OUT_7
python3 ../bin/genomic_distribution_2.py -i test_03b.bed3.gz  -b  hg19_H3K4me3.bed4,hg19_CGI.bed4,hg19_H3K27ac_with_H3K4me1.bed4,hg19_H3K27me3.bed4 -o OUT_8
 
python3 ../bin/methyl_logo.py -i test_07_450K_CH.bed -r /database/hg19.fa -o OUT_9

python3 ../bin/dmc_logit.py -i test_04_TwoGroup.tsv -g test_04_TwoGroup.grp.csv -o OUT_4_logit_quasi
python3 ../bin/dmc_logit.py -i test_04_TwoGroup.tsv -g test_04_TwoGroup.grp.csv -o OUT_4_logit -f 2 
python3 ../bin/dmc_bb.py -i test_04_TwoGroup.tsv -g test_04_TwoGroup.grp.csv -o OUT_4_bb
python3 ../bin/dmc_glm.py  -i test_05_TwoGroup.tsv -g test_05_TwoGroup.grp.csv -o OUT_5_glm


python3 ../bin/beta_jitterPlot.py -f 1 -i test_05_TwoGroup.tsv
python3 ../bin/beta_PCA.py -i test_12.dataFrame.tsv -g test_12.group.csv -n 3 -o aa
