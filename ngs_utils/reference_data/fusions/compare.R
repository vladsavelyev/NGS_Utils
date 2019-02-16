################ Fusions #################
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(purrr)


# simple_sv_annotation list is coming from FusionCatcher, so we get it from there to make sure it's most recent (e.g. 11 Feb, 2019 fusioncather has 7866 pairs versus 6527 in simple_sv_annotation):
# wget https://raw.githubusercontent.com/ndaniel/fusioncatcher/master/bin/generate_known.py
# grep "        \['" generate_known.py | sed "s#        \['##" | sed "s#','#,#" | sed "s#'\],##" | sed "s#'\]##" > fusioncatcher_pairs.txt
(fus_catcher = read_tsv("fusioncatcher_pairs.txt", col_names=c("pair")) %>% 
    separate(pair, c("H_gene", "T_gene"), sep = ",") %>% distinct()
)
# 7632 distinct pairs

(hmf_pairs = read_csv("knownFusionPairs.csv", quote = '"')
)

(hmf_prom_head = read_csv("knownPromiscuousFive.csv", quote = '"')
)

(hmf_prom_tail = read_csv("knownPromiscuousThree.csv", quote = '"')
)

(cancer_genes = read_tsv("umccr_cancer_genes.latest.tsv") %>% 
    filter(fusion == T) %>% 
    select(gene = symbol)
)

(pairs = hmf_pairs %>% 
  full_join(az %>% mutate(AZ = T), by = c("H_gene", "T_gene")))

# Are there genes that are both head and tail?
pairs %>% filter(H_gene %in% pairs$T_gene) %>% distinct(H_gene)
# 1318 such genes, e.g.:
pairs %>% filter(H_gene == 'ACPP' | T_gene == 'ACPP')
# and tend to be in the fusioncatcher list, which is gigantic compared to HMF. Try to remove from it promiscuous fusions?


# Do HMF pairs cover fusioncatcher?
fus_catcher %>% unite(fus, H_gene, T_gene, sep='&') %>% 
  filter(fus %in% (unite(hmf_pairs, fus, H_gene, T_gene, sep='&')$fus))
# 255 / 7632
fus_catcher %>% unite(fus, H_gene, T_gene, sep='&') %>% 
  filter(fus %in% (unite(hmf_pairs, fus, T_gene, H_gene, sep='&')$fus))
# - plus 174 if we swap T and H
# Do fusioncatcher cover HMF pairs?
hmf_pairs %>% unite(fus, H_gene, T_gene, sep='&') %>% 
  filter(fus %in% (unite(fus_catcher, fus, H_gene, T_gene, sep='&')$fus))
# 255 / 401

# Do HMF promiscous cover fusioncatcher?
fus_catcher %>% filter(H_gene %in% hmf_prom_head$gene | T_gene %in% hmf_prom_tail$gene | T_gene %in% hmf_prom_head$gene | H_gene %in% hmf_prom_tail$gene)
# 1659 / 7632
# So promuscous cover only to 14% of fusions, so we better stick to HMF fusions only.
# Also, do fusioncatcher cover HMF promiscous?
hmf_prom_head %>% filter(gene %in% fus_catcher$H_gene | gene %in% fus_catcher$T_gene)
# 29/30
hmf_prom_tail %>% filter(gene %in% fus_catcher$H_gene | gene %in% fus_catcher$T_gene)
# 36/36



# https://github.com/pmelsted/pizzly/issues/19
fus_catcher %>% filter(str_detect(H_gene, "IGH\\.*"))
fus_catcher %>% filter(str_detect(H_gene, "DUX4"))



################
# How about cancer genes?
cancer_genes
# 352
hmf_pairs %>% count(H_gene %in% cancer_genes$gene, T_gene %in% cancer_genes$gene)
#`H_gene %in% cancer_genes$gene` `T_gene %in% cancer_genes$gene`     n
# FALSE                           FALSE                              49
# FALSE                           TRUE                               72
# TRUE                            FALSE                              62
# TRUE                            TRUE                              218
# 218 full pairs, 62 heads only, 72 tails only match, and only 49 pairs do not match completely.
# Also:
hmf_prom_head %>% count(gene %in% cancer_genes$gene)  # 25/36
hmf_prom_tail %>% count(gene %in% cancer_genes$gene)  # 27/30
# Mostly matching, so we'll just add remaining fusions in the cancer gene list, and use HMF list of fusions in simple_sv_annotation.




















