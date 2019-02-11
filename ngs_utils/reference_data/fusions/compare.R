################ Fusions #################
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(purrr)



(fus_catcher = read_tsv("fusioncatcher_pairs.txt", col_names=c("pair")) %>% 
    separate(pair, c("H_gene", "T_gene"), sep = ",")
)

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
# and tend to be in the fusion catcher list, which is gigantic compared to HMF. Try to remove from it promiscuous fusions?

fus_catcher %>% 
  filter(H_gene %in% hmf_prom_head$gene | T_gene %in% hmf_prom_tail$gene)
# 1082 with "or", 46 with "and"
# So removing promuscous will still leave us with an enourmous amount of fusions. We better stick to HMF fusions only.


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




















