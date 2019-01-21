---
title: "Cancer genes"
output: html_notebook
---
```{r set-options, echo=FALSE, cache=FALSE}
```

```{r}
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(stringi)
library(readr)
library(purrr)
library(crayon)
library(knitr)
library(kableExtra)
```

```{r}
ngc_tsgonco <-  read_tsv("/Users/vsaveliev/git/NGS_Utils/ngs_utils/reference_data/key_genes/sources/NCG6_tsgoncogene.tsv")
ngc_cancer <- read_tsv("/Users/vsaveliev/git/NGS_Utils/ngs_utils/reference_data/key_genes/sources/NCG6_cancergenes.tsv")
umccr <- read_tsv("/Users/vsaveliev/git/NGS_Utils/ngs_utils/reference_data/key_genes/sources/umccr_cancer_genes.20181112.txt", col_names = c('symbol', 'sources'))
cancermine <- read_tsv("/Users/vsaveliev/git/NGS_Utils/ngs_utils/reference_data/key_genes/sources/cancermine_collated.tsv")
az300 <- read_tsv("/Users/vsaveliev/git/NGS_Utils/ngs_utils/reference_data/key_genes/sources/az_key_genes.300.txt", col_names = c('symbol'))
oncokb <- read_tsv("/Users/vsaveliev/git/NGS_Utils/ngs_utils/reference_data/key_genes/sources/oncoKB_cancerGeneList.txt")

genes <- umccr %>% 
  separate_rows(sources, sep = "\\|") %>% 
  filter(sources != 'AZ800') %>% 
  bind_rows(az300 %>% mutate(sources = "AZ300"))

genes1 <- genes %>% 
  bind_rows(ngc_tsgonco %>% 
              mutate(sources = "NGC_known",
                     ngc_og = NCG6_oncogene == 1,
                     ngc_ts = NCG6_tsg == 1,
                     ngc_fu = str_detect(cgc_annotation, "fusion")) %>% 
              select(symbol, sources, starts_with("ngc_")))

ngc_cancer_other <- ngc_cancer %>% 
  filter(type != "Known Cancer") %>% 
  group_by(symbol) %>%
  summarise(
    ngc_n = dplyr::n(),
    type = str_c(unique(type), collapse=", "), 
    cancer_type = str_c(unique(cancer_type), collapse=", "), 
    primary_site = str_c(unique(primary_site), collapse=", "), 
    method = str_c(unique(method), collapse=", ")
  ) %>% 
  arrange(desc(ngc_n)) %>% 
  filter(ngc_n >= 3)

genes2 <- genes1 %>% 
  bind_rows(ngc_cancer_other %>% 
              mutate(sources = "NGC_other") %>% 
              select(symbol, sources, ngc_n))

cancermine_prep <- cancermine %>% 
  filter(citation_count >= 3) %>% 
  mutate(symbol = gene_normalized,
         cm_cite = citation_count,
         cm_og = role == "Oncogene",
         cm_ts = role == "Tumor_Suppressor",
         cm_driver = role == "Driver") %>% 
  select(symbol, starts_with("cm_")) %>% 
  group_by(symbol) %>% 
  summarise(sources = "Cancermine",
            cm_pub = dplyr::n(),
            cm_total_cite = sum(cm_cite),
            cm_og = mean(cm_og),
            cm_ts = mean(cm_ts),
            cm_driver = mean(cm_driver)) %>% 
  filter(cm_pub >= 2)
# cancermine_prep %>% arrange(cm_pub, cm_total_cite)

genes3 <- genes2 %>% 
  bind_rows(cancermine_prep)

# distinct(genes3, sources)

# TSG source is unrelyable - BRCA1/2 and TP53 are oncogenes :(
# genes4 <- genes3 %>% 
#   mutate(tsg_og = ifelse(sources == "Oncogenes", T, NA),
#          tsg_ts = ifelse(sources == "TumourSuppressors", T, NA))
genes4 <- genes3 %>% filter(sources != "Oncogenes" & sources != "TumourSuppressors")

# umccr
# oncokb
# oncokb_prep <- oncokb %>% 
#   dplyr::rename(
#     ok_n = `# of occurrence within resources (Column D-J)`,
#     ok_annotated = ifelse(`OncoKB Annotated` == "Yes", T, F),
#     ok_annotated = ifelse(`OncoKB Annotated` == "Yes", T, F),
#   ) %>% 
#   select(-`Entrez Gene ID`, )
#   
# genes5 <- genes4 %>%
#   bind_rows(oncokb)

# kable() %>% kable_styling()
```

```{r summarise}
genes_sum <- genes4 %>%
  group_by(symbol) %>%
  summarise(
    n = dplyr::n(),
    sources = str_c(sources, collapse = "|"),
    ngc = "NGC_known" %in% sources,
    cgc1 = "CancerGeneCensus_Tier1" %in% sources,
    cancermine = "Cancermine" %in% sources,
    # az300 = "AZ300" %in% sources,
    ngc_ts = any(ngc_ts, na.rm = T),
    cm_ts = sum(cm_ts, na.rm = T),
    ngc_og = any(ngc_og, na.rm = T),
    cm_og = sum(cm_og, na.rm = T),
    driver = sum(cm_driver, na.rm = T),
    fusion = any(ngc_fu, na.rm = T),
    # tsg_ts = any(tsg_ts, na.rm = T),
    # tsg_og = any(tsg_og, na.rm = T)
  ) %>%
  mutate(
    tumorsuppressor = ngc_ts | cm_ts > 0,
    oncogene = ngc_og | cm_og > 0
  ) %>% 
  mutate(keep = ngc | cgc1 | cancermine | n >= 2) %>% 
  select(-cm_og, -cm_ts, -ngc_ts, -ngc_og, -starts_with("ngc_"), -driver, -ngc, -cgc1, -cancermine)
```

```{r compare_to_previous_list}
added <- genes_sum %>% filter(keep, !symbol %in% umccr$symbol)
removed <- genes_sum %>% filter(symbol %in% umccr$symbol, !keep)

# # exlore sources of added genes
# genes4 %>% filter(symbol %in% added$symbol) %>% arrange(symbol)
# 
# # exlore sources of removed genes
# genes4 %>% filter(symbol %in% removed$symbol) %>% arrange(symbol)
added
removed

genes_sum %>% filter(str_detect(sources, "AZ300")) %>% select(symbol, n, keep) %>% count(keep)

umccr %>% nrow()
genes_sum %>% filter(keep) %>% nrow()
```

```{r write_results}
genes_sum %>% filter(keep) %>% select(-keep) %>% write_tsv("/Users/vsaveliev/git/NGS_Utils/ngs_utils/reference_data/key_genes/umccr_cancer_genes.20190121.tsv")
```

```{r}
genes_sum %>% 
  filter(symbol == "MYC") %>% 
  select(-tsg_ts, -tsg_og, -cm_driver)
```
```{r}
  # group_by(symbol) %>% 
  # summarise(
  #   n = dplyr::n(),
  #   cm_ts = str_c(cm_ts),
  #   cm_og = str_c(cm_og),
  #   ngc_og = ngc_og,
  #   ngc_ts = ngc_ts,
  #   ngc_fu = ngc_fu)

# Selecing:
# - Cancermine (we pre-selected at least 2 publications with at least 3 citations)
# - NGC_known
# - CancerGeneCensus_Tier1
# - At least in 2 of of clinical panels: MSKC-IMPACT, MSKC-HEME, PMCC, TS500, TEMPUS, FoundationONE, FoundationHEME
# - At least 2 sources from CancerGeneCensus_Tier2, AZ300, OncoKB-Annotated, FamilialCancer, Oncogenes, TumourSuppressors
genes4 %>% count(sources, sort = T)
# Annotating:
# - oncogene if ngc_og or mc_og >= 0.1
# - tumor_suppressor if ngc_ts or cm_ts >= 0.1
# - fusion if ngc_fu
genes4 %>% count(ngc_ts)
# The result is a list of 1387 genes:
#   - 1114 genes in at least one clinical panel
#   - 194 genes in database that are not in any of used clinical panel
#   - 79 genes in ensemble (2+) lists that are not included above

# | any(cm_og, na.rm = T) | any(cm_og, na.rm = T))
```

```{r}
genes_sum %>% count(cm_ts > 0, ngc_ts, tsg_ts)
genes_sum %>% filter(cm_ts > 0, ngc_ts, tsg_ts)
genes_sum %>% count(cm_og > 0, ngc_og, tsg_og)
genes_sum %>% filter(symbol == 'BRCA1')
```

```{r}
cancermine_oncogene <- cancermine %>% filter(role == "Oncogene")
ncg_oncogene <- ngc_tsgonco %>% filter(NCG6_oncogene == 1)
cancermine_oncogene
ncg_oncogene
```

```{r}
ngc_cancer$symbol %>% unique() %>% length()
```

```{r}
ngc_tsgonco$symbol %>% unique() %>% length()
```

```{r}
dplyr::intersect(umccr, unique(ngc_cancer$symbol)) %>% unique() %>% length()
```
```{r}
dplyr::intersect(umccr, ngc_tsgonco$symbol %>% unique()) %>% unique() %>% length()
```
```{r}
dplyr::setdiff(unique(ngc_tsgonco$symbol), umccr) %>% unique()
```

```{r}
ngc_cancer_collapsed <- ngc_cancer %>% 
  group_by(symbol) %>%
  summarise(
    n = dplyr::n(),
    type = str_c(unique(type), collapse=", "), 
    cancer_type = str_c(unique(cancer_type), collapse=", "), 
    primary_site = str_c(unique(primary_site), collapse=", "), 
    method = str_c(unique(method), collapse=", ")
  )

ngc_cancer_collapsed %>% filter(symbol %in% umccr) %>% arrange(desc(n)) %>% filter(n <= 8)
```


```{r}
ngc_cancer %>% filter(symbol == "MKL1")
```
