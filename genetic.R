# longitudinal genetic data

# pkg
library(librarian)
librarian::shelf(tidyverse,
                 ggrepel,
                 ontologyIndex,
                 data.table,
                 scales,
                 finalfit,
                 kableExtra)

# helper fn
source("func.R")

### ONTOLOGY ------------------------------------------------------------------
# load ontology
ont_hpo <- get_ontology("hp.obo.txt", 
                        propagate_relationships = "is_a", 
                        extract_tags = "everything")

# prepare UMLS-HPO map
hpo_map <- lapply(ont_hpo$xref, function(x){
  x <- x[x %like% "UMLS:"]
  x <- sub('.*\\:', '', x)
}) 

hpo_map <- enframe(hpo_map) %>%
  unnest(value) %>%
  rename(ConceptID = value)

desc_map <- tibble(term = ont_hpo$id,
                   description = ont_hpo$name)

# prepare propagation map
prop_map <- ont_hpo$ancestors %>% 
  enframe() %>%
  rename(term = name, prop_terms = value)

### DATA ----------------------------------------------------------------------
# data: all encounters per patient, ages 0-6, grouped by gene positive / negative
df_raw <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/longitudinal_genetic.csv")

# only keep columns we are interested in
df <- df_raw %>% 
  select(PatientId, # patient ID
         ConceptID, # UMLS code for encounter
         GENEPOS_comb, # binary vector: non-genetic vs likely genetic patient
         ContactAge, # relative age at encounter
         ProcAge # age at epilepsy CPT
  )

### CONSTANTS -----------------------------------------------------------------
n_patients <- length(unique(df$PatientId))
max_age <- max(df$ContactAge)

### MISSING DATA --------------------------------------------------------------
# order by PatientId
# check for NA in ContactAge
# replace NA with previous value
# then drop missing rows
df <- df %>%
  group_by(PatientId) %>%
  arrange(desc(ContactAge)) %>%
  fill(ContactAge, .direction = c("up")) %>%
  na.omit

### SUMMARY STATS -------------------------------------------------------------
# demographic table
df_person <- df_raw %>%
  group_by(PatientId, DateOfBirth, Gender, Ethnicity, GENEPOS, GENEPOS_comb, ProcAge) %>%
  summarize(max_age = max(ContactAge, na.rm = TRUE),
            min_age = min(ContactAge, na.rm = TRUE),
            median_age = median(ContactAge, na.rm = TRUE))

ls_col = c("Gender", "Ethnicity", "ProcAge", "min_age", "median_age", "max_age")

tbl1 <- df_person %>%  
  mutate(Gender = recode(Gender, 
                         "C0086582" = "Male",
                         "C0086287" = "Female")) %>%
  mutate(Ethnicity = recode(Ethnicity, 
                            "C1518424" = "Not Hispanic or Latino",
                            "C1549625" = "Unknown",
                            "C5441846" = "Hispanic or Latino",
                            "None" = "Unknown")) %>% 
  summary_factorlist("GENEPOS_comb", ls_col, p = TRUE, na_include = TRUE)  %>%
  knitr::kable("html") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
  save_kable("out/longitudinal_demographic_tbl.png", density = 900, zoom = 1.5)

# flag plot of encounters over age
df %>%
  summarise(lower = min(ContactAge), 
            upper = max(ContactAge), 
            p = mean(ContactAge)) %>%
  ggplot(aes(x = p, xmin = lower, xmax = upper, 
             y = reorder(PatientId, upper))) +
  geom_linerange(size = 0.1) +
  ylab("Individuals") +
  xlab("Age at encounter") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(0, 35), expand = FALSE) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

# density plot of encounters over time
df_filtered <- df %>%
  filter (PatientId %in% df_person$PatientId)

ggpubr::ggdensity(df_filtered, x = "ContactAge",
                  add = "mean", rug = FALSE,
                  fill = "GENEPOS_comb", palette = c("#00AFBB", "#E7B800"))

# raincloud plot of age at diagnosis (ProcAge)
devtools::source_gist("2a1bb0133ff568cbe28d", filename = "geom_flat_violin.R")

pval <- data.frame(
  group1 = "N",
  group2 = "Y",
  label = t.test(df_filtered$ProcAge ~ df_filtered$GENEPOS_comb, paired = FALSE)$p.value,
  y.position = 6
)

df_filtered %>% 
  ggplot(aes(x = GENEPOS_comb, y = ProcAge, fill = GENEPOS_comb)) +
  geom_flat_violin(position = position_nudge(x = .1, y = 0), alpha = .8) +
  guides(fill = "none", color = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Group") +
  ylab("Age at diagnosis") +
  theme_classic() +
  coord_cartesian(xlim = c(1.5, 2)) +
  geom_boxplot(width = .1, show.legend = FALSE, outlier.shape = NA, alpha = 0.5) +
  ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), 
                             label.x = 1.5, label.y = 6)

# age at encounter
ggplot(data = df_filtered, aes(x = GENEPOS_comb, y = ContactAge, fill = GENEPOS_comb)) +
  geom_flat_violin(position = position_nudge(x = .1, y = 0), alpha = .8) +
  guides(fill = "none", color = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Group") +
  ylab("Age at encounter") +
  theme_classic() + 
  coord_cartesian(xlim = c(1.5, 2)) +
  ## point cloud
  # geom_point(aes(color = GENEPOS_comb), 
  #            position = position_jitter(width = 0.15, seed = 1),
  #            size = .1, alpha = 0.01) +
  geom_boxplot(width = .1, show.legend = FALSE, outlier.shape = NA, alpha = 0.5) +
  ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), 
                             label.x = 1.5, label.y = 6)

# number of patients and years of follow-up available
df_person %>%
  ungroup() %>%
  mutate(followup_yrs = max_age-min_age) %>%
  summarize(n_patients = length(PatientId),
            followup_total = sum(followup_yrs),
            followup_min = min(followup_yrs),
            followup_max = max(followup_yrs),
            followup_mean = mean(followup_yrs),
            followup_sd = sd(followup_yrs)) %>%
  t()

# number of absolute and unique concepts
df_raw %>%
  ungroup() %>%
  summarize(n_concepts = length(ConceptID),
            n_concepts_uniq = length(unique(ConceptID))) %>%
  t()

# number of absolute and unique concepts per patient
df_raw %>%
  group_by(PatientId) %>%
  summarize(n_concepts = length(ConceptID),
            n_concepts_uniq = length(unique(ConceptID))) %>%
  ungroup() %>%
  summarize(n_concepts_mean = mean(n_concepts),
            n_concepts_min = min(n_concepts),
            n_concepts_max = max(n_concepts),
            n_concepts_sd = sd(n_concepts),
            n_concepts_unique_mean = mean(n_concepts_uniq),
            n_concepts_unique_min = min(n_concepts_uniq),
            n_concepts_unique_max = max(n_concepts_uniq),
            n_concepts_unique_sd = sd(n_concepts_uniq)) %>%
  t()

### HPO ANALYSIS --------------------------------------------------------------
# create age bins
df_binned <- df %>%
  ungroup() %>%
  mutate(bin = cut_number(ContactAge, n = 10))

# hacky way to get numeric breaks for histogram
breaks_binned <- levels(df_binned$bin) %>%
  sapply(., function(x) {gsub("\\,", " ", x)}) %>%
  parse_number()

breaks_binned <- c(breaks_binned, max(df_binned$ContactAge))

# split into age bins and print unique concepts per bin as sanity-check
ls_binned <- df_binned %>% 
  split(.$bin)

lapply(ls_binned, function(x) length(unique(x$ConceptID))) %>%
  unlist() %>%
  tibble(mean = mean(.), sd = sd(.), min = min(.), max = max(.)) %>%
  select(-.) %>%
  tail(1) %>%
  print()

# for each bin, map to HPO terms, then map to propagated HPO terms
ls_mapped <- list()
for(i in 1:length(ls_binned)){
  ls_mapped[[i]] <- left_join(ls_binned[[i]], hpo_map, by = "ConceptID") %>%
    rename(term = name) %>%
    na.omit()
  
  ls_mapped[[i]] <- left_join(ls_mapped[[i]], prop_map, by = "term")
}

# get information content (IC)
ls_terms <- ls_mapped %>%
  rbindlist() %>%
  select(PatientId, prop_terms) %>%
  .$prop_terms %>%
  unlist() %>%
  unique()

ic <- get_term_info_content(ontology = ont_hpo, term_sets = ls_terms)

ic <- tibble(term = names(ic),
             ic = ic)

# reshape each bin to unique HPO terms
# count of gene positive (Y) and negative (N) patients
# run Fisher's test
# get n most significant terms per bin
# TODO consider whether to do distinct patient-term pairs per bin
ls_pvalues <- list()
for(i in 1:length(ls_mapped)){
  df_group <- ls_mapped[[i]] %>%
    select(PatientId, GENEPOS_comb, prop_terms) %>%
    unnest(cols = c(prop_terms)) %>%
    group_by(prop_terms, GENEPOS_comb) %>%
    count(prop_terms) %>%
    pivot_wider(names_from = GENEPOS_comb, values_from = n) %>%
    replace(is.na(.), 0) %>%
    rename(term = prop_terms)
  
  # merge in description
  df_group <- left_join(df_group, desc_map, by = "term")
  
  # df_group can also be used for enrichment plots
  df_group <- df_group %>% 
    mutate(Y_out = max(df_group$Y)-Y,
           N_out = max(df_group$N)-N) %>% 
    mutate(pvalue = fish_test_it(Y, Y_out, N, N_out, "pvalue"),
           odds = fish_test_it(Y, Y_out, N, N_out, "odds"),
           freq1 = Y/max(df_group$Y),
           freq2 = N/max(df_group$N),
           color_sig = ifelse(p.adjust(pvalue, "holm") < 0.05, "<", ">"),
           size_sel = -log10(pvalue)*4)
  
  # get n best p-values for each bin for longitudinal plot
  df_group <- df_group %>%
    ungroup() %>%
    select(term, description, pvalue) %>% 
    slice_min(order_by = pvalue, n = 8, with_ties = FALSE) # non-trivial to choose
  
  # # filter: by IC
  # df_group <- df_group %>%
  #   left_join(ic, by = "term") %>%
  #   slice_max(order_by = ic, n = 2, with_ties = FALSE)
  
  # filter: by ancestors
  min_set <- ontologyIndex::minimal_set(ont_hpo, df_group$term)
  df_group <- df_group[df_group$term %in% min_set, ]
  
  # return
  ls_pvalues[[i]] <- df_group
}

# fix bin to correspond to mean age of bin for graph
seq <- seq(1, length(breaks_binned), 1)
breaks_mean <- sapply(seq, function(i) {mean(breaks_binned[i:(i+1)])}) %>%
  na.omit

# filter: for each term, keep only the bin with the highest p-value
names(ls_pvalues) <- breaks_mean

df_pvalues <- ls_pvalues %>%
  rbindlist(idcol = "bin") %>%
  mutate(bin = as.numeric(bin))

df_pvalues <- df_pvalues %>%
  group_by(term) %>% 
  slice_min(order_by = pvalue, n = 1)

# point plot: log10(pvalue) over age bins
palette <- colorRampPalette(RColorBrewer::brewer.pal(9, name = 'RdBu'))(length(breaks_mean))

df_pvalues %>%
  ggplot(aes(x = bin, y = -log10(pvalue), fill = factor(bin, levels = breaks_mean))) +
  geom_point() +
  geom_label_repel(aes(label = description), size = 3,
                   color = "black", max.overlaps = 10,
                   force_pull = 0.5) +
  scale_fill_manual(values = palette, 
                    name = "Mean age (years)",
                    breaks = breaks_mean, 
                    labels = format(round(breaks_mean, 3), nsmall = 1),
                    guide = guide_legend(override.aes = list(label = ""))) +
  theme_classic() +
  xlab("Age (years)")

# histogram: very simple way to show bin distribution
df_binned %>%
  ggplot(aes(x = ContactAge, fill = bin)) + 
  geom_histogram(breaks = breaks_binned) +
  xlab("Age at contact") +
  ylab("") +
  labs(fill = "Age bins") +
  theme_classic()

### SUBGROUPS: PREPROCESSING ---------------------------------------------------
# data: get list of MRNs per patient
df_cdkl5 <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/cdkl5_genetic.csv")
df_scn1a <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/scn1a_genetic.csv")

# preprocessing as above
df_genes <- df_raw %>%
  mutate(cdkl5 = .$MedicalRecordNumber %in% df_cdkl5$MRN) %>%
  mutate(scn1a = .$MedicalRecordNumber %in% df_scn1a$PAT_MRN_ID) %>%
  select(PatientId, # patient ID
         ConceptID, # UMLS code for encounter
         GENEPOS_comb, # binary vector: non-genetic vs likely genetic patient
         ContactAge, # relative age at encounter
         ProcAge, # age at epilepsy CPT
         cdkl5, # bool, cdkl5?
         scn1a # bool, scn1a?
  ) %>%
  group_by(PatientId) %>%
  arrange(desc(ContactAge)) %>%
  fill(ContactAge, .direction = c("up")) %>%
  na.omit

# fix group column
df_genes <- df_genes %>%
  mutate(status = 
           ifelse(scn1a == TRUE, "scn1a", 
                  ifelse(cdkl5 == TRUE, "cdkl5", 
                         ifelse(GENEPOS_comb == "N", "nongenetic",
                                ifelse(GENEPOS_comb == "Y", "genetic",
                                       NA))))) %>%
  select(PatientId, ConceptID, ContactAge, ProcAge, status)

# map to HPO and propagate
df_genes_mapped <- df_genes %>%
  left_join(hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  left_join(prop_map, by = "term") %>%
  select(-term) %>%
  rename(term = prop_terms) %>%
  na.omit %>%
  unnest(cols = c(term))

# # may also choose not to propagate to keep a clear signal
# df_genes_mapped <- df_genes %>%
#   left_join(hpo_map, by = "ConceptID") %>%
#   rename(term = name) %>%
#   na.omit

# optional: filter by ancestor
min_set <- ontologyIndex::minimal_set(ont_hpo, df_genes_mapped$term)
df_genes_mapped <- df_genes_mapped[df_genes_mapped$term %in% min_set, ]

### SUBGROUPS: ENRICHMENT PLOTS ------------------------------------------------
# note: ggrepel does not handled propagated plots well; may need manual labels
# SCN1A
res <- df_genes_mapped %>%
  filter(status %in% c("scn1a",
                       "genetic")) %>% # can choose groups here
  mutate(group = status == "scn1a") %>%
  enrichmentPlot(., ont_hpo) 

res$plot + 
  coord_fixed(xlim = c(0, .8), ylim = c(0, .8)) +
  ggtitle("SCN1A vs Genetic") +
  theme(plot.title = element_text(hjust = 0.5))

# CDKL5
res <- df_genes_mapped %>%
  filter(status %in% c("cdkl5",
                       "genetic")) %>% 
  mutate(group = status == "cdkl5") %>%
  enrichmentPlot(., ont_hpo) 

res$plot + 
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  ggtitle("CDKL5 vs Genetic") +
  theme(plot.title = element_text(hjust = 0.5))

### SUBGROUPS: LONGITUDINAL ----------------------------------------------------
# create age bins
df_genes_binned <- df_genes %>%
  ungroup() %>%
  mutate(bin = cut_number(ContactAge, n = 10))

# hacky way to get numeric breaks for histogram
breaks_binned <- levels(df_genes_binned$bin) %>%
  sapply(., function(x) {gsub("\\,", " ", x)}) %>%
  parse_number()

breaks_binned <- c(breaks_binned, max(df_genes_binned$ContactAge))

# split into age bins and print unique concepts per bin as sanity-check
ls_genes_binned <- df_genes_binned %>% 
  split(.$bin)

lapply(ls_genes_binned, function(x) length(unique(x$ConceptID))) %>%
  unlist() %>%
  tibble(mean = mean(.), sd = sd(.), min = min(.), max = max(.)) %>%
  select(-.) %>%
  tail(1) %>%
  print()

# for each bin, map to HPO terms, then map to propagated HPO terms
ls_genes_mapped <- list()
for(i in 1:length(ls_genes_binned)){
  ls_genes_mapped[[i]] <- left_join(ls_genes_binned[[i]], hpo_map, by = "ConceptID") %>%
    rename(term = name) %>%
    na.omit()
  
  ls_genes_mapped[[i]] <- left_join(ls_genes_mapped[[i]], prop_map, by = "term")
}

# # get information content (IC)
# ls_genes_terms <- ls_genes_mapped %>%
#   rbindlist() %>%
#   select(PatientId, prop_terms) %>%
#   .$prop_terms %>%
#   unlist() %>%
#   unique()
# 
# ic <- get_term_info_content(ontology = ont_hpo, term_sets = ls_genes_terms)
# 
# ic <- tibble(term = names(ic),
#              ic = ic)

# reshape each bin to unique HPO terms
# define positive and negative groups
for(i in 1:length(ls_genes_mapped)){
  ls_genes_mapped[[i]][ls_genes_mapped[[i]]$status == "cdkl5", ]$status <- "Y"
  ls_genes_mapped[[i]][ls_genes_mapped[[i]]$status != "cdkl5" &
                         ls_genes_mapped[[i]]$status != "Y", ]$status <- "N"
}

# count of each genetic subgroup
# run Fisher's test
# get n most significant terms per bin
ls_genes_pvalues <- list()
for(i in 1:length(ls_genes_mapped)){
  df_group <- ls_genes_mapped[[i]] %>%
    select(PatientId, status, prop_terms) %>%
    unnest(cols = c(prop_terms)) %>%
    group_by(prop_terms, status) %>%
    count(prop_terms) %>%
    pivot_wider(names_from = status, values_from = n) %>%
    replace(is.na(.), 0) %>%
    rename(term = prop_terms)
  
  # merge in description
  df_group <- left_join(df_group, desc_map, by = "term")
  
  # check if tbin contains any case observations
  # this is not the case for the CDKL5 cohort, due to sample size limitations
  # if true, then skip this bin
  if(is.null(df_group$Y)){next}
  
  # df_group can also be used for enrichment plots
  df_group <- df_group %>% 
    mutate(Y_out = max(df_group$Y)-Y,
           N_out = max(df_group$N)-N) %>% 
    mutate(pvalue = fish_test_it(Y, Y_out, N, N_out, "pvalue"),
           odds = fish_test_it(Y, Y_out, N, N_out, "odds"),
           freq1 = Y/max(df_group$Y),
           freq2 = N/max(df_group$N),
           color_sig = ifelse(p.adjust(pvalue, "holm") < 0.05, "<", ">"),
           size_sel = -log10(pvalue)*4)
  
  # get n best p-values for each bin for longitudinal plot
  df_group <- df_group %>%
    ungroup() %>%
    select(term, description, pvalue) %>% 
    slice_min(order_by = pvalue, n = 8, with_ties = FALSE) # non-trivial to choose
  
  # # filter: by IC
  # df_group <- df_group %>%
  #   left_join(ic, by = "term") %>%
  #   slice_max(order_by = ic, n = 2, with_ties = FALSE)
  
  # filter: by ancestors
  min_set <- ontologyIndex::minimal_set(ont_hpo, df_group$term)
  df_group <- df_group[df_group$term %in% min_set, ]
  
  # return
  ls_genes_pvalues[[i]] <- df_group
}

# fix bin to correspond to mean age of bin for graph
seq <- seq(1, length(breaks_binned), 1)
breaks_mean <- sapply(seq, function(i) {mean(breaks_binned[i:(i+1)])}) %>%
  na.omit

# filter: for each term, keep only the bin with the highest p-value
names(ls_genes_pvalues) <- breaks_mean[1:length(ls_genes_pvalues)]

df_genes_pvalues <- ls_genes_pvalues %>%
  rbindlist(idcol = "bin") %>%
  mutate(bin = as.numeric(bin))

df_genes_pvalues <- df_genes_pvalues %>%
  group_by(term) %>% 
  slice_min(order_by = pvalue, n = 1)

# point plot: log10(pvalue) over age bins
# here we can just plot genetic vs non-genetic as sanity check for the subanalysis
palette <- colorRampPalette(RColorBrewer::brewer.pal(9, name = 'RdBu'))(length(breaks_mean))

df_genes_pvalues %>%
  ggplot(aes(x = bin, y = -log10(pvalue), fill = factor(bin, levels = breaks_mean))) +
  geom_point() +
  geom_label_repel(aes(label = description), size = 3,
                   color = "black", max.overlaps = 10,
                   force_pull = 0.5) +
  scale_fill_manual(values = palette, 
                    name = "Mean age (years)",
                    breaks = breaks_mean, 
                    labels = format(round(breaks_mean, 3), nsmall = 1),
                    guide = guide_legend(override.aes = list(label = ""))) +
  theme_classic() +
  xlab("Age (years)")

