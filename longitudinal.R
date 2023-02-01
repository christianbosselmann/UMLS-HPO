# longitudinal genetic data

# pkg
library(librarian)
librarian::shelf(tidyverse,
                 ggrepel,
                 ontologyIndex,
                 data.table,
                 scales,
                 finalfit,
                 kable,
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
# TODO fix gender/ethnicity labels once we have UMLS access

ls_col = c("Gender", "Ethnicity", "ProcAge", "min_age", "median_age", "max_age")

df_person %>%  
  summary_factorlist("GENEPOS_comb", ls_col, p = TRUE, na_include = TRUE) %>%
  knitr::kable("html") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  save_kable("out/longitudinal_demographic_tbl.pdf")

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
    slice_min(order_by = pvalue, n = 10, with_ties = FALSE)
  
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

names(ls_pvalues) <- breaks_mean

df_pvalues <- ls_pvalues %>%
  rbindlist(idcol = "bin") %>%
  mutate(bin = as.numeric(bin))

# point plot: log10(pvalue) over age bins
palette <- colorRampPalette(RColorBrewer::brewer.pal(9, name = 'RdBu'))(length(breaks_mean))

df_pvalues %>%
  ggplot(aes(x = bin, y = -log10(pvalue), fill = factor(bin, levels = breaks_mean))) +
  geom_point() +
  geom_label_repel(aes(label = description), size = 3,
                   color = "black", max.overlaps = 10,
                   force_pull = 0.5) +
  scale_fill_manual(values = palette, "Mean age (years)",
                    breaks = breaks_mean, labels = format(round(breaks_mean, 3), nsmall = 1),
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

### TODO refine
