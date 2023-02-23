### UMLS-HPO ANALYSIS OF INDIVIDUALS WITH GENETIC EPILEPSY SYNDROMES -----------
##
## Author: Christian Bosselmann, MD
##
## Date Created: 2023-02-07
##
## Copyright (c) Christian Bosselmann, 2023
## Email: bosselc@ccf.org
##
### ----------------------------------------------------------------------------

### HEADER ---------------------------------------------------------------------
# packages
library(librarian)
librarian::shelf(tidyverse,
                 ggrepel,
                 ontologyIndex,
                 data.table,
                 scales,
                 finalfit,
                 kableExtra,
                 MatchIt,
                 kernlab,
                 Spectrum,
                 CGEN,
                 egg,
                 scales,
                 qgraph,
                 igraph,
                 bnlearn,
                 ggplotify)

# functions
source("func.R")

# seed
set.seed(42)

### PARAMETERS ----------------------------------------------------------------
# matching method for matchit
flag_match <- "nearest"

### ONTOLOGY -------------------------------------------------------------------
# load ontologyIndex object
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

# prepare HPO term description map
desc_map <- tibble(term = ont_hpo$id,
                   description = ont_hpo$name)

# prepare propagation map
prop_map <- ont_hpo$ancestors %>% 
  enframe() %>%
  rename(term = name, prop_terms = value)

### DATA ----------------------------------------------------------------------
# data: all encounters per patient, ages 0-6, grouped by gene positive / negative
df_raw <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/longitudinal_genetic.csv")

# data: list of MRNs per patient
df_cdkl5 <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/cdkl5_genetic.csv")
df_scn1a <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/scn1a_genetic.csv")
df_tsc <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/tsc_genetic.csv")

# preprocessing: all
df <- df_raw %>% 
  select(PatientId, # patient ID
         ConceptID, # UMLS code for encounter
         GENEPOS_comb, # binary vector: non-genetic vs likely genetic patient
         ContactAge, # relative age at encounter
         ProcAge # age at epilepsy CPT
  ) %>%
  group_by(PatientId) %>%
  arrange(desc(ContactAge)) %>%
  # missing data imputation
  fill(ContactAge, .direction = c("up")) %>%
  na.omit

# preprocessing: gene subgroups
df_genes <- df_raw %>%
  # mutate in logical flag for gene subgroup membership
  mutate(cdkl5 = .$MedicalRecordNumber %in% df_cdkl5$MRN) %>%
  mutate(scn1a = .$MedicalRecordNumber %in% df_scn1a$PAT_MRN_ID) %>%
  mutate(tsc = .$MedicalRecordNumber %in% df_tsc$MedicalRecordNumber) %>%
  select(PatientId, ConceptID, GENEPOS_comb, ContactAge, ProcAge, cdkl5, scn1a, tsc) %>%
  group_by(PatientId) %>%
  arrange(desc(ContactAge)) %>%
  # missing data imputation
  fill(ContactAge, .direction = c("up")) %>%
  na.omit %>%
  # merging group columns, descriptive labels
  mutate(status = 
           ifelse(scn1a == TRUE, "scn1a", 
                  ifelse(cdkl5 == TRUE, "cdkl5", 
                         ifelse(tsc == TRUE, "tsc",
                                ifelse(GENEPOS_comb == "N", "nongenetic",
                                       ifelse(GENEPOS_comb == "Y", "genetic",
                                              NA)))))) %>%
  select(PatientId, ConceptID, ContactAge, ProcAge, status)

# by-patient demographic data and age statistics
df_person <- df_raw %>%
  group_by(PatientId, DateOfBirth, Gender, Ethnicity, GENEPOS, GENEPOS_comb, ProcAge) %>%
  summarize(max_age = max(ContactAge, na.rm = TRUE),
            min_age = min(ContactAge, na.rm = TRUE),
            median_age = median(ContactAge, na.rm = TRUE))

# map subgroup dataframe to HPO and propagate
df_genes_mapped <- df_genes %>%
  left_join(hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  left_join(prop_map, by = "term") %>%
  select(-term) %>%
  rename(term = prop_terms) %>%
  na.omit %>%
  unnest(cols = c(term))

# get term information content
vec_ic <- split(df_genes_mapped, df_genes_mapped$PatientId, df_genes_mapped$term) %>%
  lapply(function(x){x <- x$term}) %>%
  get_term_info_content(ont_hpo, ., patch_missing = FALSE) 

df_ic <- data.frame(term = names(vec_ic), ic = vec_ic)

# define cohorts and match by age, sex and ethnicity; maintain label
df_match <- df_genes %>%
  left_join(df_person[,c("PatientId", "Gender", "Ethnicity", "median_age")], by = "PatientId") %>%
  distinct(PatientId, Gender, Ethnicity, median_age, status)

# set covariates as factors
df_match <- df_match %>%
  mutate(Gender = as.factor(Gender)) %>%
  mutate(Ethnicity = as.factor(Ethnicity))

## Group 0: non-genetic vs. non-genetic (null)
df_match0 <- df_match %>%
  filter(status == "nongenetic") %>%
  mutate(status = sample(0:1, n(), replace = TRUE))

df_match0 <- matchit(status ~ median_age + Ethnicity + Gender, 
                     data = df_match0, ratio = 1,
                     method = flag_match, distance = "glm")

df_match0 <- match.data(df_match0)

# optional: downsample number of encounters per patient to control for group diff.
df_match0 <- downsampleMatch(df_match0, df)

df_match0 <- df_match0 %>%
  # # merge in ConceptIDs; not done after downsampling
  # left_join(df[ ,c("PatientId", "ConceptID")], by = "PatientId") %>%
  # merge in propagated HPO terms
  left_join(hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  left_join(prop_map, by = "term") %>%
  select(-term) %>%
  rename(term = prop_terms) %>%
  na.omit %>%
  unnest(cols = c(term)) %>%
  # recode for later enrichment plots
  rename(group = status) %>%
  mutate(group = as.logical(group))

## Group 1: genetic vs. non-genetic
df_match1 <- df_match %>%
  mutate(status = recode(status, 
                         "nongenetic" = 0,
                         "genetic" = 1,
                         "scn1a" = 1,
                         "cdkl5" = 1,
                         "tsc" = 1)) 

df_match1 <- matchit(status ~ median_age + Ethnicity + Gender, 
                     data = df_match1, ratio = 1,
                     method = flag_match, distance = "glm")

df_match1 <- match.data(df_match1)

# optional: downsample number of encounters per patient to control for group diff.
df_match1 <- downsampleMatch(df_match1, df)

df_match1 <- df_match1 %>%
  # # merge in ConceptIDs; not done after downsampling
  # left_join(df[ ,c("PatientId", "ConceptID")], by = "PatientId") %>%
  # merge in propagated HPO terms
  left_join(hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  left_join(prop_map, by = "term") %>%
  select(-term) %>%
  rename(term = prop_terms) %>%
  na.omit %>%
  unnest(cols = c(term)) %>%
  # recode for later enrichment plots
  rename(group = status) %>%
  mutate(group = as.logical(group))

## Group 2: SCN1A vs. Genetic
df_match2 <- df_match %>%
  filter(status %in% c("genetic", "scn1a", "cdkl5", "tsc")) %>%
  mutate(status = recode(status, 
                         "scn1a" = 1,
                         "genetic" = 0,
                         "cdkl5" = 0,
                         "tsc" = 0)) 

df_match2 <- matchit(status ~ median_age + Ethnicity + Gender, 
                     data = df_match2, ratio = 1,
                     method = flag_match, distance = "glm")

df_match2 <- match.data(df_match2)

df_match2 <- downsampleMatch(df_match2, df)

df_match2 <- df_match2 %>%
  # # merge in ConceptIDs
  # left_join(df[ ,c("PatientId", "ConceptID")], by = "PatientId") %>%
  # merge in propagated HPO terms
  left_join(hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  left_join(prop_map, by = "term") %>%
  select(-term) %>%
  rename(term = prop_terms) %>%
  na.omit %>%
  unnest(cols = c(term)) %>%
  # recode for later enrichment plots
  rename(group = status) %>%
  mutate(group = as.logical(group))

## Group 3: CDKL5 vs. Genetic
df_match3 <- df_match %>%
  filter(status %in% c("genetic", "scn1a", "cdkl5", "tsc")) %>%
  mutate(status = recode(status, 
                         "cdkl5" = 1,
                         "genetic" = 0,
                         "scn1a" = 0,
                         "tsc" = 0)) 

df_match3 <- matchit(status ~ median_age + Ethnicity + Gender, 
                     data = df_match3, ratio = 1,
                     method = flag_match, distance = "glm")

df_match3 <- match.data(df_match3)

df_match3 <- downsampleMatch(df_match3, df)

df_match3 <- df_match3 %>%
  # # merge in ConceptIDs
  # left_join(df[ ,c("PatientId", "ConceptID")], by = "PatientId") %>%
  # merge in propagated HPO terms
  left_join(hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  left_join(prop_map, by = "term") %>%
  select(-term) %>%
  rename(term = prop_terms) %>%
  na.omit %>%
  unnest(cols = c(term)) %>%
  # recode for later enrichment plots
  rename(group = status) %>%
  mutate(group = as.logical(group))

## Group 4: SCN1A vs. Non-Genetic
df_match4 <- df_match %>%
  filter(status %in% c("nongenetic", "scn1a")) %>%
  mutate(status = recode(status, 
                         "nongenetic" = 0,
                         "scn1a" = 1)) 

df_match4 <- matchit(status ~ median_age + Ethnicity + Gender, 
                     data = df_match4, ratio = 1,
                     method = flag_match, distance = "glm")

df_match4 <- match.data(df_match4)

df_match4 <- downsampleMatch(df_match4, df)

df_match4 <- df_match4 %>%
  # # merge in ConceptIDs
  # left_join(df[ ,c("PatientId", "ConceptID")], by = "PatientId") %>%
  # merge in propagated HPO terms
  left_join(hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  left_join(prop_map, by = "term") %>%
  select(-term) %>%
  rename(term = prop_terms) %>%
  na.omit %>%
  unnest(cols = c(term)) %>%
  # recode for later enrichment plots
  rename(group = status) %>%
  mutate(group = as.logical(group))

## Group 5: CDKL5 vs. Non-Genetic
df_match5 <- df_match %>%
  filter(status %in% c("nongenetic", "cdkl5")) %>%
  mutate(status = recode(status, 
                         "nongenetic" = 0,
                         "cdkl5" = 1)) 

df_match5 <- matchit(status ~ median_age + Ethnicity + Gender, 
                     data = df_match5, ratio = 1,
                     method = flag_match, distance = "glm")

df_match5 <- match.data(df_match5)

df_match5 <- downsampleMatch(df_match5, df)

df_match5 <- df_match5 %>%
  # # merge in ConceptIDs
  # left_join(df[ ,c("PatientId", "ConceptID")], by = "PatientId") %>%
  # merge in propagated HPO terms
  left_join(hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  left_join(prop_map, by = "term") %>%
  select(-term) %>%
  rename(term = prop_terms) %>%
  na.omit %>%
  unnest(cols = c(term)) %>%
  # recode for later enrichment plots
  rename(group = status) %>%
  mutate(group = as.logical(group))

## Group 6: TSC vs. Genetic
df_match6 <- df_match %>%
  filter(status %in% c("genetic", "scn1a", "cdkl5", "tsc")) %>%
  mutate(status = recode(status, 
                         "tsc" = 1,
                         "genetic" = 0,
                         "scn1a" = 0,
                         "cdkl5" = 0)) 

df_match6 <- matchit(status ~ median_age + Ethnicity + Gender, 
                     data = df_match6, ratio = 1,
                     method = flag_match, distance = "glm")

df_match6 <- match.data(df_match6)

# optional: downsample number of encounters per patient to control for group diff.
df_match6 <- downsampleMatch(df_match6, df)

df_match6 <- df_match6 %>%
  # # merge in ConceptIDs; not done after downsampling
  # left_join(df[ ,c("PatientId", "ConceptID")], by = "PatientId") %>%
  # merge in propagated HPO terms
  left_join(hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  left_join(prop_map, by = "term") %>%
  select(-term) %>%
  rename(term = prop_terms) %>%
  na.omit %>%
  unnest(cols = c(term)) %>%
  # recode for later enrichment plots
  rename(group = status) %>%
  mutate(group = as.logical(group))

## Group 7: TSC vs. Non-Genetic
df_match7 <- df_match %>%
  filter(status %in% c("nongenetic", "tsc")) %>%
  mutate(status = recode(status, 
                         "nongenetic" = 0,
                         "tsc" = 1)) 

df_match7 <- matchit(status ~ median_age + Ethnicity + Gender, 
                     data = df_match7, ratio = 1,
                     method = flag_match, distance = "glm")

df_match7 <- match.data(df_match7)

df_match7 <- downsampleMatch(df_match7, df)

df_match7 <- df_match7 %>%
  # # merge in ConceptIDs
  # left_join(df[ ,c("PatientId", "ConceptID")], by = "PatientId") %>%
  # merge in propagated HPO terms
  left_join(hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  left_join(prop_map, by = "term") %>%
  select(-term) %>%
  rename(term = prop_terms) %>%
  na.omit %>%
  unnest(cols = c(term)) %>%
  # recode for later enrichment plots
  rename(group = status) %>%
  mutate(group = as.logical(group))

## Group 8: SCN1A vs. CDKL5
df_match8 <- df_match %>%
  filter(status %in% c("scn1a", "cdkl5")) %>%
  mutate(status = recode(status, 
                         "cdkl5" = 0,
                         "scn1a" = 1)) 

# df_match8 <- matchit(status ~ median_age + Ethnicity + Gender, 
#                      data = df_match8, ratio = 1,
#                      method = flag_match, distance = "glm")
# 
# df_match8 <- match.data(df_match8)
# 
# df_match8 <- downsampleMatch(df_match8, df)

df_match8 <- df_match8 %>%
  # # merge in ConceptIDs
  left_join(df[ ,c("PatientId", "ConceptID")], by = "PatientId") %>%
  # merge in propagated HPO terms
  left_join(hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  left_join(prop_map, by = "term") %>%
  select(-term) %>%
  rename(term = prop_terms) %>%
  na.omit %>%
  unnest(cols = c(term)) %>%
  # recode for later enrichment plots
  rename(group = status) %>%
  mutate(group = as.logical(group))

### SUMMARY STATS -------------------------------------------------------------
## demographic table
tbl_person <- df_person %>%  
  # recode UMLS to sex
  mutate(Gender = recode(Gender, 
                         "C0086582" = "Male",
                         "C0086287" = "Female")) %>%
  # recode UMLS to ethnicity
  mutate(Ethnicity = recode(Ethnicity, 
                            "C1518424" = "Not Hispanic or Latino",
                            "C1549625" = "Unknown",
                            "C5441846" = "Hispanic or Latino",
                            "None" = "Unknown")) %>% 
  summary_factorlist(dependent = "GENEPOS_comb", 
                     explanatory = c("Gender", "Ethnicity", "ProcAge", "min_age", "median_age", "max_age"),
                     p = TRUE, na_include = TRUE)  %>%
  knitr::kable("html") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) 

## p1: flag plot of encounters over age
p1 <- df %>%
  summarise(lower = min(ContactAge), 
            upper = max(ContactAge), 
            p = mean(ContactAge)) %>%
  ggplot(aes(x = p, xmin = lower, xmax = upper, 
             y = reorder(PatientId, upper))) +
  geom_linerange(size = 0.1) +
  ylab("Individuals") +
  xlab("Age at encounter (years)") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(0, 35), expand = FALSE) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

# stats of length of follow-up
stats_followup <- p1$data %>%
  mutate(dur = upper-lower) %>%
  summarize(mean = mean(dur), median = median(dur),
            sd = sd(dur), min = min(dur), max = max(dur),
            iqr = IQR(dur))

# add mean age at follow-up back to plot
p1 <- p1 + 
  geom_vline(xintercept = stats_followup$mean, linetype = "dashed") +
  geom_text(aes(x = stats_followup$mean, label = "\nMean: 6.5 years", y = 400),
            colour = "black", angle = 90, size = 4)

## p2: flat violin (raincloud) plot of age at encounter
p2 <- ggplot(data = df, aes(x = GENEPOS_comb, y = ContactAge, fill = GENEPOS_comb)) +
  geom_flat_violin(position = position_nudge(x = .1, y = 0), alpha = .8) +
  guides(fill = "none", color = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Group") +
  ylab("Age at encounter") +
  theme_classic() + 
  coord_cartesian(xlim = c(1.5, 2)) +
  geom_boxplot(width = .1, show.legend = FALSE, outlier.shape = NA, alpha = 0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..), 
                             comparisons = list(c("N", "Y")),
                             label.x = 1.5, label.y = c(20))

## p3: flat violin (raincloud) plot of age at diagnosis
p3 <- ggplot(data = df, aes(x = GENEPOS_comb, y = ProcAge, fill = GENEPOS_comb)) +
  geom_flat_violin(position = position_nudge(x = .1, y = 0), alpha = .8) +
  guides(fill = "none", color = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Group") +
  ylab("Age at diagnosis") +
  theme_classic() + 
  coord_cartesian(xlim = c(1.5, 2), ylim = c(0, 7)) +
  geom_boxplot(width = .1, show.legend = FALSE, outlier.shape = NA, alpha = 0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             comparisons = list(c("N", "Y")),
                             label.x = 1.5, label.y = c(6))

## table of relative encounter frequency
df_encounters <- df %>%
  group_by(PatientId) %>%
  mutate(encounter_freq = n_distinct(ContactAge)) %>%
  distinct(PatientId, GENEPOS_comb, encounter_freq)

# t-test: mean difference of number of encounters between groups
p_encounters <- df_encounters %>%
  mutate(GENEPOS_comb = as.factor(GENEPOS_comb)) %>%
  ungroup() %>%
  summarize(pval = t.test(encounter_freq ~ GENEPOS_comb)$p.value)

# summary statistics of encounter frequency
stats_encounters <- df_encounters %>%
  group_by(GENEPOS_comb) %>%
  summarize(mean = mean(encounter_freq), median = median(encounter_freq),
            sd = sd(encounter_freq),
            min = min(encounter_freq), max = max(encounter_freq),
            n = n_distinct(PatientId))

# plot encounter frequency over time
df_encounters_freq <- df %>%
  group_by(PatientId) %>%
  select(-ConceptID) %>%
  distinct(PatientId, GENEPOS_comb, ContactAge) %>%
  count(cut_width(ContactAge, width = 1, boundary = 0, labels = F)) %>%
  ungroup() %>%
  # maintain group label
  left_join(df_encounters[ ,c("GENEPOS_comb", "PatientId")], by = "PatientId") %>%
  # fix bin label
  rename(bin = `cut_width(ContactAge, width = 1, boundary = 0, labels = F)`) %>%
  group_by(bin, GENEPOS_comb) %>%
  # get mean number of encounters per bin per group
  summarize(mean = mean(n), sd = sd(n))

# get pvalue for encounter freq during transition (age range 18-20 years)
stats_encounters_freq <- df %>%
  group_by(PatientId) %>%
  select(-ConceptID) %>%
  distinct(PatientId, GENEPOS_comb, ContactAge) %>%
  count(cut_width(ContactAge, width = 1, boundary = 0, labels = F)) %>%
  ungroup() %>%
  # maintain group label
  left_join(df_encounters[ ,c("GENEPOS_comb", "PatientId")], by = "PatientId") %>%
  # fix bin label
  rename(bin = `cut_width(ContactAge, width = 1, boundary = 0, labels = F)`) %>%
  # age filter for transition
  filter(bin %in% c(18, 19, 20)) %>%
  ungroup() %>%
  summarize(pval = t.test(n ~ GENEPOS_comb)$p.value)

p4 <- df_encounters_freq %>%
  mutate(ymin = mean-sd, ymax = mean+sd) %>%
  ggplot(aes(x = bin-1, y = mean, color = GENEPOS_comb, fill = GENEPOS_comb)) +
  geom_smooth(se = TRUE) +
  # geom_ribbon(aes(ymin = ymin, ymax = ymax, alpha = 0.1)) +
  theme_classic() +
  guides(fill = guide_legend(title = "Group"),
         color = guide_legend(title = "Group")) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(0, 25), expand = FALSE) +
  ylab("Mean encounters per year") +
  xlab("Age (years)") +
  # add horizontal lines to define transition age for next panel
  geom_vline(xintercept = 18, linetype = "dashed") +
  geom_vline(xintercept = 20, linetype = "dashed") +
  geom_text(aes(x = 18, label = "\nTransition", y = 7.25),
            colour = "black", angle = 90, size = 4)

## mean concepts per encounter over age
df_encounters_age <- df %>%
  group_by(PatientId, GENEPOS_comb, ContactAge) %>%
  summarize(count_distinct = n_distinct(ConceptID)) %>%
  group_by(GENEPOS_comb, ContactAge) %>%
  summarize(count = mean(count_distinct))

# summary statistics of unique concepts per encounter
stats_encounters_age <- df_encounters_age %>%
  group_by(GENEPOS_comb) %>%
  summarize(mean = mean(count), median = median(count),
            sd = sd(count),
            min = min(count), max = max(count))

# t-test: mean difference of unique concepts per encounters between groups
p_encounters_age <- df_encounters_age %>%
  mutate(GENEPOS_comb = as.factor(GENEPOS_comb)) %>%
  ungroup() %>%
  summarize(pval = t.test(count ~ GENEPOS_comb)$p.value)

# plot: unique concept ID over count for each group
p5 <- df_encounters_age %>%
  ggplot(aes(x = ContactAge, y = count, color = GENEPOS_comb, fill = GENEPOS_comb)) +
  geom_point() +
  geom_smooth() +
  theme_classic() +
  guides(fill = guide_legend(title = "Group"),
         color = guide_legend(title = "Group")) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(0, 25), expand = FALSE) +
  ylab("Mean concepts per encounter") +
  xlab("Age (years)")

### ENRICHMENT PLOTS -----------------------------------------------------------
## Group 0: non-genetic vs. non-genetic (null)
enrich0 <- df_match0 %>%
  enrichmentPlot(., ont_hpo, forest = TRUE)

enrich0$plot <- enrich0$plot +
  coord_fixed(xlim = c(0, .3), ylim = c(0, .3)) +
  ggtitle("Null") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

## Group 1: genetic vs. non-genetic
enrich1 <- df_match1 %>%
  enrichmentPlot(., ont_hpo, forest = TRUE)

# manual labels
enrich1$plot$data$expcat_text <- NA
enrich1$plot$data[enrich1$plot$data$description == "Abnormality of the genitourinary system", ]$expcat_text <- "Abnormality of the genitourinary system"
enrich1$plot$data[enrich1$plot$data$description == "Intracranial hemorrhage", ]$expcat_text <- "Intracranial hemorrhage"
enrich1$plot$data[enrich1$plot$data$description == "Behavioral abnormality", ]$expcat_text <- "Behavioral abnormality"

enrich1$plot <- enrich1$plot +
  coord_fixed(xlim = c(0, .2), ylim = c(0, .2)) +
  ggtitle("Genetic vs. Non-Genetic") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

## Group 2: SCN1A vs. Genetic
enrich2 <- df_match2 %>%
  enrichmentPlot(., ont_hpo, forest = TRUE)

# manual labels
enrich2$plot$data$expcat_text <- NA
enrich2$plot$data[enrich2$plot$data$description == "Seizure", ]$expcat_text <- "Seizure"
enrich2$plot$data[enrich2$plot$data$description == "Abnormality of movement", ]$expcat_text <- "Abnormality of movement"
enrich2$plot$data[enrich2$plot$data$description == "Abnormality of the cardiovascular system", ]$expcat_text <- "Abnormality of the cardiovascular system"
enrich2$plot$data[enrich2$plot$data$description == "Infection-related seizure", ]$expcat_text <- "Infection-related seizure"

enrich2$plot <- enrich2$plot +
  coord_fixed(xlim = c(0, .2), ylim = c(0, .2)) +
  ggtitle("SCN1A vs. Genetic") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

## Group 3: CDKL5 vs. Genetic
enrich3 <- df_match3 %>%
  enrichmentPlot(., ont_hpo, forest = TRUE)

# manual labels
enrich3$plot$data$expcat_text <- NA
enrich3$plot$data[enrich3$plot$data$description == "Arrhythmia", ]$expcat_text <- "Arrhythmia"
enrich3$plot$data[enrich3$plot$data$description == "Abnormality of movement", ]$expcat_text <- "Abnormality of movement"
enrich3$plot$data[enrich3$plot$data$description == "Abnormal inflammatory response", ]$expcat_text <- "Abnormal inflammatory response"
enrich3$plot$data[enrich3$plot$data$description == "Low levels of vitamin D", ]$expcat_text <- "Low levels of vitamin D"

enrich3$plot <- enrich3$plot +
  coord_fixed(xlim = c(0, .5), ylim = c(0, .5)) +
  ggtitle("CDKL5 vs. Genetic") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

## Group 4: SCN1A vs. Non-Genetic
enrich4 <- df_match4 %>%
  enrichmentPlot(., ont_hpo, forest = TRUE)

# manual labels
enrich4$plot$data$expcat_text <- NA
enrich4$plot$data[enrich4$plot$data$description == "Seizure", ]$expcat_text <- "Seizure"
enrich4$plot$data[enrich4$plot$data$description == "Morphological central nervous system abnormality", ]$expcat_text <- "Morphological central nervous system abnormality"
enrich4$plot$data[enrich4$plot$data$description == "Muscle weakness", ]$expcat_text <- "Muscle weakness"
enrich4$plot$data[enrich4$plot$data$description == "Infection-related seizure", ]$expcat_text <- "Infection-related seizure"

enrich4$plot <- enrich4$plot +
  coord_fixed(xlim = c(0, .35), ylim = c(0, .35)) +
  ggtitle("SCN1A vs. Non-Genetic") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

## Group 5: CDKL5 vs. Non-Genetic
enrich5 <- df_match5 %>%
  enrichmentPlot(., ont_hpo, forest = TRUE)

# manual labels
enrich5$plot$data$expcat_text <- NA
enrich5$plot$data[enrich5$plot$data$description == "Arrhythmia", ]$expcat_text <- "Arrhythmia"
enrich5$plot$data[enrich5$plot$data$description == "Abnormality of the nervous system", ]$expcat_text <- "Abnormality of the nervous system"
enrich5$plot$data[enrich5$plot$data$description == "Abnormal inflammatory response", ]$expcat_text <- "Abnormal inflammatory response"
enrich5$plot$data[enrich5$plot$data$description == "Low levels of vitamin D", ]$expcat_text <- "Low levels of vitamin D"

enrich5$plot <- enrich5$plot +
  coord_fixed(xlim = c(0, .4), ylim = c(0, .4)) +
  ggtitle("CDKL5 vs. Non-Genetic") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

## Group 6: TSC vs. Genetic
enrich6 <- df_match6 %>%
  enrichmentPlot(., ont_hpo, forest = TRUE)

# manual labels
enrich6$plot$data$expcat_text <- NA
enrich6$plot$data[enrich6$plot$data$description == "Abnormal cardiovascular system physiology", ]$expcat_text <- "Abnormal cardiovascular system physiology"
enrich6$plot$data[enrich6$plot$data$description == "Visual field defect", ]$expcat_text <- "Visual field defect"
enrich6$plot$data[enrich6$plot$data$description == "Abnormality of the digestive system", ]$expcat_text <- "Abnormality of the digestive system"
enrich6$plot$data[enrich6$plot$data$description == "Renal cyst", ]$expcat_text <- "Renal cyst"

enrich6$plot <- enrich6$plot +
  coord_fixed(xlim = c(0, .4), ylim = c(0, .4)) +
  ggtitle("TSC vs. Genetic") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

## Group 7: TSC vs. Non-Genetic
enrich7 <- df_match7 %>%
  enrichmentPlot(., ont_hpo, forest = TRUE)

# manual labels
enrich7$plot$data$expcat_text <- NA
enrich7$plot$data[enrich7$plot$data$description == "Abnormal cardiovascular system physiology", ]$expcat_text <- "Abnormal cardiovascular system physiology"
enrich7$plot$data[enrich7$plot$data$description == "Visual field defect", ]$expcat_text <- "Visual field defect"
enrich7$plot$data[enrich7$plot$data$description == "Abnormality of the digestive system", ]$expcat_text <- "Abnormality of the digestive system"
enrich7$plot$data[enrich7$plot$data$description == "Proteinuria", ]$expcat_text <- "Proteinuria"

enrich7$plot <- enrich7$plot +
  coord_fixed(xlim = c(0, .4), ylim = c(0, .4)) +
  ggtitle("TSC vs. Non-Genetic") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

## Group 8: SCN1A vs. CDKL5
enrich8 <- df_match8 %>%
  enrichmentPlot(., ont_hpo, forest = TRUE)

# manual labels
enrich8$plot$data$expcat_text <- NA
enrich8$plot$data[enrich8$plot$data$description == "Arrhythmia", ]$expcat_text <- "Arrhythmia"
enrich8$plot$data[enrich8$plot$data$description == "Abnormality of the nervous system", ]$expcat_text <- "Abnormality of the nervous system"
enrich8$plot$data[enrich8$plot$data$description == "Infection-related seizure", ]$expcat_text <- "Infection-related seizure"
enrich8$plot$data[enrich8$plot$data$description == "Abnormality of the immune system", ]$expcat_text <- "Abnormality of the immune system"

enrich8$plot <- enrich8$plot +
  coord_fixed(xlim = c(0, .5), ylim = c(0, .5)) +
  ggtitle("SCN1A vs. CDKL5") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

### LONGITUDINAL ANALYSIS -----------------------------------------------------
## Group 0
plong0 <- longitudinalPlot(df_genes, df_match0, odds_plot = TRUE)

## Group 1
plong1 <- longitudinalPlot(df_genes, df_match1, odds_plot = TRUE)

## Group 2
plong2 <- longitudinalPlot(df_genes, df_match2, odds_plot = TRUE)

## Group 3
plong3 <- longitudinalPlot(df_genes, df_match3, odds_plot = TRUE)

## Group 4
plong4 <- longitudinalPlot(df_genes, df_match4, odds_plot = TRUE)

## Group 5
plong5 <- longitudinalPlot(df_genes, df_match5, odds_plot = TRUE)

## Group 6 
plong6 <- longitudinalPlot(df_genes, df_match6, odds_plot = TRUE)

## Group 7
plong7 <- longitudinalPlot(df_genes, df_match7, odds_plot = TRUE)

## Longitudinal heatmap analysis
plong1f <- longitudinalPlot(df_genes, df_match1, odds_plot = TRUE, filter = FALSE)

# set legend expression (for subscript)
str_legend <- expression(log['10']*(OR))

# set terms of interest from discovery graph (plong)
vec_longterms <- c("HP:0000708", # behav. abnormality
                   "HP:0002719", # recurrent infections
                   "HP:0002311", # incoordination
                   "HP:0011968", # feeding difficulties
                   "HP:0002019", # constipation
                   "HP:0000939", # osteoporosis
                   "HP:0002664", # neoplasm incl. neurofibromas
                   "HP:0012418", # hypoxemia
                   "HP:0000083", # renal insufficiency
                   ### negative controls ###
                   "HP:0003074", # Hyperglycemia
                   "HP:0003193", # Allergic rhinitis
                   "HP:0000388") # Otitis media

## heatmap label clustering
# get list of numeric vectors of odds ratios by term
dist_ls <- plong1f$plot$data %>%
  ungroup() %>%
  filter(term %in% vec_longterms) %>%
  mutate(odds = ifelse(pvalue > 0.05, NA, odds)) %>%
  select(term, odds) %>%
  split(by = "term") %>%
  lapply(., function(x) {as.numeric(x$odds)})

# pad vector length to equal, rbind into df
dist_df <- lapply(lapply(sapply(dist_ls, unlist), "length<-", max(lengths(dist_ls))), as.list) %>%
  lapply(., function(x) {unlist(x); as.data.frame(t(x))}) %>%
  rbindlist()

# get matrix, fix NA and Inf
dist_mat <- as.matrix(dist_df)
dist_mat[which(is.na(dist_mat))] <- 0
dist_mat[dist_mat == Inf] <- 99

# calculate distance and clustering
dist_dist <- dist(dist_mat, method = "euclidean")
dist_clust <- hclust(dist_dist)
## end of heatmap label clustering

pheat1 <- plong1f$plot$data %>%
  ungroup() %>%
  # select terms of interest from discovery graph (plong)
  filter(term %in% vec_longterms) %>%
  mutate(odds = ifelse(pvalue > 0.05, NA, odds)) %>%
  # set description factor levels
  mutate(description = as.factor(description)) %>%
  mutate(description = factor(description, levels = levels(description)[tmp_clust$order])) %>%
  # plot
  ggplot(aes(x = as.factor(bin), 
             y = description,
             # y = forcats::fct_relevel(description, rev),
             fill = log10(odds))) + 
  geom_tile(color = "black") +
  geom_text(aes(label = round(odds, 1)), size = 4) +
  theme_classic() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = NA),
        panel.spacing = unit(0.2, "cm")) +
  ylab("") +
  xlab("Median age at encounter (years)") +
  scale_fill_gradient2(name = str_legend,
                       low = "#00b4fb",
                       mid = "#F5F5F5",
                       high = "#ff8422",
                       midpoint = 0, # adjust based on OR or log OR
                       na.value = "#F5F5F5",
                       limits = c(-2, 2),
                       oob = scales::oob_squish_any)

### MEDICATION ANALYSIS -------------------------------------------------------
## Data and Preprocessing
# load data: list of patient medical record IDs and prescriptions
df_asm <-  readxl::read_excel("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/Medication_history.xlsx") %>%
  filter(MED_PHARM_CLASS_DESC %like% "ANTICONVULSANT")

# load data: FDA NDC database
df_fda <- readxl::read_excel("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/ndcxls.xlsx")
df_fda <- df_fda[df_fda$PHARM_CLASSES %like% "epilept", ] %>%
  distinct(PROPRIETARYNAME, NONPROPRIETARYNAME)

# source: AES Summary of Antiseizure Medications Available in the United States: 2020 Update
# https://www.aesnet.org/docs/default-source/pdfs-clinical/2020-september-aes_summary_of_asms.pdf?sfvrsn=c1a0ed0b_2
asm_vec <- readxl::read_excel("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/ASM list.xlsx", col_names = FALSE) %>%
  pull()

# match by AES and FDA lookup tables
for(i in 1:length(asm_vec)){
  ind_asm <- agrep(pattern = asm_vec[[i]],
                   x = df_asm$MED_NAME,
                   max.distance = 1,
                   ignore.case = TRUE)
  
  df_asm[ind_asm,]$MED_NAME <- asm_vec[[i]]
}

for(i in 1:nrow(df_fda)){
  ind_asm <- agrep(pattern = df_fda$PROPRIETARYNAME[[i]],
                   x = df_asm$MED_NAME,
                   max.distance = 1,
                   ignore.case = TRUE)
  
  df_asm[ind_asm,]$MED_NAME <- df_fda$NONPROPRIETARYNAME[[i]]
}

# manual recoding
vec_asm <- df_asm %>%
  mutate(MED_NAME = case_when(str_detect(MED_NAME, regex('vimpat', ignore_case = T)) ~ 'lacosamide',
                              str_detect(MED_NAME, regex('lacosamide', ignore_case = T)) ~ 'lacosamide', 
                              str_detect(MED_NAME, regex('banzel', ignore_case = T)) ~ 'rufinamide',
                              str_detect(MED_NAME, regex('tegretol', ignore_case = T)) ~ 'carbamazepine',
                              str_detect(MED_NAME, regex('carbatrol', ignore_case = T)) ~ 'carbamazepine',
                              str_detect(MED_NAME, regex('topamax', ignore_case = T)) ~ 'topiramate',
                              str_detect(MED_NAME, regex('topiram', ignore_case = T)) ~ 'topiramate',
                              str_detect(MED_NAME, regex('qudexy', ignore_case = T)) ~ 'topiramate',
                              str_detect(MED_NAME, regex('trokendi', ignore_case = T)) ~ 'topiramate',
                              str_detect(MED_NAME, regex('eprontia', ignore_case = T)) ~ 'topiramate',
                              str_detect(MED_NAME, regex('spritam', ignore_case = T)) ~ 'levetiracetam',
                              str_detect(MED_NAME, regex('keppra', ignore_case = T)) ~ 'levetiracetam',
                              str_detect(MED_NAME, regex('levetir', ignore_case = T)) ~ 'levetiracetam',
                              str_detect(MED_NAME, regex('sympazan', ignore_case = T)) ~ 'clobazam',
                              str_detect(MED_NAME, regex('onfi', ignore_case = T)) ~ 'clobazam',
                              str_detect(MED_NAME, regex('loraze', ignore_case = T)) ~ 'lorazepam',
                              str_detect(MED_NAME, regex('midazol', ignore_case = T)) ~ 'midazolam',
                              str_detect(MED_NAME, regex('nayzilam', ignore_case = T)) ~ 'midazolam',
                              str_detect(MED_NAME, regex('klonopin', ignore_case = T)) ~ 'clonazepam',
                              str_detect(MED_NAME, regex('lyrica', ignore_case = T)) ~ 'pregabalin',
                              str_detect(MED_NAME, regex('pregabalin', ignore_case = T)) ~ 'pregabalin',
                              str_detect(MED_NAME, regex('zonis', ignore_case = T)) ~ 'zonisamide',
                              str_detect(MED_NAME, regex('zoneg', ignore_case = T)) ~ 'zonisamide',
                              str_detect(MED_NAME, regex('valtoco', ignore_case = T)) ~ 'diazepam',
                              str_detect(MED_NAME, regex('diastat', ignore_case = T)) ~ 'diazepam',
                              str_detect(MED_NAME, regex('valpro', ignore_case = T)) ~ 'valproic acid',
                              str_detect(MED_NAME, regex('depak', ignore_case = T)) ~ 'valproic acid',
                              str_detect(MED_NAME, regex('lacosamide', ignore_case = T)) ~ 'lacosamide',
                              str_detect(MED_NAME, regex('fycompa', ignore_case = T)) ~ 'perampanel',
                              str_detect(MED_NAME, regex('perampanel', ignore_case = T)) ~ 'perampanel',
                              str_detect(MED_NAME, regex('fintepla', ignore_case = T)) ~ 'fintepla',
                              str_detect(MED_NAME, regex('fenfluramine', ignore_case = T)) ~ 'fenfluramine',
                              str_detect(MED_NAME, regex('lamot', ignore_case = T)) ~ 'lamotrigine',
                              str_detect(MED_NAME, regex('oxc', ignore_case = T)) ~ 'oxcarbazepine',
                              str_detect(MED_NAME, regex('trileptal', ignore_case = T)) ~ 'oxcarbazepine',
                              str_detect(MED_NAME, regex('loraz', ignore_case = T)) ~ 'lorazepam',
                              str_detect(MED_NAME, regex('vigabatr', ignore_case = T)) ~ 'vigabatrin',
                              str_detect(MED_NAME, regex('tiagab', ignore_case = T)) ~ 'tiagabine',
                              str_detect(MED_NAME, regex('rufinam', ignore_case = T)) ~ 'rufinamide',
                              str_detect(MED_NAME, regex('primid', ignore_case = T)) ~ 'primidone',
                              str_detect(MED_NAME, regex('phenyt', ignore_case = T)) ~ 'phenytoin',
                              str_detect(MED_NAME, regex('phenobarb', ignore_case = T)) ~ 'phenobarbital',
                              str_detect(MED_NAME, regex('gabapent', ignore_case = T)) ~ 'gabapentin',
                              str_detect(MED_NAME, regex('felbamat', ignore_case = T)) ~ 'felbamate',
                              str_detect(MED_NAME, regex('everol', ignore_case = T)) ~ 'everolimus',
                              str_detect(MED_NAME, regex('ethosux', ignore_case = T)) ~ 'ethosumixide',
                              str_detect(MED_NAME, regex('zaront', ignore_case = T)) ~ 'ethosumixide',
                              str_detect(MED_NAME, regex('cenoba', ignore_case = T)) ~ 'cenobamate',
                              str_detect(MED_NAME, regex('xcopri', ignore_case = T)) ~ 'cenobamate',
                              str_detect(MED_NAME, regex('cannab', ignore_case = T)) ~ 'cannabidiol',
                              str_detect(MED_NAME, regex('epidiolex', ignore_case = T)) ~ 'cannabidiol',
                              str_detect(MED_NAME, regex('briv', ignore_case = T)) ~ 'brivaracetam',
                              str_detect(MED_NAME, regex('stiri', ignore_case = T)) ~ 'stiripentol',
                              str_detect(MED_NAME, regex('celontin', ignore_case = T)) ~ 'mesuximide',
                              str_detect(MED_NAME, regex('eslicarbazepine', ignore_case = T)) ~ 'eslicarbazepine acetate')) %>%
  pull(MED_NAME)

# merge matches, drop non-matched rows and unneeded columns
df_asm$MED_NAME <- vec_asm

df_asm <- df_asm[df_asm$MED_NAME %in% asm_vec, ] %>%
  select(MedicalRecordNumber, MED_NAME, MED_START_DATE, ORD_MODE_DESC)

# merge in Patient Ids and DOB by MRN, fix date
df_lookup <- df_raw %>%
  distinct(PatientId, MedicalRecordNumber, DateOfBirth)

df_med <- df_asm %>%
  left_join(df_lookup[ ,c("PatientId", "MedicalRecordNumber", "DateOfBirth")], by = "MedicalRecordNumber") %>%
  mutate(DateOfBirth = lubridate::mdy(DateOfBirth)) %>%
  mutate(MED_START_DATE = lubridate::as_date(MED_START_DATE)) %>%
  mutate(AgePrescription = difftime(MED_START_DATE, DateOfBirth, units = "days")) %>%
  mutate(MonthsPrescription = as.numeric(round(AgePrescription/30, 0))) %>%
  mutate(YearsPrescription = as.numeric(AgePrescription/365.2425))
# mutate(YearsPrescription = as.numeric(round(AgePrescription/365, 0)))

## Descriptive stats
df_med <- df_med %>%
  group_by(PatientId) %>%
  mutate(n_unique_asm = n_distinct(MED_NAME)) %>% # number of unique ASMs per patient
  mutate(n_all_prescriptions = n_distinct(MonthsPrescription)) # number of all prescriptions per patient

## Group analysis: Heatmap (OR)
# define ASM groups, cf. doi.org/10.1007/s40263-021-00827-8
asm_map <- tibble(MED_NAME = asm_vec,
                  MED_GROUP = NA) %>%
  mutate(MED_GROUP = case_when(
    str_detect(MED_NAME, regex('felbam', ignore_case = T)) ~ 'Mixed/unknown',
    str_detect(MED_NAME, regex('valpr', ignore_case = T)) ~ 'Mixed/unknown',
    str_detect(MED_NAME, regex('hormon', ignore_case = T)) ~ 'Mixed/unknown',
    str_detect(MED_NAME, regex('zonis', ignore_case = T)) ~ 'Mixed/unknown',
    str_detect(MED_NAME, regex('rufin', ignore_case = T)) ~ 'Mixed/unknown',
    str_detect(MED_NAME, regex('cannab', ignore_case = T)) ~ 'Mixed/unknown',
    str_detect(MED_NAME, regex('cenob', ignore_case = T)) ~ 'Mixed/unknown',
    str_detect(MED_NAME, regex('topira', ignore_case = T)) ~ 'Mixed/unknown',
    str_detect(MED_NAME, regex('phenyt', ignore_case = T)) ~ 'Na',
    str_detect(MED_NAME, regex('carbamaz', ignore_case = T)) ~ 'Na',
    str_detect(MED_NAME, regex('carbaz', ignore_case = T)) ~ 'Na',
    str_detect(MED_NAME, regex('lamotr', ignore_case = T)) ~ 'Na',
    str_detect(MED_NAME, regex('lacosam', ignore_case = T)) ~ 'Na',
    str_detect(MED_NAME, regex('suxim', ignore_case = T)) ~ 'Ca',
    str_detect(MED_NAME, regex('gabap', ignore_case = T)) ~ 'Ca',
    str_detect(MED_NAME, regex('pregabal', ignore_case = T)) ~ 'Ca',
    str_detect(MED_NAME, regex('phenobarb', ignore_case = T)) ~ 'GABA',
    str_detect(MED_NAME, regex('primidon', ignore_case = T)) ~ 'GABA',
    str_detect(MED_NAME, regex('stirip', ignore_case = T)) ~ 'GABA',
    str_detect(MED_NAME, regex('azol', ignore_case = T)) ~ 'GABA',
    str_detect(MED_NAME, regex('azepam', ignore_case = T)) ~ 'GABA',
    str_detect(MED_NAME, regex('clobaz', ignore_case = T)) ~ 'GABA',
    str_detect(MED_NAME, regex('tiaga', ignore_case = T)) ~ 'GABA',
    str_detect(MED_NAME, regex('vigabatr', ignore_case = T)) ~ 'GABA',
    str_detect(MED_NAME, regex('peramp', ignore_case = T)) ~ 'AMPA',
    str_detect(MED_NAME, regex('brivara', ignore_case = T)) ~ 'SV2A',
    str_detect(MED_NAME, regex('levetir', ignore_case = T)) ~ 'SV2A',
    str_detect(MED_NAME, regex('fenfl', ignore_case = T)) ~ '5-HT',
    str_detect(MED_NAME, regex('everol', ignore_case = T)) ~ 'MTOR'
  ) )

# take strict matched case-control set
df_heatmap <- df_match1 %>%
  # merge in patient ID and group label
  distinct(PatientId, group) %>%
  left_join(df_med, by = "PatientId") %>%
  na.omit %>%
  # define bins
  mutate(YearsPrescription = cut_number(YearsPrescription, 4)) %>%
  # count ASM prescription per group; for each age bin (year)
  group_by(group, MED_NAME, YearsPrescription) %>%
  summarize(test = n()) %>%
  # prepare Fisher's test
  pivot_wider(names_from = group, values_from = test) %>%
  rename(N = `FALSE`, Y = `TRUE`) %>%
  # force complete all cases
  complete(MED_NAME, YearsPrescription) %>%
  replace(is.na(.), 0) %>%
  mutate(N_out = max(N)-N, 
         Y_out = max(Y)-Y) %>%
  rowwise() %>%
  # do Fisher's test
  mutate(P = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$p.value,
         OR = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$estimate,
         CI1 = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$conf.int[[1]],
         CI2 = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$conf.int[[2]]) %>%
  # adjust for multiple testing
  mutate(P = p.adjust(P, method = "bonferroni")) %>%
  # force OR with insignificant pvalues to be NA
  mutate(OR = ifelse(P > 0.05, NA, OR)) %>%
  # merge in group description for facets
  left_join(asm_map, by = "MED_NAME") %>%
  # define factor levels
  mutate(MED_NAME = factor(MED_NAME, levels = sort(unique(asm_vec)))) %>%
  # frequency filter
  ungroup() %>%
  mutate(freqY = Y/max(Y)) %>%
  mutate(freqN = N/max(N)) %>%
  filter(freqY > 0.01 & freqN > 0.01)

# plot
p_asm <- df_heatmap %>%
  mutate(MED_NAME = str_to_title(MED_NAME)) %>%
  ggplot(aes(x = YearsPrescription, 
             y = MED_NAME,
             fill = log10(OR))) + 
  geom_tile(color = "black") +
  geom_text(aes(label = round(OR, 1)), size = 4) +
  theme_classic() +
  facet_grid(MED_GROUP~., scales = "free_y", space = "free", 
             switch = "y", drop = FALSE, margins = FALSE) +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = NA),
        panel.spacing = unit(0.2, "cm")) +
  ylab("") +
  xlab("Age at prescription (years)") +
  scale_x_discrete(labels = c("0-2", "2-4", "4-8", ">8")) +
  scale_fill_gradient2(name = str_legend,
                       low = "#00b4fb",
                       mid = "#F5F5F5",
                       high = "#ff8422",
                       midpoint = 0, # adjust based on OR or log OR
                       na.value = "#F5F5F5",
                       limits = c(-2, 2),
                       oob = scales::oob_squish_any)

### NON-HPO CONCEPT ANALYSIS ---------------------------------------------------
# the basic concept: find out how many UMLS concepts do not match to HPO terms
# find patterns in these concepts, e.g. healthcare utilization or procedures

## data
# load UMLS concepts
# source: nlm.nih.gov/research/umls/
umls_map <- read.delim("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/MRCONSO.RRF",
                       sep = "|", header = FALSE)

umls_map <- umls_map %>%
  as_tibble() %>%
  filter(V2 == "ENG") %>%
  distinct(V1, V15) %>%
  rename(ConceptID = V1, ConceptDesc = V15)

# only keep one description per term (multiple vocabularies)
umls_map <- umls_map %>%
  group_by(ConceptID) %>%
  filter(row_number() == 1)

# anti-join UMLS vs HPO
df_concepts <- df_genes %>%
  ungroup %>%
  anti_join(hpo_map, by = "ConceptID")

# get concept descriptions
# note: MRCONSO.RRF represents only a common subset of concepts; the full set of
# the metathesaurus is impractically large for the purpose of this analysis
# some manual annotation downstream will be necessary
df_concepts <- df_concepts %>%
  left_join(umls_map, by = "ConceptID")

# save total number of concepts and concepts mapped to HPO
stats_concepts <- tibble(n_all = nrow(df_genes),
                         n_nonhpo = nrow(df_concepts)) %>%
  mutate(ratio = n_nonhpo/n_all, diff = n_all-n_nonhpo)

# find and save list of most common concepts
df_commonconcepts <- df_concepts %>%
  count(ConceptID, ConceptDesc) %>%
  arrange(desc(n))

# subset to matched cohort; cross-sectional for now
df_conceptmatch <- df_match1 %>%
  distinct(PatientId, group) %>%
  left_join(df_concepts[, c("PatientId", "ConceptID", "ConceptDesc")], by = "PatientId") %>%
  distinct()

# count by group, then do Fisher's test
df_conceptmatch <- df_conceptmatch %>%
  group_by(group) %>%
  count(ConceptID) %>%
  pivot_wider(names_from = group, values_from = n) %>%
  rename(N = `FALSE`, Y = `TRUE`) %>%
  replace(is.na(.), 0) %>%
  mutate(N_out = max(N)-N, Y_out = max(Y)-Y) %>%
  rowwise() %>%
  # do Fisher's test
  mutate(P = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$p.value,
         OR = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$estimate,
         CI1 = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$conf.int[[1]],
         CI2 = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$conf.int[[2]]) %>%
  # adjust for multiple testing, keep original values
  ungroup() %>%
  mutate(P_i = P) %>%
  mutate(P = p.adjust(P, method = "bonferroni"))

## QQ plot
pqq <- gg_qqplot(df_conceptmatch$P_i) +
  theme_classic() +
  annotate(geom = "text", x = -Inf, y = Inf,
           hjust = -0.15, vjust = 1 + 0.15 * 3,
           label = sprintf("  λ    = %.2f", inflation(df_conceptmatch$P_i)),
           size = 5) +
  coord_equal() +
  theme(aspect.ratio = 1) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

# keep significant associations
df_conceptmatch <- df_conceptmatch %>%
  filter(P < 0.05)

# get descriptions
df_conceptmatch <- df_conceptmatch %>%
  left_join(umls_map, by = "ConceptID")

# manual annotation
df_conceptmatch[df_conceptmatch$ConceptID == "C0478107", ]$ConceptDesc <- NA
df_conceptmatch[df_conceptmatch$ConceptID == "C0476431", ]$ConceptDesc <- NA
# df_conceptmatch[df_conceptmatch$ConceptID == "C0478107", ]$ConceptDesc <- "Other specified chromosome abnormalities"
# df_conceptmatch[df_conceptmatch$ConceptID == "C0476431", ]$ConceptDesc <- "Abnormal karyotype"
df_conceptmatch[df_conceptmatch$ConceptID == "C2875116", ]$ConceptDesc <- "Generalized epilepsy and epileptic syndromes, intractable"
df_conceptmatch[df_conceptmatch$ConceptID == "C2910620", ]$ConceptDesc <- "Screening for cardiovascular disorders"
df_conceptmatch[df_conceptmatch$ConceptID == "C3161331", ]$ConceptDesc <- "Unspecified intellectual disabilities"
df_conceptmatch[df_conceptmatch$ConceptID == "C0341102", ]$ConceptDesc <- "Gastroesophageal reflux disease"
df_conceptmatch[df_conceptmatch$ConceptID == "C2911172", ]$ConceptDesc <- "Other specified health status"
df_conceptmatch[df_conceptmatch$ConceptID == "C2911188", ]$ConceptDesc <- "Other long term drug therapy"

# forest plot
p_forest_nonhpo <- df_conceptmatch %>%
  na.omit %>%
  slice_max(order_by = OR, n = 8, with_ties = FALSE) %>%
  ggplot(aes(y = reorder(ConceptDesc, OR))) +
  geom_point(aes(x = OR), shape = 15, size = 3) +
  geom_linerange(aes(xmin = CI1, xmax = CI2)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(trans = 'log10',
                     oob = scales::oob_squish_infinite) +
  scale_y_discrete(labels = label_wrap(30)) +
  expand_limits(x = 1) +
  theme_classic() +
  ylab("") +
  xlab("Odds ratio (95% CI, log scale)")

### TRANSITION ANALYSIS --------------------------------------------------------
# the idea: find what drives the increase in encounters for likely genetic
# patients during the transition to adult care (ages 18-20 years)

## get data
df_trans <- df_genes %>%
  filter(status %in% c("genetic", "nongenetic")) %>%
  distinct(PatientId, ConceptID, ContactAge, status) %>%
  mutate(age_group = case_when(ContactAge > 18 & ContactAge < 20 ~ 1,
                               ContactAge < 18 ~ 0,
                               TRUE ~ NA_real_)) %>%
  na.omit

# analyze
df_trans <- df_trans %>%
  group_by(age_group, status) %>%
  count(ConceptID) %>%
  pivot_wider(names_from = age_group, values_from = n) %>%
  rename(N = `0`, Y = `1`) %>%
  replace(is.na(.), 0) %>%
  mutate(N_out = max(N)-N, Y_out = max(Y)-Y) %>%
  rowwise() %>%
  # do Fisher's test
  mutate(P = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$p.value,
         OR = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$estimate,
         CI1 = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$conf.int[[1]],
         CI2 = fisher.test(matrix(c(Y, Y_out, N, N_out), nrow = 2, ncol = 2))$conf.int[[2]]) %>%
  # adjust for multiple testing
  ungroup() %>%
  mutate(P = p.adjust(P, method = "bonferroni"))

# keep significant associations
df_trans <- df_trans %>%
  filter(P < 0.05)

# get descriptions
df_trans <- df_trans %>%
  left_join(umls_map, by = "ConceptID")
# note: compare to non-genetic cohort, sorted by p-value or frequency in the positive transition group

# manual annotation
df_trans[df_trans$ConceptID == "C0260698", ]$ConceptDesc <- "Other postprocedural status"
df_trans[df_trans$ConceptID == "C0036421", ]$ConceptDesc <- "Systemic Scleroderma"
df_trans[df_trans$ConceptID == "C0260860", ]$ConceptDesc <- "Encounter due to Unspecified general medical examination"
df_trans[df_trans$ConceptID == "C0477590", ]$ConceptDesc <- "Other overlap syndromes"
df_trans[df_trans$ConceptID == "C2900579", ]$ConceptDesc <- "Age-related osteoporosis without current pathological fracture"
df_trans[df_trans$ConceptID == "C0260545", ]$ConceptDesc <- "examination; infant or child"
df_trans[df_trans$ConceptID == "C2886562", ]$ConceptDesc <- "Unspecified child maltreatment, confirmed, initial encounter"
df_trans[df_trans$ConceptID == "C2852675", ]$ConceptDesc <- "(...), initial encounter for closed fracture"
df_trans[df_trans$ConceptID == "C2863970", ]$ConceptDesc <- NA
df_trans[df_trans$ConceptID == "C2852166", ]$ConceptDesc <- NA
df_trans[df_trans$ConceptID == "C2868158", ]$ConceptDesc <- NA
df_trans[df_trans$ConceptID == "C0159791", ]$ConceptDesc <- NA
df_trans[df_trans$ConceptID == "C0159906", ]$ConceptDesc <- NA
df_trans[df_trans$ConceptID == "C2868158", ]$ConceptDesc <- NA
df_trans[df_trans$ConceptID == "C2868158", ]$ConceptDesc <- NA
df_trans[df_trans$ConceptID == "C0152101", ]$ConceptDesc <- "Hypoplastic left heart syndrome"
df_trans[df_trans$ConceptID == "C0375114", ]$ConceptDesc <- "Diabetes mellitus, type I"
df_trans[df_trans$ConceptID == "C0494284", ]$ConceptDesc <- NA
df_trans[df_trans$ConceptID == "C0260698", ]$ConceptDesc <- NA
df_trans[df_trans$ConceptID == "C0004352", ]$ConceptDesc <- "Autistic disorder"
df_trans[df_trans$ConceptID == "C2911178", ]$ConceptDesc <- "Encounter for long-term use of anticoagulants"
df_trans[df_trans$ConceptID == "C0004352", ]$ConceptDesc <- "Autistic disorder"

# forest plot: genetic group
pt1 <- df_trans %>%
  filter(status == "genetic") %>%
  na.omit %>%
  slice_max(order_by = Y, n = 4, with_ties = FALSE) %>%
  ggplot(aes(y = reorder(ConceptDesc, OR))) +
  geom_linerange(aes(xmin = CI1, xmax = CI2)) +
  geom_point(aes(x = OR), shape = 15, size = 3, color = "#d95f02") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(trans = 'log10',
                     limits = c(1, 35),
                     oob = scales::oob_squish_infinite) +
  scale_y_discrete(labels = label_wrap(30)) +
  expand_limits(x = 1) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("Odds ratio (95% CI, log scale)")

# forest plot: nongenetic group
df_trans[df_trans$ConceptID == "C2887465", ]$ConceptDesc <- NA
df_trans[df_trans$ConceptID == "C0042870", ]$ConceptDesc <- "Vitamin D Deficiency"
df_trans[df_trans$ConceptID == "C2910447", ]$ConceptDesc <- "Encounter for general adult medical examination without abnormal findings"
df_trans[df_trans$ConceptID == "C2911563", ]$ConceptDesc <- "Vitamin D Deficiency"
df_trans[df_trans$ConceptID == "C0042870", ]$ConceptDesc <- "Other specified postprocedural states"
df_trans[df_trans$ConceptID == "C0154714", ]$ConceptDesc <- "Localization-related epilepsy, without mention of intractable epilepsy "
df_trans[df_trans$ConceptID == "C0042870", ]$ConceptDesc <- NA

df_trans[df_trans$ConceptID == "C0155886", ]$ConceptDesc <- "Asthma"
df_trans[df_trans$ConceptID == "C2910447", ]$ConceptDesc <- "Medical examination w/o abnormal findings"
df_trans[df_trans$ConceptID == "C0154714", ]$ConceptDesc <- "Focal epilepsy, non-intractable"

pt2 <- df_trans %>%
  filter(status == "nongenetic") %>%
  na.omit %>%
  slice_max(order_by = Y, n = 4, with_ties = FALSE) %>%
  ggplot(aes(y = reorder(ConceptDesc, OR))) +
  geom_linerange(aes(xmin = CI1, xmax = CI2)) +
  geom_point(aes(x = OR), shape = 15, size = 3, color = "#1b9e77") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(trans = 'log10',
                     limits = c(1, 35),
                     oob = scales::oob_squish_infinite) +
  scale_y_discrete(labels = label_wrap(30)) +
  expand_limits(x = 1) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("Odds ratio (95% CI)")

pt <- cowplot::plot_grid(pt1 + theme(axis.title.x = element_blank()), 
                         pt2, 
                         ncol = 1, align = "v")

### OTHER SUB-ANALYSIS ---------------------------------------------------------
## For HPO OR plot, find the terms that make up geniturourinary system abnormality
# define descendants
vec_desc <- get_descendants(ont_hpo, "HP:0000119")

df_desc <- df_match1 %>%
  # remove TSC patients, to see if association remains strong
  filter(PatientId %nin% df_tsc$PatientID) %>%
  # filter cases by those containing descendant terms
  filter(term %in% vec_desc) %>%
  filter(group == TRUE) %>%
  group_by(term) %>%
  count(sort = TRUE) %>%
  # get descriptions and ic
  left_join(desc_map, by = "term") %>%
  left_join(df_ic, by = "term") %>%
  ungroup() %>%
  # get pvalue for each term, only keep if it's independently significant
  left_join(enrich1$data[, c("term", "pvalue")], by = "term") %>%
  filter(pvalue < 0.05) %>%
  # sort by IC
  slice_max(n, prop = 0.9) %>%
  arrange(desc(ic)) 

# reduce to minimal set
vec_min <- minimal_set(ont_hpo, df_desc$term)
df_desc <- df_desc[df_desc$term %in% vec_min, ]

# restrict subgraph by propagating back up from our minimal set
vec_desc <- propagate_relations(ont_hpo, df_desc$term, "parents") %>% unique()

# make igraph, convert to bn
librarian::shelf(igraph)
parents <- ont_hpo$parents
self <- rep(names(parents), lengths(parents))
g <- igraph::make_graph(rbind(unlist(parents), self))
bng <- as.bn(g)

# generate HPO subgraph of node of interest
arcs <- data.frame(bng$arcs)
arcs <- arcs[arcs$X1 %in% vec_desc & arcs$X2 %in% vec_desc, ]
colnames(arcs) <- c("from", "to")
graph_hpo <- empty.graph(vec_desc)
arcs(graph_hpo) <- arcs

# revert to igraph and format graph
g2 <- as.igraph(graph_hpo)
## note: complete term name labels are just impossible to format with igraph
# V(g2)$label <- NA
# V(g2)[vec_min]$label <- paste0("\n", "\n", "\n", "\n", "\n", "\n", "\n", df_desc$description)
# V(g2)["HP:0003244"]$label <- paste0("\n","\n","\n","\n","\n", "\n", "\n", "\n", "\n", "\t","Penile hypospadia")
# V(g2)["HP:0000028"]$label <- paste0("\n", "\n", "\n", "\n", "\n", "\t","Cryptorchidism")
# V(g2)["HP:0000086"]$label <- paste0("\n", "\n", "\t","Ectopic kidney")
V(g2)$label <- NA
V(g2)[vec_min]$label <- rep(LETTERS)[1:length(V(g2)[vec_min]$label)]
V(g2)$color <- "gray"
V(g2)[vec_min]$color <- "red"
V(g2)$size <- 8
vec_size <- -log10(df_desc$pvalue)
vec_size[vec_size < 7] <- 8
V(g2)[vec_min]$size <- vec_size
E(g2)$arrow.mode <- 2

# use ggplotify to get a grob-able object
pqg <- as.ggplot(expression(plot(g2, 
                                 vertex.frame.color = "black",
                                 vertex.label.color = "white",
                                 vertex.label.family = "Helvetica",
                                 vertex.label.font = 1,
                                 vertex.label.cex = .9,
                                 # vertex.label.dist = 2.5,
                                 edge.arrow.size = .5,
                                 layout = layout_as_tree))) # Reingold-Tilford graph

pqg <- pqg + 
  theme(plot.margin = unit(c(-50, -20, -50, -50), "pt"))

### GENERATE REPORT ------------------------------------------------------------
## Figure 1: Descriptive statistics of the study cohort.
Fig1 <- cowplot::plot_grid(p1,
                           p2 + scale_x_discrete(labels=c("Non-genetic", "Likely genetic")),
                           p3 + scale_x_discrete(labels=c("Non-genetic", "Likely genetic")),
                           p5 + theme(legend.position = "none"),
                           p4 + theme(legend.position = "none"),
                           pt,
                           nrow = 2, labels = "AUTO", align = "none")

pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig1.pdf",
    width = 12,
    height = 8)

Fig1

dev.off()

## Figure 2: Genetic vs. Non-Genetic
Fig2 <- cowplot::plot_grid(pqq, 
                           p_forest_nonhpo,
                           enrich1$forest,
                           enrich1$plot + 
                             ggtitle("") + 
                             theme_set(theme_classic()) +
                             coord_cartesian(xlim = c(0, 0.15), ylim = c(0, 0.15)) +
                             ylab("Frequency, likely genetic patient encounters") +
                             xlab("Frequency, non-genetic patients encounters"),
                           enrich8$plot + 
                             ggtitle("") + 
                             theme_set(theme_classic()) +
                             coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 0.5)) +
                             ylab("Frequency, SCN1A patient encounters") +
                             xlab("Frequency, CDKL5 patient encounters"),
                           pqg,
                           nrow = 2, labels = "AUTO", align = "none")

pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig2.pdf",
    width = 12,
    height = 8,
    encoding = 'CP1253.enc') # to draw Lambda on panel A

Fig2

dev.off()

# Fig 3: Longitudinal heatmap of genetic vs. non-genetic patients
pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig3.pdf",
    width = 12,
    height = 4)

pheat1 

dev.off()

### Figure 4: Prescription patterns of genetic vs. non-genetic patients
pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig4.pdf",
    width = 8,
    height = 8)

p_asm 

dev.off()
