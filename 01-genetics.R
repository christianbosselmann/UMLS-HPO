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
                 scales)

# functions
source("func.R")

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
  xlab("Age at encounter") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(0, 35), expand = FALSE) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

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
  count(cut_width(ContactAge, width = 1, boundary = 0, labels = F)) 

df_encounters_freq <- df_encounters_freq %>%
  ungroup() %>%
  # maintain group label
  left_join(df_encounters[ ,c("GENEPOS_comb", "PatientId")], by = "PatientId") %>%
  # fix bin label
  rename(bin = `cut_width(ContactAge, width = 1, boundary = 0, labels = F)`) %>%
  group_by(bin, GENEPOS_comb) %>%
  # get mean number of encounters per bin per group
  summarize(mean = mean(n), sd = sd(n))

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
  coord_cartesian(xlim = c(0, 24), expand = FALSE) +
  ylab("Mean encounters per year") +
  xlab("Age (years)")

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
## Group 1: genetic vs. non-genetic
enrich1 <- df_match1 %>%
  enrichmentPlot(., ont_hpo, forest = TRUE)

# manual labels
enrich1$plot$data$expcat_text <- NA
enrich1$plot$data[enrich1$plot$data$description == "Abnormality of the genitourinary system", ]$expcat_text <- "Abnormality of the genitourinary system"
enrich1$plot$data[enrich1$plot$data$description == "Abnormality of the skeletal system", ]$expcat_text <- "Abnormality of the skeletal system"
enrich1$plot$data[enrich1$plot$data$description == "Intracranial hemorrhage", ]$expcat_text <- "Intracranial hemorrhage"

enrich1$plot <- enrich1$plot +
  coord_fixed(xlim = c(0, .3), ylim = c(0, .3)) +
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

### LONGITUDINAL ANALYSIS -----------------------------------------------------
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
                              str_detect(MED_NAME, regex('qudexy', ignore_case = T)) ~ 'topiramate',
                              str_detect(MED_NAME, regex('trokendi', ignore_case = T)) ~ 'topiramate',
                              str_detect(MED_NAME, regex('eprontia', ignore_case = T)) ~ 'topiramate',
                              str_detect(MED_NAME, regex('spritam', ignore_case = T)) ~ 'levetiracetam',
                              str_detect(MED_NAME, regex('keppra', ignore_case = T)) ~ 'levetiracetam',
                              str_detect(MED_NAME, regex('sympazan', ignore_case = T)) ~ 'clobazam',
                              str_detect(MED_NAME, regex('onfi', ignore_case = T)) ~ 'clobazam',
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
                              str_detect(MED_NAME, regex('viga', ignore_case = T)) ~ 'vigabatrin',
                              str_detect(MED_NAME, regex('tiagab', ignore_case = T)) ~ 'tiagabine',
                              str_detect(MED_NAME, regex('rufinam', ignore_case = T)) ~ 'rufinamide',
                              str_detect(MED_NAME, regex('primid', ignore_case = T)) ~ 'primidone',
                              str_detect(MED_NAME, regex('phenyt', ignore_case = T)) ~ 'phenytoin',
                              str_detect(MED_NAME, regex('phenobarb', ignore_case = T)) ~ 'phenobarbital',
                              str_detect(MED_NAME, regex('gaba', ignore_case = T)) ~ 'gabapentin',
                              str_detect(MED_NAME, regex('felb', ignore_case = T)) ~ 'felbamate',
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
  mutate(YearsPrescription = as.numeric(round(AgePrescription/365, 0)))

## Descriptive stats
df_med <- df_med %>%
  group_by(PatientId) %>%
  mutate(n_unique_asm = n_distinct(MED_NAME)) %>% # number of unique ASMs per patient
  mutate(n_all_prescriptions = n_distinct(MonthsPrescription)) # number of all prescriptions per patient

## Group analysis: Heatmap (absolute)
df_match2 %>%
  distinct(PatientId, group) %>%
  left_join(df_med, by = "PatientId") %>%
  filter(group == TRUE) %>%
  na.omit %>%
  # count the number of patients who have received each ASM by year
  group_by(MED_NAME, YearsPrescription) %>%
  mutate(freq = n_distinct(PatientId)) %>%
  ungroup() %>%
  # plot
  ggplot(aes(x = YearsPrescription, 
             y = factor(MED_NAME, levels = asm_vec[asm_vec %in% df_asm$MED_NAME]), 
             fill = freq)) + 
  scale_y_discrete(drop = FALSE) +
  geom_tile(color = "black") +
  theme_classic() +
  ylab("") +
  xlab("Age (years)") +
  scale_fill_gradientn(colors = hcl.colors(20, "Temps")) +
  coord_fixed()

### WIP
## Group analysis: Heatmap (OR)
df_match1 %>%
  # merge in patient ID and group label
  distinct(PatientId, group) %>%
  left_join(df_med, by = "PatientId") %>%
  na.omit %>%
  # # adjust binwidth here; line optional
  mutate(YearsPrescription = cut_number(YearsPrescription, 5)) %>%
  # mutate(YearsPrescription = cut(YearsPrescription, 
  #                                breaks=c(0, 6, 12, 18, 99), include.lowest=TRUE)) %>%
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
  # # force OR with insignificant pvalues to be NA
  # mutate(OR = ifelse(P > 0.05, NA, OR)) %>%
  # plot
  ggplot(aes(x = YearsPrescription, 
             y = factor(MED_NAME, levels = asm_vec[asm_vec %in% df_asm$MED_NAME]), 
             fill = log10(OR))) + 
  scale_y_discrete(drop = FALSE) +
  geom_tile(color = "black") +
  theme_classic() +
  ylab("") +
  xlab("Age (years)") +
  # scale_fill_gradientn(colors = hcl.colors(20, "Temps"),
  #                      na.value = "darkgrey") +
  scale_fill_gradient2(low = "#00b4fb",
                       mid = "#F5F5F5",
                       high = "#ff8422",
                       midpoint = 0, # adjust based on OR or log OR
                       na.value = "#F5F5F5",
                       limits = c(-2, 2),
                       oob = scales::squish) +
  coord_fixed()

### GENERATE REPORT -----------------------------------------------------------
## Figure 1: Descriptive statistics of the study cohort.
p_tmp <- cowplot::plot_grid(p2 + scale_x_discrete(labels=c("Non-genetic", "Genetic")), 
                            p3 + scale_x_discrete(labels=c("Non-genetic", "Genetic")), 
                            p4 + theme(legend.position = "none"), 
                            p5 + theme(legend.position = "none"), 
                            nrow = 2, labels = c("B", "C", "D", "E"))
Fig1 <- cowplot::plot_grid(p1, p_tmp, rel_widths = c(1/2, 1/2), labels = c("A", ""))

pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig1.pdf",
    width = 12,
    height = 6)

Fig1

dev.off()

## Figure 2: Genetic vs. Non-Genetic
p_tmp <- cowplot::plot_grid(enrich1$plot + 
                              ggtitle("") + 
                              theme_set(theme_classic()) +
                              coord_cartesian(xlim = c(0, 0.15), ylim = c(0, 0.15)) +
                              ylab("Genetic patients") +
                              xlab("Non-genetic patients"), 
                            enrich1$forest,
                            nrow = 1, labels = "AUTO")

Fig2 <- cowplot::plot_grid(p_tmp, plong1$plot, 
                           nrow = 2,
                           labels = c("", "C"))

pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig2.pdf",
    width = 12,
    height = 12)

Fig2

dev.off()

## Figure 3: SCN1A vs. Non-Genetic
p_tmp <- cowplot::plot_grid(enrich4$plot + 
                              ggtitle("") + 
                              theme_set(theme_classic()) +
                              coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
                              ylab("SCN1A") +
                              xlab("Non-genetic patients"), 
                            enrich4$forest,
                            nrow = 1, labels = "AUTO")

Fig3 <- cowplot::plot_grid(p_tmp, plong4$plot, 
                           nrow = 2,
                           labels = c("", "C"))

pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig3.pdf",
    width = 12,
    height = 12)

Fig3

dev.off()

## Figure 3.1: SCN1A vs. Genetic
p_tmp <- cowplot::plot_grid(enrich2$plot + 
                              ggtitle("") + 
                              theme_set(theme_classic()) +
                              coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
                              ylab("SCN1A") +
                              xlab("Genetic patients"), 
                            enrich2$forest,
                            nrow = 1, labels = "AUTO")

Fig3.1 <- cowplot::plot_grid(p_tmp, plong2$plot, 
                             nrow = 2,
                             labels = c("", "C"))

pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig3.1.pdf",
    width = 12,
    height = 12)

Fig3.1

dev.off()

## Figure 4: CDKL5 vs. Non-Genetic
p_tmp <- cowplot::plot_grid(enrich5$plot + 
                              ggtitle("") + 
                              theme_set(theme_classic()) +
                              coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
                              ylab("CDKL5") +
                              xlab("Non-genetic patients"), 
                            enrich5$forest,
                            nrow = 1, labels = "AUTO")

Fig4 <- cowplot::plot_grid(p_tmp, plong5$plot + coord_cartesian(xlim = c(0, 12)), 
                           nrow = 2,
                           labels = c("", "C"))

pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig4.pdf",
    width = 12,
    height = 12)

Fig4

dev.off()

## Figure 4.1: CDKL5 vs. Genetic
p_tmp <- cowplot::plot_grid(enrich3$plot + 
                              ggtitle("") + 
                              theme_set(theme_classic()) +
                              coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
                              ylab("CDKL5") +
                              xlab("Genetic patients"), 
                            enrich3$forest,
                            nrow = 1, labels = "AUTO")

Fig4.1 <- cowplot::plot_grid(p_tmp, plong3$plot + coord_cartesian(xlim = c(0, 12)), 
                             nrow = 2,
                             labels = c("", "C"))

pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig4.1.pdf",
    width = 12,
    height = 12)

Fig4.1

dev.off()

## Figure 5: TSC vs. Non-Genetic
p_tmp <- cowplot::plot_grid(enrich7$plot + 
                              ggtitle("") + 
                              theme_set(theme_classic()) +
                              coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
                              ylab("TSC") +
                              xlab("Non-genetic patients"), 
                            enrich7$forest,
                            nrow = 1, labels = "AUTO")

Fig5 <- cowplot::plot_grid(p_tmp, plong7$plot, 
                           nrow = 2,
                           labels = c("", "C"))

pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig5.pdf",
    width = 12,
    height = 12)

Fig5

dev.off()

## Figure 5.1: TSC vs. Genetic
p_tmp <- cowplot::plot_grid(enrich7$plot + 
                              ggtitle("") + 
                              theme_set(theme_classic()) +
                              coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
                              ylab("TSC") +
                              xlab("Genetic patients"), 
                            enrich7$forest,
                            nrow = 1, labels = "AUTO")

Fig5.1 <- cowplot::plot_grid(p_tmp, plong7$plot, 
                             nrow = 2,
                             labels = c("", "C"))

pdf(file = "/Users/cbosselmann/Desktop/GitHub/UMLS-HPO/out/pub_genetic/Fig5.1.pdf",
    width = 12,
    height = 12)

Fig5.1

dev.off()
