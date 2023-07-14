### UMLS-HPO ANALYSIS OF INDIVIDUALS WITH GENETIC EPILEPSY SYNDROMES -----------
## Validation cohort
##
## Author: Christian Bosselmann, MD
##
## Date Created: 2023-06-19
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
                 ggplotify,
                 pheatmap,
                 survminer,
                 RColorBrewer)

# functions
source("func.R")

# lookup tables
source("dict.R")

# seed
set.seed(42)

### PARAMETERS ----------------------------------------------------------------
# matching method for matchit
flag_match <- "nearest"

# color palette of analogous colours for the groups of the validity analysis
pal_val <- c("#1B9E35", "#D90210")

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

### DATA -----------------------------------------------------------------------
# data: all encounters per patient, ages 0-6, grouped by gene positive / negative
df_raw <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/longitudinal_genetic.csv")

# preprocessing
df <- df_raw %>% 
  select(PatientId, MedicalRecordNumber, ConceptID, GENEPOS_comb, ContactAge, ProcAge) %>%
  group_by(PatientId) %>%
  arrange(desc(ContactAge)) %>%
  fill(ContactAge, .direction = c("up")) %>%
  na.omit

# by-patient demographic data and age statistics
df_person <- df_raw %>%
  group_by(PatientId, DateOfBirth, Gender, Ethnicity, GENEPOS, GENEPOS_comb, ProcAge) %>%
  summarize(max_age = max(ContactAge, na.rm = TRUE),
            min_age = min(ContactAge, na.rm = TRUE),
            median_age = median(ContactAge, na.rm = TRUE))

### GROUPING AND MATCHING ------------------------------------------------------
# strict definition: controls also must have received genetic testing
df_cpt <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/cpts_genetic_all_pts.csv")

df_strict <- df %>%
  filter(MedicalRecordNumber %in% df_cpt$PAT_MRN_ID | GENEPOS_comb == "Y")

# strict case definition: exclude genetic individuals (cases) with VUS
df_rev <- readxl::read_excel("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/chartreview_2023-04-28.xlsx")
vec_rev <- df_rev %>% filter(is_genetic == FALSE) %>% pull(MedicalRecordNumber)

df_strict <- df_strict %>%
  filter(!MedicalRecordNumber %in% vec_rev)

# flag case and control explicitly
df_genes <- df_strict %>% mutate(status = ifelse(GENEPOS_comb == "N", 0,
                                                 ifelse(GENEPOS_comb == "Y", 1, NA)))

# map to HPO and propagate
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

# do matching
df_match <- matchit(status ~ median_age + Ethnicity + Gender, 
                    data = df_match, ratio = 1,
                    method = flag_match, distance = "glm")

df_match1 <- match.data(df_match)

df_match1 <- downsampleMatch(df_match1, df)

df_match1 <- df_match1 %>%
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

# subset original dataframe to this strict cohort
df <- df %>%
  filter(PatientId %in% df_match1$PatientId)

### DEMOGRAPHICS ---------------------------------------------------------------
tbl_person <- df_person %>%  
  filter(PatientId %in% df_match1$PatientId) %>%
  mutate(Gender = recode(Gender, 
                         "C0086582" = "Male",
                         "C0086287" = "Female")) %>%
  mutate(Ethnicity = recode(Ethnicity, 
                            "C1518424" = "Not Hispanic or Latino",
                            "C1549625" = "Unknown",
                            "C5441846" = "Hispanic or Latino",
                            "None" = "Unknown")) %>% 
  summary_factorlist(dependent = "GENEPOS_comb", 
                     explanatory = c("Gender", "Ethnicity", "ProcAge", 
                                     "min_age", "median_age", "max_age"),
                     p = TRUE, na_include = TRUE)  %>%
  knitr::kable("html") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) 

### VALIDITY ANALYSIS ----------------------------------------------------------
# repeat key analysis for this stricter case-control subset, cf. pipeline.R
# ...

### REPLICATION ANALYSIS -------------------------------------------------------
## save hypotheses (test, p-value, OR estimate) for cross-sectional phenotypes,
## longitudinal phenotypes, prescription data
write_csv(enrich1$data, "rep_cross-sectional_group2.csv")
write_csv(plong1f$plot$data, "rep_longitudinal_group2.csv")
write_csv(df_heatmap, "rep_medical_group2.csv")


