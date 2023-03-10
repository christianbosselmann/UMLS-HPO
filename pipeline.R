### UMLS-HPO ANALYSIS OF INDIVIDUALS WITH GENETIC EPILEPSY SYNDROMES -----------
## Main analysis pipeline
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

# data: list of MRNs per patient subgroup
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

## Group: genetic vs. non-genetic
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

## Group: SCN1A vs. CDKL5
df_match8 <- df_match %>%
  filter(status %in% c("scn1a", "cdkl5")) %>%
  mutate(status = recode(status, 
                         "cdkl5" = 0,
                         "scn1a" = 1)) 

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
                     explanatory = c("Gender", "Ethnicity", "ProcAge", 
                                     "min_age", "median_age", "max_age"),
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
  geom_text(aes(x = stats_followup$mean, label = "\nMean: 6·5 years", y = 400),
            colour = "black", angle = 90, size = 4)

## p2: flat violin (raincloud) plot of age at encounter
p2 <- ggplot(data = df, aes(x = GENEPOS_comb, y = ContactAge, fill = GENEPOS_comb)) +
  geom_flat_violin(position = position_nudge(x = .1, y = 0), alpha = .8) +
  guides(fill = "none", color = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Group") +
  ylab("Age at all encounters") +
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
## Group: genetic vs. non-genetic
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

## Group: SCN1A vs. CDKL5
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
plong1f <- longitudinalPlot(df_genes, df_match1, odds_plot = TRUE, filter = FALSE)

# set legend expression (for subscript)
str_legend <- expression(log['10']*(OR))

# set case terms
vec_caseterms <- c("HP:0000708", # behav. abnormality
                   "HP:0002719", # recurrent infections
                   "HP:0002311", # incoordination
                   "HP:0011968", # feeding difficulties
                   "HP:0002019", # constipation
                   "HP:0001944", # dehydration
                   "HP:0002415", # leukodystrophy
                   "HP:0000939", # osteoporosis
                   "HP:0002664", # neoplasm incl. neurofibromas
                   "HP:0012418", # hypoxemia
                   "HP:0002098", # respiratory distress
                   "HP:0000083") # renal insufficiency

# set control terms
vec_controlterms <- c("HP:0025234", # parasomnia
                      "HP:0001287", # meningitis
                      "HP:0001342", # cerebral hemorrhage
                      "HP:0003074", # Hyperglycemia
                      "HP:0003193", # Allergic rhinitis
                      "HP:0000388") # Otitis media

vec_longterms <- c(vec_caseterms, vec_controlterms)

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
dist_df <- lapply(lapply(sapply(dist_ls, unlist), "length<-", 
                         max(lengths(dist_ls))), as.list) %>%
  lapply(., function(x) {unlist(x); as.data.frame(t(x))}) %>%
  rbindlist()

# get matrix, fix NA and Inf
dist_mat <- as.matrix(dist_df)
dist_mat[which(is.na(dist_mat))] <- 0
dist_mat[dist_mat == Inf] <- 99

# calculate distance and clustering
dist_dist <- dist(dist_mat, method = "euclidean")
dist_clust <- hclust(dist_dist)

pheat1 <- plong1f$plot$data %>%
  ungroup() %>%
  # get each combination of bin and term
  mutate(bin = replace_na(bin, 99)) %>%
  mutate(bin = as.factor(bin)) %>%
  complete(bin, description) %>%
  group_by(description) %>%
  mutate(term = unique(term[!is.na(term)])) %>%
  ungroup() %>%
  # select terms of interest from discovery graph (plong) and filter by significance
  filter(term %in% vec_longterms) %>%
  mutate(odds = ifelse(pvalue > 0.05, NA, odds)) %>%
  # set facet groups
  mutate(fct_group = case_when(
    term %in% vec_caseterms ~ "Likely genetic",
    term %in% vec_controlterms ~ "Non-genetic")) %>%
  # set description factor levels
  mutate(description = as.factor(description)) %>%
  mutate(description = factor(description, levels = levels(description)[dist_clust$order])) %>%
  # plot
  ggplot(aes(x = as.factor(bin), 
             y = description,
             fill = log10(odds))) + 
  geom_tile(color = "black") +
  geom_text(aes(label = round(odds, 1)), size = 4) +
  theme_classic() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = NA),
        panel.spacing = unit(0.2, "cm")) +
  facet_grid(fct_group~., scales = "free_y", space = "free", 
             switch = "y", drop = FALSE, margins = FALSE) +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = NA),
        panel.spacing = unit(0.2, "cm")) +
  ylab("") +
  xlab("Age at encounter (years)") +
  scale_x_discrete(labels = c("0-2", "2-12", "12-18", ">18")) +
  scale_fill_gradient2(name = "OR",
                       labels = c("0", "", "1", "", "Inf"),
                       low = "#00b4fb",
                       mid = "#F5F5F5",
                       high = "#ff8422",
                       midpoint = 0, 
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
vec_asm <- recodeASM(df_asm)

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

## Descriptive stats
df_med <- df_med %>%
  group_by(PatientId) %>%
  mutate(n_unique_asm = n_distinct(MED_NAME)) %>% # number of unique ASMs per patient
  mutate(n_all_prescriptions = n_distinct(MonthsPrescription)) # number of all prescriptions per patient

## Group analysis: Heatmap (OR)
# define ASM groups
asm_map <- groupASM(asm_vec)

# take strict matched case-control set
df_heatmap <- df_match1 %>%
  # merge in patient ID and group label
  distinct(PatientId, group) %>%
  left_join(df_med, by = "PatientId") %>%
  na.omit %>%
  ## fixed bin width by ILAE age categories
  mutate(YearsPrescription = cut(YearsPrescription, breaks = c(0, 2, 12, 18, Inf))) %>%
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
  # # frequency filter
  # ungroup() %>%
  # mutate(freqY = Y/max(Y)) %>%
  # mutate(freqN = N/max(N)) %>%
  # filter(freqY > 0.01 & freqN > 0.01)
  # another frequency filter
  ungroup() %>%
  mutate(freqY = Y/max(Y)) %>%
  mutate(freqN = N/max(N)) %>%
  mutate(OR = ifelse(freqY > 0.01 & freqN > 0.01, OR, NA))

# plot
p_asm <- df_heatmap %>%
  mutate(MED_NAME = str_to_title(MED_NAME)) %>%
  # filter rows with no significant associations
  group_by(MED_NAME) %>%
  filter(any(!is.na(OR))) %>%
  ungroup() %>%
  # plot
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
  scale_x_discrete(labels = c("0-2", "2-12", "12-18", ">18")) +
  scale_fill_gradient2(name = "OR",
                       labels = c("0", "", "1", "", "Inf"),
                       low = "#00b4fb",
                       mid = "#F5F5F5",
                       high = "#ff8422",
                       midpoint = 0, 
                       na.value = "#F5F5F5",
                       limits = c(-2, 2),
                       oob = scales::oob_squish_any)

### NON-HPO CONCEPT ANALYSIS ---------------------------------------------------
# find out how many UMLS concepts do not match to HPO terms
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
  left_join(df_concepts[, c("PatientId", "ConceptID", "ConceptDesc")], 
            by = "PatientId") %>%
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

# manual annotation (labels missing in MRCONSO.RFF)
df_conceptmatch[df_conceptmatch$ConceptID == "C0478107", ]$ConceptDesc <- NA
df_conceptmatch[df_conceptmatch$ConceptID == "C0476431", ]$ConceptDesc <- NA
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
# find what drives the increase in encounters for likely genetic
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
# note: compare to non-genetic cohort
# sorted by p-value or frequency in the positive transition group

# manual annotation for plot and recoding of synonymous concepts
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
# manual annotation for plot and recoding of synonymous concepts
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

### MORTALITY ------------------------------------------------------------------
## data
df_deaths <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/Deathdates.csv",
                      col_names = c("MedicalRecordNumber", "Status", "DateOfDeath"))

# map to PatientId
mrn_map <- df_raw %>% distinct(PatientId, MedicalRecordNumber)

df_deaths <- df_deaths %>%
  left_join(mrn_map, by = "MedicalRecordNumber")

# get age at death
dob_map <- df_raw %>% distinct(PatientId, DateOfBirth)

df_deaths <- df_deaths %>%
  left_join(dob_map, by = "PatientId")

df_deaths <- df_deaths %>%
  mutate(DateOfBirth = as.Date(DateOfBirth, format = "%m/%d/%y")) %>%
  mutate(DateOfDeath = as.Date(DateOfDeath, format = "%m/%d/%Y")) %>%
  mutate(AgeAtDeath = DateOfDeath - DateOfBirth)

# add death dates to matched cohort; only keep distinct patient-death pairs
df_death1 <- df_match1 %>%
  left_join(df_deaths[ ,c("PatientId", "AgeAtDeath")], by = "PatientId")

# change ContactAge to days for survival analysis
df_death1 <- df_death1 %>%
  mutate(ContactAge = ContactAge*365)

# if a patient is deceased, set their age at last contact for the survival analysis
df_death1 <- df_death1 %>%
  mutate(isDead = if_else(!is.na(AgeAtDeath), 1, 0))

df_surv <- df_death1 %>%
  filter(isDead == 1) %>%
  group_by(PatientId) %>%
  slice_max(order_by = ContactAge, n = 1) %>%
  mutate(survivalFlag = 1)

df_surv <- left_join(df_death1, df_surv) %>%
  mutate(survivalFlag = if_else(is.na(survivalFlag),0,1))

# reduce to single encounters
df_surv <- df_surv %>%
  distinct(PatientId, group, ContactAge, survivalFlag)

# preprocessing for suvival analysis
df_surv <- df_surv %>%
  mutate(group = as.integer(group)) %>%
  mutate(ContactAge = as.numeric(ContactAge)/365)

# for each patient, keep only their last known status
df_surv <- df_surv %>%
  group_by(PatientId) %>%
  slice_max(order_by = ContactAge, n = 1)

km_fit <- survfit(Surv(ContactAge, survivalFlag) ~ group, data = df_surv)

p_surv <- survminer::ggsurvplot(km_fit, data = df_surv,
                                pval = TRUE,
                                risk.table = "nrisk_cumcensor",
                                pval.coord = c(2, 0.55),
                                ggtheme = theme_classic(),
                                xlab = c("Age (years)"),
                                xlim = c(0, 25), ylim = c(0.50, 1.0),
                                legend = c(0.8, 0.15),
                                legend.labs = c("Non-genetic", "Likely genetic"),
                                palette = "Dark2") 

### INPATIENT / OUTPATIENT STATS -----------------------------------------------
## data: use prescription and admission data to find inpatient/outpatient encounters
df_stays <- df_med %>%
  distinct(PatientId, AgePrescription, ORD_MODE_DESC)

df_nstays <- df_stays %>%
  mutate(ORD_MODE_DESC = factor(ORD_MODE_DESC, levels = c("INPATIENT", "OUTPATIENT"))) %>%
  group_by(PatientId, .drop = FALSE) %>%
  count(ORD_MODE_DESC)

# find PatientIds who had multiple inpatient encounters
df_inpatient <- df_nstays %>%
  filter(ORD_MODE_DESC == "INPATIENT" & n > 1) %>%
  mutate(hasInpatient = TRUE)

### UTILIZATION ANALYSIS -------------------------------------------------------
## data: list of patient specialist encounters from Alina Ivaniuk, 2023-03-07
# TODO: update input data. current table omits rows as max file length has been reached.
df_util <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/encounter-info-epilepsy.csv")

# subset by case-control cohort
df_util <- df_util %>%
  filter(PatientID %in% df_match1$PatientId)

# get date of birth from df_lookup and calculate age at encounter
df_util <- df_util %>%
  rename(PatientId = PatientID) %>%
  left_join(df_lookup[ ,c("PatientId", "DateOfBirth")], by = "PatientId") %>%
  mutate(DateOfBirth = lubridate::mdy(DateOfBirth)) %>%
  mutate(ENC_DT = lubridate::as_date(ENC_DT)) %>%
  mutate(AgeAtEncounter = difftime(ENC_DT, DateOfBirth, units = "days")) %>%
  mutate(MonthsEncounter = as.numeric(round(AgeAtEncounter/30, 0))) %>%
  mutate(YearsEncounter = as.numeric(AgeAtEncounter/365.2425))

# get group label
map_match <- df_match1 %>%
  distinct(PatientId, group)

df_util <- df_util %>%
  left_join(map_match)

# get count of unique specialty descriptions per patient
# TODO: use number of unique specialties seen for PheIndex score instead of current assumption
stats_util <- df_util %>%
  group_by(PatientId) %>%
  distinct(PatientId, group, SPECIALTY_DESC) %>%
  count()

# idea: we could look at ENC_TYPE_DESC between groups and during 2019-2022 vs before. Telehealth during COVID?

### ER VISITS -----------------------------------------------------------------
## data: ER admissions for all patients
df_er <- readxl::read_excel("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/ER_Visits.xlsx")

# get group label
df_er <- df_er %>%
  rename(PatientId = PatientID) %>%
  left_join(map_match) %>%
  na.omit

# count admissions for each patient
p_er <- df_er %>%
  add_count(PatientId, AdmissionType) %>%
  summarize(pval = t.test(n ~ group)$p.value)

stats_er <- df_er %>%
  group_by(group) %>%
  add_count(PatientId, AdmissionType) %>%
  summarize(mean = mean(n), median = median(n),
            sd = sd(n),
            min = min(n), max = max(n),
            n = n_distinct(PatientId))

### PHEINDEX SCORING -----------------------------------------------------------
## ref: https://www.medrxiv.org/content/10.1101/2023.01.27.23285056v1.full.pdf
## data: list of PatientIds, ContactAges, group membership and concepts
df_pheindex <- df_concepts[df_concepts$PatientId %in% df_match1$PatientId, ] %>%
  distinct(PatientId, ConceptID) %>%
  mutate(score = 0)

## C0: list of dataframes for scoring
ls_pheindex <- list()

## C1: Prolonged stay in the neonatal intensive care unit
# data not available

## C2: Prolonged or multiple hospitalizations after discharged from birth
# hasInpatient flag from df_inpatient
ls_pheindex[[2]]  <- df_inpatient %>%
  distinct(PatientId, hasInpatient) %>%
  filter(hasInpatient == TRUE) %>%
  mutate(score = 3)

## C3: Visits or consults with multiple specialists other than general pediatricians.
# True for all patients (due to cohort definition)
ls_pheindex[[3]]  <- df_pheindex %>%
  distinct(PatientId) %>%
  mutate(score = 3)

## C4: Multiple emergency room (ER) visits.
# NA

## C5: Feeding support (Gastrostomy tube).
ls_pheindex[[5]] <- rbind(
  # df_pheindex[df_pheindex$ConceptID == "C0699815", ], # Feeding difficulties and mismanagement
  # df_pheindex[df_pheindex$ConceptID == "C0159023", ], # Feeding problems in newborn
  # df_pheindex[df_pheindex$ConceptID == "C5539211", ], # Other feeding difficulties
  df_pheindex[df_pheindex$ConceptID == "C0260683", ], # Gastrostomy status
  # df_pheindex[df_pheindex$ConceptID == "C0270273", ], # Slow feeding in newborn
  # df_pheindex[df_pheindex$ConceptID == "C5539209", ], # Feeding difficulties, unspecified
  # df_pheindex[df_pheindex$ConceptID == "C0478153", ], # Other symptoms and signs concerning food and fluid intake
  df_pheindex[df_pheindex$ConceptID == "C0260761", ] # Encounter for attention to gastrostomy
) %>%
  distinct(PatientId, score) %>%
  mutate(score = 2)

## C6: Respiratory support (tracheostomy and mechanical ventilation outside of surgery).
ls_pheindex[[6]] <- rbind(
  # df_pheindex[df_pheindex$ConceptID == "C0348712", ], # Other disorders of lung
  # df_pheindex[df_pheindex$ConceptID == "C0431510", ], # Other anomalies of larynx, trachea, and bronchus
  # df_pheindex[df_pheindex$ConceptID == "C0029601", ], # Other respiratory anomalies
  df_pheindex[df_pheindex$ConceptID == "C0260682", ], # Tracheostomy status
  # df_pheindex[df_pheindex$ConceptID == "C0748355", ], # Acute respiratory distress
  # df_pheindex[df_pheindex$ConceptID == "C0456017", ], # Chronic respiratory disease in perinatal period
  df_pheindex[df_pheindex$ConceptID == "C2911575", ] # Dependence on respirator [ventilator] status
  # df_pheindex[df_pheindex$ConceptID == "C2977073", ] # Respiratory failure, unspecified
) %>%
  distinct(PatientId, score) %>%
  mutate(score = 2)

## C7: Imaging.
# NA

## C8: Genetic diagnostic tests.
# True for all likely genetic patients, false for all non-genetic patients.
ls_pheindex[[8]] <- df_match1 %>%
  distinct(PatientId, group) %>%
  filter(group == TRUE) %>%
  mutate(score = 1)

## C9: Metabolic diagnostic tests
ls_pheindex[[9]] <- rbind(
  df_pheindex[df_pheindex$ConceptID == "C0494356", ], # Hypo-osmolality and hyponatremia
  df_pheindex[df_pheindex$ConceptID == "C0020645", ], # Hyposmolality and/or hyponatremia
  df_pheindex[df_pheindex$ConceptID == "C0029481", ] # Other abnormal blood chemistry
) %>%
  distinct(PatientId, score) %>%
  mutate(score = 1)

## C10: In-hospital death
# True for all deaths in survival analysis, false otherwise
ls_pheindex[[10]] <- df_surv %>%
  filter(survivalFlag == 1) %>%
  distinct(PatientId) %>%
  mutate(score = 3)

## C11: Developmental delay.
ls_pheindex[[11]] <- rbind(
  df_pheindex[df_pheindex$ConceptID == "C0424605", ], # Developmental delay
  df_pheindex[df_pheindex$ConceptID == "C0878706", ], # Lack of normal physiological development, unspecified
  df_pheindex[df_pheindex$ConceptID == "C0476241", ], # Delayed developmental milestones
  df_pheindex[df_pheindex$ConceptID == "C0878753", ], # Unspecified lack of expected normal physiological development in childhood
  df_pheindex[df_pheindex$ConceptID == "C0154633", ], # Other developmental speech or language disorder
  df_pheindex[df_pheindex$ConceptID == "C2830458", ], # Delayed milestone in childhood
  df_pheindex[df_pheindex$ConceptID == "C0011757", ], # Developmental Coordination Disorder
  df_pheindex[df_pheindex$ConceptID == "C0349324", ], # Other developmental disorders of speech and language
  df_pheindex[df_pheindex$ConceptID == "C0236826", ], # Developmental expressive language disorder
  df_pheindex[df_pheindex$ConceptID == "C3161331", ] # Unspecified intellectual disabilities
) %>%
  distinct(PatientId, score) %>%
  mutate(score = 1)

## C12: Diagnosis codes corresponding to metabolic diseases with ≥ 2 encounters
ls_pheindex[[12]] <- rbind(
  df_pheindex[df_pheindex$ConceptID == "C0025517", ], # Metabolic Diseases
  df_pheindex[df_pheindex$ConceptID == "C0268641", ] # Amino acid transport disorder
) %>%
  distinct(PatientId, score) %>%
  mutate(score = 3)

## C13: Heart surgeries
ls_pheindex[[13]] <- rbind(
  df_pheindex[df_pheindex$ConceptID == "C2921289", ], # Personal history of (corrected) congenital malformations of heart and circulatory system 
  df_pheindex[df_pheindex$ConceptID == "C0477999", ] # Other specified congenital malformations of heart
) %>%
  distinct(PatientId, score) %>%
  mutate(score = 3)

## get group membership for each patient, calculate pheindex scores
df_pheindex <- lapply(ls_pheindex, function(x){x <- x[,c("PatientId", "score")]}) %>%
  rbindlist(idcol = "id")

map_match <- df_match1 %>%
  distinct(PatientId, group)

df_pheindex <- df_pheindex %>%
  group_by(PatientId) %>%
  summarize(score = sum(score)) %>%
  left_join(map_match) %>%
  na.omit

stats_pheindex <- df_pheindex %>%
  group_by(group) %>%
  summarize(mean = mean(score),
            sd = sd(score),
            min = min(score),
            max = max(score),
            iqr = IQR(score))

stats_p_pheindex <- df_pheindex %>%
  summarize(pval = t.test(score ~ group)$p.value)

## recode
df_pheindex <- df_pheindex %>%
  mutate(group = as.integer(group))

df_pheindex$label <- NA
df_pheindex[df_pheindex$group == 0, ]$label <- "Non-genetic"
df_pheindex[df_pheindex$group == 1, ]$label <- "Likely genetic"

## visualization: PhenIndex by group
p_pheindex_violin <-df_pheindex %>%
  mutate(label = as.factor(label)) %>%
  ggplot(aes(y = score, x = label, fill = label)) +
  geom_flat_violin(position = position_nudge(x = .1, y = 0), alpha = .8) +
  guides(fill = "none", color = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Group") +
  ylab("PheIndex score") +
  theme_classic() + 
  coord_cartesian(xlim = c(1.5, 2)) +
  geom_boxplot(width = .1, show.legend = FALSE, outlier.shape = NA, alpha = 0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..), 
                             comparisons = list(c("Non-genetic", "Likely genetic")),
                             label.x = 1.5, label.y = c(15))

## visualization: PhenIndex categories by group (id)
df_phenindex <- lapply(ls_pheindex, function(x){x <- x[,c("PatientId", "score")]}) %>%
  rbindlist(idcol = "id") %>%
  group_by(PatientId, id) %>%
  summarize(score = sum(score)) %>%
  left_join(map_match) %>%
  na.omit

p_pheindex_bar <- df_phenindex %>%
  ggplot(aes(x = as.factor(id), fill = group)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab("PheIndex criteria (#)") +
  ylab("Individuals (n)") +
  guides(fill = "none", color = "none") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2", labels = c("Non-genetic", "Likely genetic")) +
  theme_classic() +
  coord_cartesian(expand = FALSE) 

# PheIndex criteria pairwise similarity heatmap
df_phenindex <- df_phenindex %>%
  mutate(id = as.factor(id))

ls_phenindex <- split(df_phenindex, df_phenindex$id)
ls_phenindex <- lapply(ls_phenindex, function(x){x <- x$PatientId})

df_prop <- ls_phenindex %>% 
  map_dfr(~ .x %>% as_tibble(), .id = "name") 

df_prop <- reshape2::dcast(df_prop, name ~ value, length)

df_prop <- df_prop[,-1]
df_prop[df_prop > 0] <- 1 

mat_pheno <- proxy::simil(x = df_prop,
                          method = "cosine") # "jaccard", "euclidean", "cosine"

mat_pheno <- as.matrix(mat_pheno)
diag(mat_pheno) <- 1

colnames(mat_pheno) <- levels(df_phenindex$id)
rownames(mat_pheno) <- levels(df_phenindex$id)

coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)

# plot
p_pheindex_heatmap <- pheatmap(mat_pheno,
                               display_numbers = TRUE,
                               symm = TRUE,
                               col = coul,
                               fontsize_number = 6,
                               xlab = "PheIndex criteria (#)",
                               ylab = "PheIndex criteria (#)")

p_pheindex_heatmap <- ggplotify::as.ggplot(p_pheindex_heatmap) 

# # CX: Preterm
# df_pheindex[df_pheindex$ConceptID == "C2909946", ] # Preterm newborn, unspecified weeks of gestation
# df_pheindex[df_pheindex$ConceptID == "C0029713", ] # Other preterm infants
# df_pheindex[df_pheindex$ConceptID == "C3264533", ] # Extreme immaturity of newborn

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
                                 edge.arrow.size = .5,
                                 layout = layout_as_tree))) # Reingold-Tilford 

pqg <- pqg + 
  theme(plot.margin = unit(c(-50, -20, -50, -50), "pt"))

### GENERATE REPORT ------------------------------------------------------------
# set global option of decimal point to midline for Lancet
options(OutDec = "·")

## Figure 1: Descriptive statistics of the study cohort.
Fig1 <- cowplot::plot_grid(p1,
                           p2 + scale_x_discrete(labels=c("Non-genetic", "Likely genetic")),
                           p3 + scale_x_discrete(labels=c("Non-genetic", "Likely genetic")),
                           p5 + theme(legend.position = "none"),
                           p4 + theme(legend.position = "none"),
                           pt,
                           nrow = 2, labels = "AUTO", align = "none")

pdf(file = "Fig1.pdf",
    width = 12,
    height = 8)

Fig1

dev.off()

## Figure 2: Genetic vs. Non-Genetic
Fig2 <- cowplot::plot_grid(pqq, 
                           p_forest_nonhpo,
                           enrich1$forest,
                           pqg,
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
                           nrow = 2, labels = "AUTO", align = "none")

pdf(file = "Fig2.pdf",
    width = 12,
    height = 8,
    encoding = 'CP1253.enc') # to draw Lambda on panel A

Fig2

dev.off()

## Figure 3: Nicer heatmaps of longitudinal phenotypes and prescription patterns
Fig3 <- cowplot::plot_grid(pheat1 + theme(legend.position = "none"), 
                           p_asm,
                           nrow = 1, labels = "AUTO", align = "none")

pdf(file = "Fig3.pdf",
    width = 12,
    height = 8)

Fig3

dev.off()

## Figure S1: Kaplan-Meier plot
pdf(file = "FigS1.pdf",
    width = 8,
    height = 6)

p_surv

dev.off()

## Figure S2: PheIndex plots
FigS2 <- cowplot::plot_grid(p_pheindex_violin, p_pheindex_bar, p_pheindex_heatmap,
                            nrow = 1, labels = "AUTO", align = "hv")

pdf(file = "FigS2.pdf",
    width = 12,
    height = 4)

FigS2

dev.off()

### CHART REVIEW --------------------------------------------------------------
# pull patient MRNs for manual chart review
df_match1 %>%
  filter(term == "HP:0000083") %>% # renal insufficiency
  distinct(PatientId) %>%
  left_join(mrn_map) %>%
  write_csv("/Users/cbosselmann/Desktop/hp_0000083.csv")

df_match1 %>%
  filter(term == "HP:0004383") %>% # hypoplastic left heart syndrome
  distinct(PatientId) %>%
  left_join(mrn_map) %>%
  write_csv("/Users/cbosselmann/Desktop/hp_0004383.csv")

df_match1 %>%
  filter(term == "HP:0000028") %>% # cryptorchidism
  distinct(PatientId) %>%
  left_join(mrn_map) %>%
  write_csv("/Users/cbosselmann/Desktop/hp_0000028.csv")

df_match1 %>%
  filter(term == "HP:0000047") %>% # hypospadia
  distinct(PatientId) %>%
  left_join(mrn_map) %>%
  write_csv("/Users/cbosselmann/Desktop/hp_0000047.csv")

df_match1 %>%
  filter(term == "HP:0000939") %>% # osteoporosis
  distinct(PatientId) %>%
  left_join(mrn_map) %>%
  write_csv("/Users/cbosselmann/Desktop/hp_0000939.csv")

### MISC ----------------------------------------------------------------------
## get the number of SCN1A patients in the cohort and their mean follow-up
tmp_scn1a <- df_scn1a %>%
  rename(MedicalRecordNumber = PAT_MRN_ID) %>%
  left_join(mrn_map)

tmp_scn1a_fu <- p1$data %>%
  filter(PatientId %in% tmp_scn1a$PatientId) %>%
  mutate(dur = upper-lower) %>%
  summarize(mean = mean(dur), median = median(dur),
            sd = sd(dur), min = min(dur), max = max(dur),
            iqr = IQR(dur))

## DL analysis 08/03/2023
# wants to see group differences between non-mapped terms in a list provided by Mark
tmp_xl <- readxl::read_excel("~/Desktop/EpilepsyHPO Mappings.xlsx", sheet = 2)
tmp_xl$ConceptID

# find non-hpo analysis output for these concepts
df_conceptmatch %>%
  filter(ConceptID %in% tmp_xl$ConceptID) %>%
  left_join(tmp_xl[ ,c("ConceptID", "ConceptDescription")]) %>%
  write_csv("~/Desktop/nonhpo_2023-03-08.csv")






