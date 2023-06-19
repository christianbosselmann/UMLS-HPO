### UMLS-HPO ANALYSIS OF INDIVIDUALS WITH GENETIC EPILEPSY SYNDROMES -----------
## Validity analysis
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
  filter(MedicalRecordNumber %in% df_cpt$PAT_MRN_ID)

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

### VALIDITY ANALYSIS ----------------------------------------------------------
# repeat key analysis for this stricter case-control subset, cf. pipeline.R

## p1: flag plot of encounters over age
p1 <- df %>%
  # optional: subset to matched cohort
  filter(PatientId %in% df_match1$PatientId) %>%
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
  # optional: subset to matched cohort
  filter(PatientId %in% df_match1$PatientId) %>%
  mutate(dur = upper-lower) %>%
  summarize(mean = mean(dur), median = median(dur),
            sd = sd(dur), min = min(dur), max = max(dur),
            iqr = IQR(dur))

# can also compare start of follow-up between groups
p1$data %>%
  filter(PatientId %in% df_match1$PatientId) %>%
  left_join(df_match1[ ,c("PatientId", "group")]) %>%
  mutate(dur = upper-lower) %>%
  summarize(pval = t.test(lower ~ group)$p.value)

p1$data %>%
  filter(PatientId %in% df_match1$PatientId) %>%
  left_join(df_match1[ ,c("PatientId", "group")]) %>%
  mutate(dur = upper-lower) %>%
  group_by(group) %>%
  summarize(mean = mean(lower), median = median(lower),
            sd = sd(lower), min = min(lower), max = max(lower),
            iqr = IQR(lower))

# add mean age at follow-up back to plot
p1 <- p1 + 
  geom_vline(xintercept = stats_followup$mean, linetype = "dashed") +
  geom_text(aes(x = stats_followup$mean, label = "\nMean: 8·18 years", y = 100),
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
enrich1$plot$data[enrich1$plot$data$description == "Simple febrile seizure", ]$expcat_text <- "Simple febrile seizure"
enrich1$plot$data[enrich1$plot$data$description == "Hyperactivity", ]$expcat_text <- "Hyperactivity"

enrich1$plot <- enrich1$plot +
  coord_fixed(xlim = c(0, .095), ylim = c(0, .095)) +
  scale_x_continuous( breaks=pretty_breaks()) +
  scale_y_continuous( breaks=pretty_breaks()) +
  ggtitle("Genetic vs. Non-Genetic") +
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
  ungroup() %>%
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

## subanalysis: do likely genetic patients receive more ASMs?
df_med_sub <- df_med %>%
  filter(PatientId %in% df_match1$PatientId) %>%
  left_join(df_match1[ ,c("PatientId", "group")]) %>%
  distinct(PatientId, AgePrescription, MED_NAME, group)

stats_med_sub <- df_med_sub %>%
  group_by(PatientId, group) %>%
  summarize(asm = n_distinct(MED_NAME)) %>%
  group_by(group) %>%
  summarize(mean = mean(asm), sd = sd(asm), min = min(asm), max = max(asm))

p_med_sub <- df_med_sub %>%
  group_by(PatientId, group) %>%
  summarize(asm = n_distinct(MED_NAME)) %>%
  ungroup() %>% 
  summarize(pval = t.test(asm ~ group)$p.value)

## subanalysis: Do likely genetic patients receive more rescue medication?
df_med_resc <- df_med_sub %>%
  mutate(isRescue = if_else(MED_NAME %in% c("lorazepam", "clonazepam", "clobazam", "diazepam", "midazolam"), T, F))

stats_med_resc <- df_med_resc %>%
  filter(isRescue) %>%
  group_by(PatientId, group) %>%
  summarize(asm = n_distinct(AgePrescription)) %>%
  group_by(group) %>%
  summarize(mean = mean(asm), sd = sd(asm), min = min(asm), max = max(asm))

p_med_resc <- df_med_resc %>%
  filter(isRescue) %>%
  group_by(PatientId, group) %>%
  summarize(asm = n_distinct(AgePrescription)) %>%
  ungroup() %>% 
  summarize(pval = t.test(asm ~ group)$p.value)

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
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05/length(unique(df_match1$ConceptID))), linetype = "dashed", col = "red") 

# keep significant associations: FDR
df_conceptmatch <- df_conceptmatch %>%
  filter(P < 0.05)

# get descriptions; note: to discuss in text
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
  filter(status %in% c("0", "1")) %>%
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
df_trans[df_trans$ConceptID == "C2911575", ]$ConceptDesc <- "Dependence on respirator [ventilator] status"
df_trans[df_trans$ConceptID == "C0477972", ]$ConceptDesc <- "Congenital malformations of brain"

# forest plot: genetic group
pt1 <- df_trans %>%
  filter(status == "1") %>%
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

df_trans[df_trans$ConceptID == "C0477972", ]$ConceptDesc <- "Other specified congenital malformations of brain"
df_trans[df_trans$ConceptID == "C2875146", ]$ConceptDesc <- "Migraine without aura, not intractable"
df_trans[df_trans$ConceptID == "C0155886", ]$ConceptDesc <- "Asthma"
df_trans[df_trans$ConceptID == "C2910447", ]$ConceptDesc <- "Medical examination w/o abnormal findings"
df_trans[df_trans$ConceptID == "C0154714", ]$ConceptDesc <- "Focal epilepsy, non-intractable"

df_trans[df_trans$ConceptID == "C2875148", ]$ConceptDesc <- "Migraine without aura, intractable"
df_trans[df_trans$ConceptID == "C0028768", ]$ConceptDesc <- "Obsessive-Compulsive Disorder"

pt2 <- df_trans %>%
  filter(status == "0") %>%
  na.omit %>%
  slice_max(order_by = Y, n = 4, with_ties = FALSE) %>%
  ggplot(aes(y = reorder(ConceptDesc, OR))) +
  geom_linerange(aes(xmin = CI1, xmax = CI2)) +
  geom_point(aes(x = OR), shape = 15, size = 3, color = "#1b9e77") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(trans = 'log10',
                     limits = c(1, 35),
                     oob = scales::oob_squish) +
  scale_y_discrete(labels = label_wrap(30)) +
  expand_limits(x = 1) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("Odds ratio (95% CI)")

pt <- cowplot::plot_grid(pt1 + theme(axis.title.x = element_blank()), 
                         pt2, 
                         ncol = 1, align = "v")

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

# stats: are likely genetic patients more likely to be seen inpatient?
map_match <- df_match1 %>%
  distinct(PatientId, group)

p_inpatient <- df_inpatient %>%
  left_join(map_match) %>%
  na.omit %>%
  ungroup() %>%
  summarize(pval = t.test(n ~ group)$p.value)

stats_inpatient <- df_inpatient %>%
  left_join(map_match) %>%
  na.omit %>%
  group_by(group) %>%
  summarize(mean = mean(n), median = median(n),
            sd = sd(n),
            min = min(n), max = max(n),
            n = n_distinct(PatientId))

### ER VISITS -----------------------------------------------------------------
## data: ER admissions for all patients
df_er <- readxl::read_excel("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/ER_Visits.xlsx") 

# get group label
df_er <- df_er %>%
  rename(PatientId = PatientID) %>%
  left_join(map_match) %>%
  na.omit

# subset to strict cohort criteria
df_er <- df_er %>%
  filter(PatientId %in% df$PatientId)

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

df_er <- df_er %>%
  group_by(group) %>%
  add_count(PatientId, AdmissionType) %>%
  distinct(PatientId, n)

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


