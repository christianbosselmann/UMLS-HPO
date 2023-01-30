# surgery data

# pkg
library(librarian)
librarian::shelf(tidyverse,
                 ggrepel,
                 ontologyIndex,
                 data.table,
                 scales)

# helper fn
source("func.R")

# load ontology
ont_hpo <- get_ontology("hp.obo.txt", 
                        propagate_relationships = "is_a", 
                        extract_tags = "everything")

# data: diagnosis from all encounters for two years before day of most recent surgery
df <- read_csv("~/Desktop/CCF/EMR cohort study/Surgery cohort/data/SurgeryData.csv") %>%
  select(Surgery, ConceptID)

# manual map
hpo_map <- lapply(ont_hpo$xref, function(x){
  x <- x[x %like% "UMLS:"]
  x <- sub('.*\\:', '', x)
}) 

hpo_map <- enframe(hpo_map) %>%
  unnest(value) %>%
  rename(ConceptID = value)

desc_map <- tibble(term = ont_hpo$id,
                   description = ont_hpo$name)

### CLEAN DATA ----------------------------------------------------------------
# fix broken labels
df <- df %>%
  dplyr::mutate(Surgery = recode(Surgery,
                                 "Temporal Lobectomy" = "Temporal lobectomy",
                                 "Cental resection" = "Central resection",
                                 "other" = "Other",
                                 "SD/grids" = "SD grids", 
                                 "VNS lead revision" = "VNS revision",
                                 "VNS lead removal" = "VNS removal",
                                 "VNS Removal" = "VNS removal",
                                 "VNS implantation" = "VNS insertion",
                                 "Depth electrodes" = "SEEG"))

# list for HPO analysis
ls <- split(df, df$Surgery, df$ConceptID)

### DESCRIPTIVE STATS ---------------------------------------------------------
# absolute number and unique concepts per group
df_stats <- tibble(surgery = names(ls),
                   terms_absolute = sapply(ls, function(x) {nrow(x)}),
                   terms_unique = sapply(ls, function(x) {x$ConceptID %>% unique %>% length}))

# custom fill groups
df_stats <- df_stats %>% mutate(surgery = ifelse(surgery %in% c("Temporal lobectomy"), "Temporal lobectomy",
                                        ifelse(surgery %in% c("Frontal lobectomy", "Occipital lobectomy", "Parietal lobectomy"), "Extratemporal lobectomy",
                                        ifelse(surgery %in% c("Insular resection", "Central resection", "Multilobar resection"), "Resection, other",
                                        ifelse(surgery %in% c("Laser ablation"), "Laser ablation",
                                        ifelse(surgery %in% c("Hemispherectomy"), "Hemispherectomy",
                                        ifelse(surgery %in% c("Callosotomy"), "Callosotomy",
                                        ifelse(surgery %in% c("VNS", "VNS insertion", "VNS removal", "VNS revision"), "VNS",
                                        ifelse(surgery %in% c("DBS", "Generator implantation", "Explantation of NeuroPace", "NeuroPace"), "Neurostimulation, other", 
                                        ifelse(surgery %in% c("SD grids", "Other - Plates", "SEEG"), "SD/SEEG",
                                        ifelse(surgery %in% c("Other", "Other - Shunt", "Other-Ventriculostomy"), "Other", "Other"
                                        ))))))))))) %>%
  group_by(surgery) %>%
  dplyr::summarise(n_abs = sum(terms_absolute), n_uniq = sum(terms_unique))

# plot
ggplot(df_stats, aes(x = reorder(surgery, -n_abs), y = n_abs, fill = surgery)) + 
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, max(df_stats$n_abs)), expand = FALSE) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(label = comma, "Concept IDs (absolute)") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplot(df_stats, aes(x = reorder(surgery, -n_uniq), y = n_uniq, fill = surgery)) + 
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, max(df_stats$n_uniq)), expand = FALSE) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(label = comma, "Concept IDs (unique)") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# note: abs/unique can be a surrogate for phenotypic complexity, i.e. how well characterized

### FIGURE 1: TLE vs ExTLE ----------------------------------------------------
# join
df_map <- left_join(df, hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  na.omit

# different filter before defining groups; first group is positive
df_map <- df_map %>% 
  filter(Surgery %in% c("Temporal lobectomy",
                        "Frontal lobectomy",
                        "Central resection",
                        "Cental resection", 
                        "Extratemporal resection",
                        "Insular resection",
                        "Occipital lobectomy",
                        "Parietal lobectomy"))

# define positive group
df_map <- df_map %>%
  mutate(group = Surgery == "Temporal lobectomy")

# run fn
res <- enrichmentPlot(df_map, ont_hpo)
# note: coord system can be forced:
# res$plot + coord_fixed(xlim = c(0, 0.25), ylim = c(0, 0.25))

### FIGURE 2: VNS vs Resection -------------------------------------------------
# join
df_map <- left_join(df, hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  na.omit

# define groups
df_map <- df_map %>% 
  filter(Surgery %in% c("Temporal lobectomy",
                        "Frontal lobectomy",
                        "Central resection",
                        "Cental resection", 
                        "Extratemporal resection",
                        "Insular resection",
                        "Occipital lobectomy",
                        "Parietal lobectomy",
                        "VNS insertion",
                        "VNS removal",
                        "VNS Removal",
                        "VNS revision",
                        "VNS lead removal",
                        "VNS lead revision",
                        "VNS implantation"))

df_map <- df_map %>%
  mutate(group = Surgery %in% c("VNS insertion",
                                "VNS removal",
                                "VNS Removal",
                                "VNS revision",
                                "VNS lead removal",
                                "VNS lead revision",
                                "VNS implantation"))

# run fn
res <- enrichmentPlot(df_map, ont_hpo)

### FIGURE 3: Hemispherectomy vs Resection--------------------------------------
# join
df_map <- left_join(df, hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  na.omit

# define groups
df_map <- df_map %>%
  filter(Surgery %in% c("Temporal lobectomy",
                        "Frontal lobectomy",
                        "Central resection",
                        "Cental resection",
                        "Extratemporal resection",
                        "Insular resection",
                        "Occipital lobectomy",
                        "Parietal lobectomy",
                        "Hemispherectomy"))

df_map <- df_map %>%
  mutate(group = Surgery %in% c("Hemispherectomy"))

# run fn
res <- enrichmentPlot(df_map, ont_hpo)

### FIGURE 4: Callosotomy vs Resection -----------------------------------------
# join
df_map <- left_join(df, hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  na.omit

# define groups
df_map <- df_map %>%
  filter(Surgery %in% c("Temporal lobectomy",
                        "Frontal lobectomy",
                        "Central resection",
                        "Cental resection",
                        "Extratemporal resection",
                        "Insular resection",
                        "Occipital lobectomy",
                        "Parietal lobectomy",
                        "Corpus callosotomy"))

df_map <- df_map %>%
  mutate(group = Surgery %in% c("Corpus callosotomy"))

# run fn
res <- enrichmentPlot(df_map, ont_hpo)

### FIGURE 5: Multilobar vs TLE ------------------------------------------------
# join
df_map <- left_join(df, hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  na.omit

# define groups
df_map <- df_map %>%
  filter(Surgery %in% c("Temporal lobectomy",
                        "Multilobar resection"))

df_map <- df_map %>%
  mutate(group = Surgery %in% c("Multilobar resection"))

# run fn
res <- enrichmentPlot(df_map, ont_hpo)

### FIGURE 6: VNS in vs VNS out ------------------------------------------------
# join
df_map <- left_join(df, hpo_map, by = "ConceptID") %>%
  rename(term = name) %>%
  na.omit

# define groups
df_map <- df_map %>%
  filter(Surgery %in% c("VNS insertion",
                        "VNS removal"))

df_map <- df_map %>%
  mutate(group = Surgery %in% c("VNS removal"))

# run fn
res <- enrichmentPlot(df_map, ont_hpo)
