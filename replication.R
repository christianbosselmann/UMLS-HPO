### UMLS-HPO ANALYSIS OF INDIVIDUALS WITH GENETIC EPILEPSY SYNDROMES -----------
## Replication
##
## Author: Christian Bosselmann, MD
##
## Date Created: 2023-06-19
##
## Copyright (c) Christian Bosselmann, 2023
## Email: bosselc@ccf.org
##
### ----------------------------------------------------------------------------
# this script analyzes the reproducibility of associations and their effect size between the two stages

### DATA -----------------------------------------------------------------------
## hypotheses for cross-sectional clinical characteristics, longitudinal
## clinical characteristics, and prescription patterns
ls_cs <- list(read_csv("rep_cross-sectional_group1.csv"),
              read_csv("rep_cross-sectional_group2.csv")) %>%
  rbindlist(idcol = "group")

ls_lg <- list(read_csv("rep_longitudinal_group1.csv"),
              read_csv("rep_longitudinal_group2.csv")) %>%
  rbindlist(idcol = "group")

ls_md <- list(read_csv("rep_medical_group1.csv"),
              read_csv("rep_medical_group2.csv")) %>%
  rbindlist(idcol = "group")

### ANALYSIS: CROSS-SECTIONAL --------------------------------------------------
## T1: How many significant hypotheses in stage 2 also occur in stage 1, with the same effect direction?
ls_cs_t1 <- ls_cs %>%
  filter(pvalue < 0.05) %>%
  mutate(direction = ifelse(odds <= 1, -1, 1)) %>%
  select(group, term, direction)

# number of sig hypotheses in stage 1
n_cs_t1_1 <- nrow(ls_cs_t1[ls_cs_t1$group == "1", ])

# number of sig hypotheses in stage 2 that have the same direction
n_cs_t1_2 <- intersect(ls_cs_t1[ls_cs_t1$group == "1", -1], ls_cs_t1[ls_cs_t1$group == "2", -1]) %>%
  nrow()

## T2: Median original and replication effect sizes
# convert effect size to r (which is 0,1 bounded and nicely interpretable)
# cf. https://easystats.github.io/effectsize/reference/d_to_r.html
ls_cs_t2 <- ls_cs %>%
  mutate(effect = effectsize::oddsratio_to_r(odds, Y, N)) %>%
  na.omit

# effect size estimates between stages
p_cs_t2 <- ls_cs_t2 %>%
  group_by(group) %>%
  summarize(mean = mean(effect), sd = sd(effect), median = median(effect))

# p-value: t-test to compare effect size estimates between stages
stats_cs_t2 <- ls_cs_t2 %>%
  summarize(pval = t.test(effect ~ group)$p.value)

## T3: Percent original effect size within replication 95% CI
ls_cs_t3 <- ls_cs %>%
  filter(group == "1") %>%
  mutate(is_within = FALSE)

for(i in 1:nrow(ls_cs_t3)){
  row_stage1 <- ls_cs[i ,]
  row_stage2 <- ls_cs[ls_cs$term == row_stage1$term & ls_cs$group == "2", ]
  
  # get OR 95% CI; fails for very small sample sizes; if so, skip
  possibleError <- tryCatch(
    expr = {
      row_stage1 <- getCI(row_stage1)
      row_stage2 <- getCI(row_stage2)
    },
    error = function(e){ 
      e
    }
  )
  
  if(inherits(possibleError, "error")) next
  
  # check if OR of group 1 is within CI1 and CI2 of group 2
  is_within <- (row_stage1$odds > row_stage2$CI1) & (row_stage1$odds < row_stage2$CI2)
  
  # store in dataframe of just group 1 hypotheses
  ls_cs_t3[i ,]$is_within <- is_within
}

stats_cs_t3 <- nrow(ls_cs_t3[ls_cs_t3$is_within == "TRUE", ])/nrow(ls_cs_t3[ls_cs_t3$is_within == "FALSE", ])

## PLOT
ls_cs_dfp <- ls_cs_t2 %>%
  select(term, group, effect) %>%
  pivot_wider(names_from = "group", values_from = "effect") %>%
  # keep p-value
  left_join(ls_cs_t3[ls_cs_t3$group == "1", ][ ,c("term", "pvalue")]) %>%
  mutate(is_sig = ifelse(pvalue < 0.05, 1, 0)) %>%
  mutate(is_sig = as.factor(is_sig))

ls_cs_p <- ls_cs_dfp %>%
  ggplot(aes(x = `1`, y = `2`, color = is_sig)) +
  geom_point() +
  scale_x_continuous(limits = c(-.4, .4)) +
  scale_y_continuous(limits = c(-.4, .4)) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set1") +
  coord_fixed() +
  xlab("Stage 1 Effect Size\nCross-sectional phenotype") +
  ylab("Stage 2 Effect Size\nCross-sectional phenotype") +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(method = "pearson", 
           geom = "label",
           label.npc.y = 0.01, 
           label.x.npc = 0.01)

ls_cs_p <- ggMarginal(ls_cs_p, 
                      type = "histogram",
                      groupFill = TRUE,
                      groupColour = TRUE)

### ANALYSIS: LONGITUDINAL --------------------------------------------------
## T1: How many significant hypotheses in stage 2 also occur in stage 1, with the same effect direction?
ls_lg_t1 <- ls_lg %>%
  filter(pvalue < 0.05) %>%
  mutate(direction = ifelse(odds <= 1, -1, 1)) %>%
  select(group, bin, term, direction)

# number of sig hypotheses in stage 1
n_lg_t1_1 <- nrow(ls_lg_t1[ls_lg_t1$group == "1", ])

# number of sig hypotheses in stage 2 that have the same direction
n_lg_t1_2 <- intersect(ls_lg_t1[ls_lg_t1$group == "1", -1], ls_lg_t1[ls_lg_t1$group == "2", -1]) %>%
  nrow()

## T2: Median original and replication effect sizes
# convert effect size to r (which is 0,1 bounded and nicely interpretable)
# cf. https://easystats.github.io/effectsize/reference/d_to_r.html
ls_lg_t2 <- ls_lg %>%
  mutate(effect = effectsize::oddsratio_to_r(odds, Y, N)) %>%
  na.omit

# effect size estimates between stages
p_lg_t2 <- ls_lg_t2 %>%
  group_by(group) %>%
  summarize(mean = mean(effect), sd = sd(effect), median = median(effect))

# p-value: t-test to compare effect size estimates between stages
stats_lg_t2 <- ls_lg_t2 %>%
  summarize(pval = t.test(effect ~ group)$p.value)

## T3: Percent original effect size within replication 95% CI
ls_lg_t3 <- ls_lg %>%
  filter(group == "1") %>%
  mutate(is_within = FALSE)

for(i in 1:nrow(ls_lg_t3)){
  row_stage1 <- ls_lg[i ,]
  row_stage2 <- ls_lg[ls_lg$term == row_stage1$term & 
                        ls_lg$group == "2" &
                        ls_lg$bin == row_stage1$bin, ]
  
  # get OR 95% CI; fails for very small sample sizes; if so, skip
  possibleError <- tryCatch(
    expr = {
      row_stage1 <- getCI(row_stage1)
      row_stage2 <- getCI(row_stage2)
    },
    error = function(e){ 
      e
    }
  )
  
  if(inherits(possibleError, "error")) next
  
  # check if OR of group 1 is within CI1 and CI2 of group 2
  is_within <- (row_stage1$odds > row_stage2$CI1) & (row_stage1$odds < row_stage2$CI2)
  
  # store in dataframe of just group 1 hypotheses
  ls_lg_t3[i ,]$is_within <- is_within
}

stats_lg_t3 <- nrow(ls_lg_t3[ls_lg_t3$is_within == "TRUE", ])/nrow(ls_lg_t3[ls_lg_t3$is_within == "FALSE", ])

## PLOT
ls_lg_dfp <- ls_lg_t2 %>%
  select(term, bin, group, effect) %>%
  pivot_wider(names_from = "group", values_from = "effect") %>%
  # keep p-value
  left_join(ls_lg_t3[ls_lg_t3$group == "1", ][ ,c("term", "bin", "pvalue")]) %>%
  mutate(is_sig = ifelse(pvalue < 0.05, 1, 0)) %>%
  mutate(is_sig = as.factor(is_sig))

ls_lg_p <- ls_lg_dfp %>%
  ggplot(aes(x = `1`, y = `2`, color = is_sig)) +
  geom_point() +
  scale_x_continuous(limits = c(-.4, .4)) +
  scale_y_continuous(limits = c(-.4, .4)) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set1") +
  coord_fixed() +
  xlab("Stage 1 Effect Size\nLongitudinal phenotype") +
  ylab("Stage 2 Effect Size\nLongitudinal phenotype") +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(method = "pearson", 
           geom = "label",
           label.npc.y = 0.01, 
           label.x.npc = 0.01)

ls_lg_p <- ggMarginal(ls_lg_p, 
                      type = "histogram",
                      groupFill = TRUE,
                      groupColour = TRUE)

### ANALYSIS: MEDICAL ----------------------------------------------------------
# preprocessing: recover odds ratio, which was previously removed for plot
tmp <- getCI(ls_md)[ ,15:17]
ls_md[ ,9:11] <- tmp

## T1: How many significant hypotheses in stage 2 also occur in stage 1, with the same effect direction?
ls_md_t1 <- ls_md %>%
  filter(P < 0.05) %>%
  mutate(direction = ifelse(OR <= 1, -1, 1)) %>%
  select(group, YearsPrescription, MED_NAME, direction)

# number of sig hypotheses in stage 1
n_md_t1_1 <- nrow(ls_md_t1[ls_md_t1$group == "1", ])

# number of sig hypotheses in stage 2 that have the same direction
n_md_t1_2 <- intersect(ls_md_t1[ls_md_t1$group == "1", -1], ls_md_t1[ls_md_t1$group == "2", -1]) %>%
  nrow()

## T2: Median original and replication effect sizes
# convert effect size to r (which is 0,1 bounded and nicely interpretable)
# cf. https://easystats.github.io/effectsize/reference/d_to_r.html
ls_md_t2 <- ls_md %>%
  mutate(effect = effectsize::oddsratio_to_r(OR, Y, N)) %>%
  na.omit

# effect size estimates between stages
p_md_t2 <- ls_md_t2 %>%
  group_by(group) %>%
  summarize(mean = mean(effect), sd = sd(effect), median = median(effect))

# p-value: t-test to compare effect size estimates between stages
stats_md_t2 <- ls_md_t2 %>%
  summarize(pval = t.test(effect ~ group)$p.value)

## T3: Percent original effect size within replication 95% CI
ls_md_t3 <- ls_md %>%
  filter(group == "1") %>%
  mutate(is_within = FALSE)

for(i in 1:nrow(ls_md_t3)){
  row_stage1 <- ls_md[i ,]
  row_stage2 <- ls_md[ls_md$MED_NAME == row_stage1$MED_NAME & 
                        ls_md$group == "2" &
                        ls_md$YearsPrescription == row_stage1$YearsPrescription, ]
  
  # get OR 95% CI; fails for very small sample sizes; if so, skip
  possibleError <- tryCatch(
    expr = {
      row_stage1 <- getCI(row_stage1)
      row_stage2 <- getCI(row_stage2)
    },
    error = function(e){ 
      e
    }
  )
  
  if(inherits(possibleError, "error")) next
  
  # check if OR of group 1 is within CI1 and CI2 of group 2
  is_within <- (row_stage1$OR > row_stage2$CI1) & (row_stage1$OR < row_stage2$CI2)
  
  # store in dataframe of just group 1 hypotheses
  ls_md_t3[i ,]$is_within <- is_within
}

stats_md_t3 <- nrow(ls_md_t3[ls_md_t3$is_within == "TRUE", ])/nrow(ls_md_t3[ls_md_t3$is_within == "FALSE", ])

## PLOT
ls_md_dfp <- ls_md_t2 %>%
  select(MED_NAME, YearsPrescription, group, effect) %>%
  pivot_wider(names_from = "group", values_from = "effect") %>%
  # keep p-value
  left_join(ls_md_t3[ls_md_t3$group == "1", ][ ,c("MED_NAME", "YearsPrescription", "P")]) %>%
  mutate(is_sig = ifelse(P < 0.05, 1, 0)) %>%
  mutate(is_sig = as.factor(is_sig))

ls_md_p <- ls_md_dfp %>%
  ggplot(aes(x = `1`, y = `2`, color = is_sig)) +
  geom_point() +
  scale_x_continuous(limits = c(-.4, .4)) +
  scale_y_continuous(limits = c(-.4, .4)) +
  theme(legend.position = c(.8, .2)) +
  scale_color_brewer(palette = "Set1", 
                     
                     labels = c("P > 0.05", "P < 0.05"), "") +
  coord_fixed() +
  xlab("Stage 1 Effect Size\nMedical treatment") +
  ylab("Stage 2 Effect Size\nMedical treatment") +
  geom_smooth(method = lm, se = FALSE) +
  stat_cor(method = "pearson", 
           geom = "label",
           label.npc.y = 0.01, 
           label.x.npc = 0.01) +
  guides(color = guide_legend(override.aes = list(linetype = c(1, 1),
                                                  alpha = 0,
                                                  shape = c(10, 10))))

ls_md_p <- ggMarginal(ls_md_p, 
                      type = "histogram",
                      groupFill = TRUE,
                      groupColour = TRUE)

### EXPORT ---------------------------------------------------------------------
## effect size correlation plots
ls_p <- cowplot::plot_grid(ls_cs_p, ls_lg_p, ls_md_p, nrow = 1, labels = "AUTO")

pdf(file = "FigRep.pdf",
    width = 12,
    height = 4)

ls_p

dev.off()
