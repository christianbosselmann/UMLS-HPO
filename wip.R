# downsampling

# define cohorts and match by age, sex and ethnicity; maintain label
df_match <- df_genes %>%
  left_join(df_person[,c("PatientId", "Gender", "Ethnicity", "median_age")], by = "PatientId") %>%
  distinct(PatientId, Gender, Ethnicity, median_age, status) %>%
  # set covariates as factors
  mutate(Gender = as.factor(Gender)) %>%
  mutate(Ethnicity = as.factor(Ethnicity))

## Group 1: genetic vs. non-genetic
df_match1 <- df_match %>%
  mutate(status = recode(status, 
                         "nongenetic" = 0,
                         "genetic" = 1,
                         "scn1a" = 1,
                         "cdkl5" = 1)) 

df_match1 <- matchit(status ~ median_age + Ethnicity + Gender, 
                     data = df_match1, ratio = 1,
                     method = flag_match, distance = "glm")

df_match1 <- match.data(df_match1)

###HERE
df_ss <- df_match1 %>%
  left_join(df[ ,c("PatientId", "ContactAge", "ConceptID")], by = "PatientId") 

# recode ContactAge as unique values (encounters) per PatientId
df_ss <- df_ss %>%
  group_by(PatientId, ContactAge) %>%
  nest(ConceptID = c(ConceptID))

# find the number of unique encounters per patient
vec_imb <- df_ss %>%
  ungroup() %>%
  group_by(status, PatientId) %>%
  mutate(mean = mean(n_distinct(ContactAge)))

# find the mean number of unique encounters per group
vec_imb <- vec_imb %>%
  group_by(status) %>%
  summarize(mean = mean(mean))

# downsampling: find majority group and ratio
label_maj <- vec_imb$status[which.max(vec_imb$mean)]
ratio_imb <- vec_imb$mean[[which.min(vec_imb$mean)]]/vec_imb$mean[[which.max(vec_imb$mean)]]

# sample a fraction of encounters per patient from the majority group
df_ss_min <- df_ss %>% 
  ungroup() %>%
  filter(status == label_maj) %>%
  slice_sample(prop = ratio_imb)

# rowbind back with the minority group
df_ss <- df_ss %>%
  ungroup() %>%
  filter(status != label_maj) %>%
  rbind(df_ss_min)

###
df_match1 <- df_match1 %>%
  # merge in ConceptIDs
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


