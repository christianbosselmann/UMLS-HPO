
library(rstatix)
library(coin)
library(ggpubr)
library(finalfit)
library(librarian)
librarian::shelf(tidyverse,
                 data.table,
                 ontologyIndex,
                 ontologySimilarity,
                 plotly,
                 ggrepel,
                 DBI,
                 odbc,
                 readxl,
                 survminer)


# SQL connect -------------------------------------------------------------

db <- DBI::dbConnect(odbc::odbc(), "GeneDatabaseIRB22147")

new.df <- DBI::dbGetQuery(db,"SELECT * FROM [IRB22147_2].[PeopleCPTBeforeAge6_2]")

gene.query <- "SELECT [PAT_MRN_ID],[GENE]
  FROM [CCBSDT02].[MARK].[MARK_MRN_POS]
  WHERE [PAT_MRN_ID] = 'CDLK5'"

CDLK5.df <- DBI::dbGetQuery(db,gene.query)

DBI::dbDisconnect(db)

# Clean Data  -----------------------------------------------

new.df <- new.df %>%
	dplyr::mutate(GENEPOS = as.factor(GENEPOS))%>% #GENEPOS 0 not likely, 1 reviewed and confirmed, 2 to be reviewed
      dplyr::mutate(GENEPOS_comb = recode(GENEPOS,
                                '0' = 'N',
                                '1' = 'Y',
                                '2' = 'Y'))%>%
	dplyr::mutate(ContactAge = as.numeric(round(difftime(
			ContactDate,DateOfBirth, units = "weeks") / 
			52,digits = 1)))%>%
	dplyr::mutate(ProcAge = as.numeric(round(difftime(
			ProcedureDateMin,DateOfBirth, units = "weeks") / 
			52,digits = 1)))


person.df <- new.df %>%
      group_by (PatientId,DateOfBirth,Gender,Ethnicity,GENEPOS,GENEPOS_comb,ProcAge) %>%
	summarize(max_age = max(ContactAge, na.rm = TRUE),
		    min_age = min(ContactAge, na.rm = TRUE),
		    median_age = median(ContactAge, na.rm = TRUE)) # %>%
	# dplyr::filter (!(GENEPOS == '0' & max_age < 10)) %>% # filter out people in non gene epilepsy group with no encounter over age 10 
	# dplyr::filter (!(max_age < 10)) # filter out everyone without an encounter over age 10

# Calculate summary statistics -------------------------------------------------------------

collist = c("Gender","Ethnicity","ProcAge","min_age","median_age","max_age")

person.sum <- person.df %>%  
  summary_factorlist("GENEPOS_comb", collist, 
                     p=TRUE, na_include=TRUE)

table(person.df$GENEPOS_comb)

# Encounters over time --------------------------------------------------

	
filtered.df <- new.df %>%
   filter (PatientId %in% person.df$PatientId)

ggpubr::ggdensity(filtered.df, x = "ContactAge",
            add = "mean", rug = FALSE,
            fill = "GENEPOS_comb", palette = c("#00AFBB", "#E7B800"))


# Input HPO data -------------------------------------------------------------

## UMLS to HPO direct mapping file

db2 <- DBI::dbConnect(odbc::odbc(), "MedLanguage")

query.hpo <- "SELECT * FROM [CCBSDT01].[HPO].[HPOtoUMLS]"

hpo.input <- DBI::dbGetQuery(db2,query.hpo)

DBI::dbDisconnect(db2)

# HPO ontology file

fileloc <- "C:/Users/TAYLORS54/Documents/Scripts/GeneticEpilepsy/hp.obo.txt"
ont_hpo <- ontologyIndex::get_ontology(fileloc, 
                                       propagate_relationships = "is_a", 
                                       extract_tags = "minimal") # "everything" gets all properties incl. cross-references to other ontologies

# Clean data --------------------------------------------------------------

## need to reformat hpo.df so if multiple ConceptIds (sep by commas) are turn into multiple rows

hpo.df <- hpo.input %>%
    dplyr::mutate(ConceptId = strsplit(as.character(UMLSCodes), ",")) %>% 
    unnest(ConceptId)

## get list of patients by group, excluding concepts

groups.df <- new.df %>%
  dplyr::select(PatientId,GENEPOS_comb)%>%
  dplyr::distinct()

# values for demographic table sample size and Fisher tests

Y_tot <- groups.df %>% filter(GENEPOS_comb == "Y") %>% nrow() #get total count of patients in Y gene group
N_tot <- groups.df %>% filter(GENEPOS_comb == "N") %>% nrow() #get total count of patients in N gene group

## calculate age at encounter, filter encounters here if necessary
## end with one row per patient per concept
## keeping a count of that concept, min and max age they had an encounter with that concept

pat_concept.df <- new.df %>%
  dplyr::mutate(ConceptId = gsub(" ", "", ConceptID)) %>%
  dplyr::group_by(PatientId,GENEPOS_comb,ConceptId)%>%
  summarise (date_earlist = min(ContactAge),
             date_latest = max(ContactAge),
             n = n())


# Map UMLS to HPO, by group --------------------------------------------------------------

# get umls counts by GENE group (number of patients in that group with at least one)

group_umls.df <- pat_concept.df %>%
  dplyr::select(PatientId,GENEPOS_comb,ConceptId) %>% 
  dplyr::distinct() %>%
  dplyr::group_by(GENEPOS_comb, ConceptId)%>%
  dplyr::summarise(group_count = n())%>%
  pivot_wider (values_from = group_count, names_from = GENEPOS_comb)

# directly map UMLS to HPO, still by group, using SQL mapping table

group_hpo.df <- group_umls.df %>%
  dplyr::filter (ConceptId != "None") %>%
  distinct(ConceptId) %>%
  merge (., hpo.df, by = "ConceptId", all.x = TRUE)

## test for mapping completeness

n_distinct(group_hpo.df$ConceptId)
n_distinct(group_hpo.df$HPOCode)
sum(is.na(group_hpo.df$HPOCode)) ## number of UMLS did not map


# Fisher Test function --------------------------------------------------------------

fish_test_it <- function(g1,g1_out,g2,g2_out,label){
  
  pvalue <- c()
  odds <- c()
  
  for(i in 1:length(g1)){
    
    fish_out <- matrix(c(g1[i],g1_out[i],g2[i],g2_out[i]),ncol =2) %>% fisher.test()
    
    print(fish_out)
    
    pvalue <- c(pvalue,fish_out$p.value)
    odds <- c(odds,fish_out$estimate)
    
  }  
  
  if(label == "pvalue"){
    return(pvalue)
  }else{
    return(odds)
  }
  
  
}

# Map UMLS to HPO at patient level -------------------------------------------------------------------------


## start from pat_concept.df - check filter to see if any filtering logic in place

hpo_concept.df <- pat_concept.df %>%
  dplyr::filter (ConceptId != "None") %>%
  merge (., hpo.df, by = "ConceptId") 

write.csv(hpo_concept.df,"hpo_list_patient.csv", row.names = FALSE)

group_hpo_concept.df <- hpo_concept.df %>%
  dplyr::select(PatientId,GENEPOS_comb, HPOCode, HPOCodeDescription) %>% 
  dplyr::distinct() %>%
  dplyr::group_by(GENEPOS_comb, HPOCode, HPOCodeDescription)%>%
  dplyr::summarise(group_count = n())%>%
  pivot_wider (values_from = group_count, names_from = GENEPOS_comb)%>% 
  replace(is.na(.), 0) 


write.csv(group_hpo_concept.df,"hpo_list_group_1.25.csv", row.names = FALSE)

# Conduct Fisher Test --------------------------------------------------------------

n_tests = nrow(group_hpo_concept.df)

concept_vis_input.df2 <- group_hpo_concept.df %>% 
  mutate(Y_out =Y_tot-Y,
         N_out = N_tot-N) %>% 
  mutate(pvalue = fish_test_it(Y,Y_out,N,N_out,"pvalue"),
         odds = fish_test_it(Y,Y_out,N,N_out,"odds"),
         freq1 = Y/Y_tot,
         freq2 = N/N_tot,
         color_sig = ifelse(pvalue < (0.05/n_tests),"< Bonferroni corr","> Bonferroni corr"),
         size_sel = -log10(pvalue)*4) #add size of dots in the plot 

max_freq <- c(concept_vis_input.df2$freq1,concept_vis_input.df2$freq2) %>% max() 

print(concept_vis_input.df2[concept_vis_input.df2$pvalue < (0.05/n_tests),],n = 50)

# Graph Fisher Test Results --------------------------------------------------------------

concept_vis_input.df2 %>% 
  mutate(expcat_text = ifelse(pvalue < (0.05/n_tests),HPOCodeDescription,NA)) %>% 
  ggplot(aes(x = freq2, y= freq1, color = color_sig))+
    geom_point(aes(size = size_sel),show.legend = FALSE)+
    theme_classic(base_size = 20)+
    coord_cartesian(xlim = c(0,max_freq), ylim = c(0,max_freq))+
    geom_abline(slope = 1, linetype = "dashed")+
    scale_color_manual(values = c("red","black"))+
    labs(y = "Likely Genetic Epilepsy",
         x = "Control Epilepsy Group")+
    geom_label_repel(aes(label = expcat_text), color = "black", max.overlaps = 14)+
    theme(axis.text = element_text(color = "black"),
          axis.line = element_line(color = "black"))+
    guides(color = "none")


# Propagate HPO --------------------------------------------------------------

# fully propagate terms
propagateList <- function(set){
  lapply(set, ontologyIndex::propagate_relations, 
         ontology = ont_hpo, relations = "parents")
}

## map UMLS to HPO by concept by patient

hpo_all.df <- pat_concept.df %>%
  dplyr::filter (ConceptId != "None") %>%
  merge (., hpo.df, by = "ConceptId") %>% 
  group_by(ConceptId) %>% 
  mutate(freq = n()) %>% 
  ungroup() %>% 
  select(-freq)

## propagate each HPO term

group_hpo_anc.df <- hpo_all.df %>%
  dplyr::select(PatientId,GENEPOS_comb, HPOCode,HPOCodeDescription,date_earlist) %>% 
  mutate(anc = propagateList(HPOCode)) %>%
  unnest(anc)%>%
  dplyr::group_by(GENEPOS_comb, anc)%>%
  dplyr::distinct(PatientId) %>% 
  dplyr::summarise(group_count = n())%>%
  pivot_wider (values_from = group_count, names_from = GENEPOS_comb)%>% 
  replace(is.na(.), 0)%>%
  mutate(anc_terms = ont_hpo$name[anc])

write.csv(group_hpo_anc.df,"hpo_anc_group_1.25.csv", row.names = FALSE)


# Conduct Fisher Test --------------------------------------------------------------

n_tests = nrow(group_hpo_anc.df)
n_tests

concept_vis_input.df3 <- group_hpo_anc.df %>% 
  mutate(Y_out =Y_tot-Y,
         N_out = N_tot-N) %>% 
  mutate(pvalue = fish_test_it(Y,Y_out,N,N_out,"pvalue"),
         odds = fish_test_it(Y,Y_out,N,N_out,"odds"),
         freq1 = Y/Y_tot,
         freq2 = N/N_tot,
         color_sig = ifelse(pvalue < (0.05/n_tests),"< Bonferroni corr","> Bonferroni corr"),
         size_sel = -log10(pvalue)*4) #add size of dots in the plot 

max_freq <- c(concept_vis_input.df3$freq1,concept_vis_input.df3$freq2) %>% max() 

print(concept_vis_input.df3[concept_vis_input.df3$pvalue < (0.05/n_tests),],n = 50)

# Graph Fisher Test Results --------------------------------------------------------------

top20sig <- head(sort(concept_vis_input.df3$pvalue,decreasing=FALSE),n=20)

concept_vis_input.df3 %>% 
  mutate(expcat_text = ifelse(pvalue %in% top20sig,anc_terms,NA)) %>% # label only those with 20 lowest p values
  ggplot(aes(x = freq2, y= freq1, color = color_sig))+
    geom_point(aes(size = size_sel),show.legend = FALSE)+
    theme_classic(base_size = 20)+
    coord_cartesian(xlim = c(0,max_freq), ylim = c(0,max_freq))+
    geom_abline(slope = 1, linetype = "dashed")+
    scale_color_manual(values = c("red","black"))+
    labs(y = "Likely Genetic Epilepsy",
         x = "Control Epilepsy Group")+
    geom_label_repel(aes(label = expcat_text), color = "black", max.overlaps = 14)+
    theme(axis.text = element_text(color = "black"),
          axis.line = element_line(color = "black"))+
    guides(color = "none")

# KM graph - seizure --------------------------------------------------------------

## format data for KM graph

HPOlist_seiz <- get_descendants(ont_hpo, "HP:0001250", exclude_roots = FALSE) #HPO code for seizure and code for all descendents

hpo_seizure.df <- pat_concept.df %>%
  dplyr::filter (ConceptId != "None") %>%
  merge (., hpo.df, by = "ConceptId") %>% 
  mutate(hpo_seizure = HPOCode %in% HPOlist_seiz) # add a column to code in logic whether HPOCode is related to seizure or not

true.df <- hpo_seizure.df %>%
  group_by(PatientId,GENEPOS_comb) %>%
  filter(hpo_seizure == TRUE)%>%
  summarize (time = min(date_earlist))%>%
  mutate(phen_occ = 0)

## graph age of earliest seizure encounter for each group - histogram

ggpubr::ggdensity(true.df, x = "time",
            add = "mean", rug = FALSE,
            fill = "GENEPOS_comb", palette = c("#00AFBB", "#E7B800"),
            xlab = "age of first encounter with seizure tag")

false.df <- hpo_seizure.df %>%
  group_by(PatientId,GENEPOS_comb) %>%
  filter(!PatientId %in% true.df$PatientId)%>%
  summarize (time = max(date_latest))%>%
  mutate(phen_occ = 1)
  
km_input.df <- rbind(true.df,false.df)%>%
  select(GENEPOS_comb,time,phen_occ)

fit <- survfit(Surv(time, phen_occ) ~ GENEPOS_comb, data = km_input.df)

survminer::ggsurvplot(
   fit,                     # survfit object with calculated statistics.
   data = km_input.df,               # data used to fit survival curves.
   risk.table = TRUE,       # show risk table.
   pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for 
                            # point estimates of survival curves.
 #  xlim = c(0,500),         # present narrower X axis, but not affect
                            # survival estimates.
   xlab = "Age in Years",   # customize X axis label.
   ylab = "Proability no tag mapping to seizure in HPO",
 #  break.time.by = 100,     # break X axis in time intervals by 500.
   ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE) # show bars instead of names in text annotations
                            # in legend of risk table

# KM graph - Neurodevelopmental disorder--------------------------------------------------------------

## format data for KM graph

HPOlist_phen <- get_descendants(ont_hpo, "HP:0012758", exclude_roots = FALSE) #HPO code for neurodevelopmental and code for all descendents

hpo_phen.df <- pat_concept.df %>%
  dplyr::filter (ConceptId != "None") %>%
  merge (., hpo.df, by = "ConceptId") %>% 
  mutate(hpo_phen = HPOCode %in% HPOlist_phen) # add a column to code in logic whether HPOCode is related to seizure or not

true.df <- hpo_phen.df %>%
  group_by(PatientId,GENEPOS_comb) %>%
  filter(hpo_phen == TRUE)%>%
  summarize (time = min(date_earlist))%>%
  mutate(phen_occ = 0)

## graph age of earliest seizure encounter for each group - histogram

ggpubr::ggdensity(true.df, x = "time",
            add = "mean", rug = FALSE,
            fill = "GENEPOS_comb", palette = c("#00AFBB", "#E7B800"),
            xlab = "age of first encounter with NDD tag")

false.df <- hpo_phen.df %>%
  group_by(PatientId,GENEPOS_comb) %>%
  filter(!PatientId %in% true.df$PatientId)%>%
  summarize (time = max(date_latest))%>%
  mutate(phen_occ = 1)
  
km_input.df <- rbind(true.df,false.df)%>%
  select(GENEPOS_comb,time,phen_occ)

fit <- survfit(Surv(time, phen_occ) ~ GENEPOS_comb, data = km_input.df)

survminer::ggsurvplot(
   fit,                     # survfit object with calculated statistics.
   data = km_input.df,               # data used to fit survival curves.
   risk.table = TRUE,       # show risk table.
   pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for 
                            # point estimates of survival curves.
 #  xlim = c(0,500),         # present narrower X axis, but not affect
                            # survival estimates.
   xlab = "Age in Years",   # customize X axis label.
   ylab = "Proability no tag mapping to NDD in HPO",
 #  break.time.by = 100,     # break X axis in time intervals by 500.
   ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE) # show bars instead of names in text annotations
                            # in legend of risk table

# KM graph - Behavioral abnormality--------------------------------------------------------------

## format data for KM graph

HPOlist_phen <- get_descendants(ont_hpo, "HP:0000708", exclude_roots = FALSE) #HPO code for behavioral abnormality and code for all descendents

hpo_phen.df <- pat_concept.df %>%
  dplyr::filter (ConceptId != "None") %>%
  merge (., hpo.df, by = "ConceptId") %>% 
  mutate(hpo_phen = HPOCode %in% HPOlist_phen) # add a column to code in logic whether HPOCode is related to seizure or not

true.df <- hpo_phen.df %>%
  group_by(PatientId,GENEPOS_comb) %>%
  filter(hpo_phen == TRUE)%>%
  summarize (time = min(date_earlist))%>%
  mutate(phen_occ = 0)

## graph age of earliest seizure encounter for each group - histogram

ggpubr::ggdensity(true.df, x = "time",
            add = "mean", rug = FALSE,
            fill = "GENEPOS_comb", palette = c("#00AFBB", "#E7B800"),
            xlab = "age of first encounter with behavioral abnormality tag")

false.df <- hpo_phen.df %>%
  group_by(PatientId,GENEPOS_comb) %>%
  filter(!PatientId %in% true.df$PatientId)%>%
  summarize (time = max(date_latest))%>%
  mutate(phen_occ = 1)
  
km_input.df <- rbind(true.df,false.df)%>%
  select(GENEPOS_comb,time,phen_occ)

fit <- survfit(Surv(time, phen_occ) ~ GENEPOS_comb, data = km_input.df)

survminer::ggsurvplot(
   fit,                     # survfit object with calculated statistics.
   data = km_input.df,               # data used to fit survival curves.
   risk.table = TRUE,       # show risk table.
   pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for 
                            # point estimates of survival curves.
 #  xlim = c(0,500),         # present narrower X axis, but not affect
                            # survival estimates.
   xlab = "Age in Years",   # customize X axis label.
   ylab = "Proability no tag mapping to NDD in HPO",
 #  break.time.by = 100,     # break X axis in time intervals by 500.
   ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE) # show bars instead of names in text annotations
                            # in legend of risk table

