### UMLS-HPO ANALYSIS OF INDIVIDUALS WITH GENETIC EPILEPSY SYNDROMES -----------
## Helper functions
##
## Author: Christian Bosselmann, MD
##
## Date Created: 2023-02-07
##
## Copyright (c) Christian Bosselmann, 2023
## Email: bosselc@ccf.org
##
### ----------------------------------------------------------------------------

### visualize relative term frequency as enrichment plot
#' @param data data.frame of Surgery, ConceptID (UMLS), term (HPO), and two groups to compare (TRUE; FALSE)
#' @param ontology HPO ontologyIndex object
#' @param forest logical flag; whether to calculate OR and draw a Forest plot
#' @param qq logical flag; whether to plot a qq plot with CGEN
#' @return res list of data.frame and plot object
enrichmentPlot <- function(data,
                           ontology,
                           forest = FALSE,
                           qq = TRUE){
  
  # object to return
  res <- list()
  
  # read term names from ontologyIndex object
  desc_map <- tibble(term = ontology$id,
                     description = ontology$name)
  
  # reshape to dataframe of terms and columns Y/N for term count in group/!group
  df_group <- data %>% 
    group_by(term) %>%
    count(group) %>%
    pivot_wider(names_from = group, values_from = n) %>%
    rename(N = `FALSE`, Y = `TRUE`) %>%
    replace(is.na(.), 0)
  
  # add in term descriptions for graph
  df_group <- left_join(df_group, desc_map, by = "term")
  
  # # filter: by ancestors, as in longitudinal plot
  # min_set <- ontologyIndex::minimal_set(ontology, df_group$term)
  # df_group <- df_group[df_group$term %in% min_set, ]
  
  # data characteristics
  Y_tot <- max(df_group$Y) # total number of observations in group
  N_tot <- max(df_group$N) # total number of observations in !group
  n_tests <- nrow(df_group) # number of multiple tests
  
  # setup plot
  concept_vis_input.df3 <- df_group %>% 
    mutate(Y_out = Y_tot-Y,
           N_out = N_tot-N) %>% 
    mutate(pvalue = fish_test_it(Y, Y_out,N, N_out, "pvalue"),
           odds = fish_test_it(Y, Y_out, N, N_out, "odds"),
           freq1 = Y/Y_tot,
           freq2 = N/N_tot,
           color_sig = ifelse(p.adjust(pvalue, "bonferroni") < 0.001, "<", ">"),
           size_sel = -log10(pvalue)*4) %>%
    ungroup() %>%
    mutate(pvalue = p.adjust(pvalue, "bonferroni"))
  # %>%
  #   filter(freq1 > 0.05 | freq2 > 0.05) # minimum term frequency filter
  
  # keep dataframe
  res$data <- concept_vis_input.df3
  
  max_freq <- c(concept_vis_input.df3$freq1, concept_vis_input.df3$freq2) %>% max() 
  
  # n of top p-values to label
  top_sig <- head(sort(concept_vis_input.df3$pvalue, decreasing = FALSE), n = 15) 
  
  # plot
  res$plot <- concept_vis_input.df3 %>%
    mutate(expcat_text = ifelse(pvalue %in% top_sig, description, NA)) %>%
    ggplot(aes(x = freq2, y = freq1, color = color_sig)) +
    geom_point(aes(size = size_sel), show.legend = FALSE) +
    theme_classic(base_size = 20) +
    # coord_fixed(xlim = c(0, 0.25), ylim = c(0, 0.25)) +
    coord_cartesian(xlim = c(0, max_freq), ylim = c(0, max_freq)) +
    geom_abline(slope = 1, linetype = "dashed") +
    scale_color_manual(values = c("red", "black")) +
    labs(y = "Case",
         x = "Control") +
    geom_label_repel(aes(label = expcat_text), 
                     color = "black", 
                     max.overlaps = 8, 
                     size = 3, 
                     force_pull = 0.2,
                     min.segment.length = 0) +
    theme(axis.text = element_text(color = "black"),
          axis.line = element_line(color = "black")) +
    guides(color = "none")
  
  # optional: return a Forest plot
  if(forest == TRUE){
    # redo OR with confidence interval, assign to preallocated dataframe
    concept_odds <- data.frame(OR = 1:nrow(concept_vis_input.df3),
                               CI1 = 1:nrow(concept_vis_input.df3),
                               CI2 = 1:nrow(concept_vis_input.df3))
    
    for(i in 1:nrow(concept_vis_input.df3)){
      row <- concept_vis_input.df3[i, ] %>%
        ungroup() %>%
        select(Y, Y_out, N, N_out) %>%
        data.matrix()
      mat <- matrix(data = row, nrow=2)
      fish <- fisher.test(mat)
      concept_odds$OR[i] <- fish$estimate
      concept_odds$CI1[i] <- fish$conf.int[[1]] # lower bound
      concept_odds$CI2[i] <- fish$conf.int[[2]] # upper bound
    }
    
    # merge with full df
    df_concept <- cbind(concept_vis_input.df3, concept_odds)
    
    df_concept <- df_concept %>%
      ungroup() 
    # # only keep significant OR
    # filter(CI1 > 1) %>%
    # # only keep significant observations after correction
    # filter(color_sig == "<") 
    
    # # filter: by ancestors
    # min_set <- ontologyIndex::minimal_set(ont_hpo, df_concept$term)
    # df_concept <- df_concept[df_concept$term %in% min_set, ]
    
    # filter: custom terms to display on Forest plot
    vec_terms <- c("Abnormality of metabolism/homeostasis",
                   "Abnormality of the immune system",
                   "Abnormality of the genitourinary system", 
                   "Abnormality of the skeletal system",
                   "Abnormality of the cardiovascular system",
                   "Abnormality of the nervous system",
                   "Abnormality of the digestive system",
                   "Abnormal respiratory system physiology"
    )
    
    # plot
    res$forest <- df_concept %>%
      # # keep n best
      # slice_max(order_by = OR, n = 8, with_ties = FALSE) %>%
      filter(description %in% vec_terms) %>%
      ggplot(aes(y = description)) +
      # ggplot(aes(y = reorder(description, OR))) +
      geom_point(aes(x = OR), shape = 15, size = 3) +
      geom_linerange(aes(xmin = CI1, xmax = CI2)) +
      geom_vline(xintercept = 1, linetype = "dashed") +
      scale_x_continuous(trans = 'log10') +
      scale_y_discrete(labels = label_wrap(30)) +
      expand_limits(x = 1) +
      theme_classic() +
      ylab("") +
      xlab("Odds ratio (95% CI, log scale)")
  }
  
  # diagnostics: QQ plot
  concept_vis_input.df3 %>%
    ungroup() %>%
    pull(pvalue) %>%
    QQ.plot(.)
  
  abline(v = -log10(0.05), col = "blue")
  
  return(res)
}

### visualize kernel matrix
#' @param K kernel matrix to be visualized
#' @param hc.order boolean, see hc.order attribute of ggcorrplot fn
#' @return ggplot2 object
kernelVisualization <- function(K, hc.order = TRUE){
  library(ggplot2)
  library(ggcorrplot)
  
  ggcorrplot(K,
             hc.order = hc.order,
             hc.method = "average",
             outline.color = NA,
             legend.title = "  r",
             colors = c("#6D9EC1", "white", "#E46726"))
}

### negate in
`%nin%` = Negate(`%in%`)

### normalize kernel matrix
# cf. Kernel Methods for Pattern Analysis, Algorithm 5.1
#' @param K kernel matrix to be normalized
#' @return normalized kernel matrix
kernelNormalisation <- function(K){
  # min_eigenvalue <- min(eigen(K)$values)
  # K <- sqrt(diag(K) + min_eigenvalue)
  D <- diag(1/sqrt(diag(K)))
  K <- D %*% K %*% D
  return(K)
}

### center in feature space
# cf. Kernel Methods for Pattern Analysis, Algorithm 5.3
#' @param K kernel matrix to be centered
#' @return kernel matrix centered in feature space
kernelCentering <- function(K){
  ell <- dim(K)[1]
  D <- colSums(K)/ell # row vector storing the column averages of K
  E <- sum(D)/ell # average of all the entries of K
  J <- matrix(1, ell, 1) %*% D
  Jt <- Conj(t(J)) # complex conjugate transpose of J
  K <- K - J - Jt + E * matrix(1, ell, ell)
}

### check if matrix is symmetric and psd
#' @param K kernel matrix to be checked
#' @return prints matrix properties to console
#' adapted from base and matrixcalc
kernelCheck <- function(K, tol = 1e-08){
  if(!isSymmetric(K)){stop("Argument is not a symmetric matrix.")}
  if(isSymmetric(K)){print("Argument is a symmetric matrix.")}
  
  ev <- eigen(K)[[1]]
  n <- nrow(K)
  
  for (i in 1:n) {
    if (abs(ev[i]) < tol) {
      ev[i] <- 0
    }
  }
  if (any(ev < 0)) {
    stop("Argument is not a psd matrix.")
  }
  print("Argument is a psd matrix.")
}

### standard kernel preprocessing
#' @param K kernel matrix to be normalized and centered
#' @return kernel matrix 
kernelPreparation <- function(K){
  K <- kernelNormalisation(K)
  K <- kernelCentering(K)
  K <- round(K, 10)
  return(K)
}

### integrated pairwise phenotypic similarity
#' @params term_list list of character vectors of HPO terms, where each element in the list is a patient or variant
#' @params ontology ontologyIndex object
#' @params method similarity measure, choice of c("jaccard", "lin", "resnik", "euclidean", "cosine")
#' @returns mat_pheno a phenotypic similarity kernel matrix (psd)
pairwiseSimilarity <- function(term_list, ontology, method){
  # pkg
  library(ontologyIndex)
  library(ontologySimilarity)
  library(proxy)
  library(klic)
  library(tidyverse)
  
  # check method
  if (method %nin% c("jaccard", "lin", "resnik", "euclidean", "cosine")) stop("Invalid method.")
  
  if (method %in% c("jaccard", "euclidean", "cosine")) {
    # get pairwise phenotypic similarity: jaccard
    ls_prop <- lapply(term_list, 
                      propagate_relations, 
                      ontology = ontology, 
                      relations = "parents")
    
    df_prop <- ls_prop %>% 
      map_dfr(~ .x %>% as_tibble(), .id = "name")
    
    df_prop <- reshape2::dcast(df_prop, name ~ value, length)
    
    df_prop <- df_prop[,-1] # remove id col
    df_prop[df_prop > 0] <- 1 # enforce binary
    
    mat_pheno <- proxy::simil(x = df_prop,
                              method = method) # "jaccard", "euclidean", "cosine"
    
    mat_pheno <- as.matrix(mat_pheno)
    diag(mat_pheno) <- 1
  }else{
    # get pairwise phenotypic similarity: lin/resnik
    ic <- get_term_info_content(ontology, term_list, patch_missing = FALSE)
    mat_pheno <- get_sim_grid(ontology = ontology, 
                              information_content = ic,
                              term_sim_method = method, # i.e. "lin", "resnik"
                              term_sets = term_list)
    
    mat_pheno <- klic::spectrumShift(mat_pheno, coeff = 1.2) # nearest psd
  }
  
  # sanity check
  kernelCheck(mat_pheno)
  
  # return phenotypic kernel matrix
  return(mat_pheno)
}

#' this function takes a dataframe of variables for bnlearn and
#' removes constant and/or perfectly correlated columns
#' importantly, perfectly correlated columns are "grouped"
#' the colnames of removed columns are preserved
#' @param vars df of binary variables
#' @param remove_constants boolean; if TRUE, remove constant columns
#' @param remove_correlated boolean; if TRUE, remove/grouped corr. columns
#' @param threshold for remove_correlated, the correlation of threshold
#' @param as_factor boolean; if TRUE, coerce variables to factors for BDs score
preprocessVariables <- function(vars, 
                                remove_constants = TRUE,
                                remove_correlated = TRUE,
                                threshold = 0.95,
                                as_factor = TRUE){
  
  # pkg
  library(igraph)
  library(bnlearn)
  
  if(remove_constants){
    vars <- vars[sapply(vars, function(x) length(unique(na.omit(x)))) > 1]
  }
  
  if(remove_correlated){
    df_cor <- cor(as.matrix(vars), method = "pearson")
    var_cor <- df_cor*lower.tri(df_cor)
    df_cor <- which(df_cor >= threshold, arr.ind = TRUE)
    graph_cor <- igraph::graph.data.frame(df_cor, directed = FALSE)
    groups_cor <- split(unique(as.vector(df_cor)), 
                        clusters(graph_cor)$membership)
    groups_cor <- lapply(groups_cor, function(x) {rownames(var_cor)[x]})
    
    for (i in 1:length(groups_cor)){
      # keep the first element, discard the rest but keep their HPO IDs in colname
      str_cols <- groups_cor[[i]]
      to_keep <- str_cols[1]
      to_drop <- str_cols[-1]
      if(is_empty(to_drop)) next
      grp_name <- paste(str_cols, collapse = " ")
      vars <- vars[, -which(names(vars) %in% to_drop)] 
      names(vars)[names(vars) == to_keep] <- grp_name
    }
  }
  
  if(as_factor){
    dim <- ncol(vars)
    vars[,1:dim] <- lapply(vars[,1:dim], function(x) factor(x, level = c(0, 1)))
  }
  
  return(vars)
}

# Fisher's test function from Sara
fish_test_it <- function(g1,g1_out,g2,g2_out,label){
  pvalue <- c()
  odds <- c()
  
  for(i in 1:length(g1)){
    fish_out <- matrix(c(g1[i],g1_out[i],g2[i],g2_out[i]),ncol =2) %>% fisher.test()
    pvalue <- c(pvalue,fish_out$p.value)
    odds <- c(odds,fish_out$estimate)
  }  
  if(label == "pvalue"){
    return(pvalue)
  }else{
    return(odds)
  }
}

# flat violin plot from GIST 2a1bb0133ff568cbe28d (Ben Marwick)
library(ggplot2)
library(dplyr)

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )

#' Plot function for p-value of HPO terms over age, for 01-genetics.R
#' @param df_genes df with cols PatientId, ConceptID, ContactAge, ProcAge, status
#' @param df_match1 df subset after matching procedure with column group for cohort membership
#' @param show_legend ggplot2 legend.position, allowed are: “left”, “top”, “right”, “bottom”.
#' @param fix_x int; fixed x-axis length
#' @param odds_plot logical flag; if true, generate the plot with OR y-axis
#' @param filter logical flag; if true, only keep the term in the bin with the highest p-value
#' @return res list of longitudinal plot
longitudinalPlot <- function(df_genes, df_match1, 
                             show_legend = "none", fix_x = 25,
                             odds_plot = FALSE,
                             filter = TRUE){
  res <- list()
  
  # create age bins
  df_g1 <- df_genes %>%
    # subgroup by matched patients
    filter(PatientId %in% df_match1$PatientId) %>%
    # get group label
    left_join(df_match1[ ,c("PatientId", "group")] %>% unique, by = "PatientId") %>%
    ungroup() %>%
    ## bin width based on encounter frequency
    # mutate(bin = cut_number(ContactAge, n = 10))
    ## fixed bin width
    mutate(bin = cut(ContactAge, breaks = c(0, 2, 12, 18, Inf)))
  
  # get numeric breaks for histogram
  breaks_binned <- levels(df_g1$bin) %>%
    sapply(., function(x) {gsub("\\,", " ", x)}) %>%
    parse_number()
  
  breaks_binned <- c(breaks_binned, max(df_g1$ContactAge))
  
  # split into age bins
  ls_g1 <- df_g1 %>% 
    split(.$bin)
  
  # for each bin, map to HPO terms, then map to propagated HPO terms
  ls_m1 <- list()
  for(i in 1:length(ls_g1)){
    ls_m1[[i]] <- left_join(ls_g1[[i]], hpo_map, by = "ConceptID") %>%
      rename(term = name) %>%
      na.omit()
    
    ls_m1[[i]] <- left_join(ls_m1[[i]], prop_map, by = "term")
  }
  
  # count of each genetic subgroup
  # run Fisher's test
  # get n most significant terms per bin
  ls_p1 <- list()
  for(i in 1:length(ls_m1)){
    # check if bin contains observations for both groups
    # this is not the case, skip this bin
    if(length(table(ls_m1[[i]]$group)) == 1){next}
    
    df_group <- ls_m1[[i]] %>%
      select(PatientId, group, prop_terms) %>%
      unnest(cols = c(prop_terms)) %>%
      group_by(prop_terms, group) %>%
      count(prop_terms) %>%
      pivot_wider(names_from = group, values_from = n) %>%
      replace(is.na(.), 0) %>%
      rename(term = prop_terms,
             Y = `TRUE`,
             N = `FALSE`)
    
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
             color_sig = ifelse(p.adjust(pvalue, "bonferroni") < 0.05, "<", ">"),
             size_sel = -log10(pvalue)*4) %>%
      ungroup() %>%
      mutate(pvalue = p.adjust(pvalue, "bonferroni"))
    
    if(filter){
      # filter by positive ORs: we only want observations for the cases
      df_group <- df_group %>%
        filter(odds > 1)
    }
    
    # # commented out to keep all raw data 
    # # get n best p-values for each bin for longitudinal plot
    # df_group <- df_group %>%
    #   ungroup() %>%
    #   select(term, description, pvalue, odds) 
    # # %>% 
    # #   slice_min(order_by = pvalue, n = 8, with_ties = FALSE) # non-trivial to choose
    
    if(filter){
      # filter: by ancestors
      min_set <- ontologyIndex::minimal_set(ont_hpo, df_group$term)
      df_group <- df_group[df_group$term %in% min_set, ]
    }
    
    # return
    ls_p1[[i]] <- df_group
  }
  
  # fix bin to correspond to mean age of bin for graph
  seq <- seq(1, length(breaks_binned), 1)
  breaks_mean <- sapply(seq, function(i) {mean(breaks_binned[i:(i+1)])}) %>%
    na.omit
  
  names(ls_p1) <- breaks_mean[1:length(ls_p1)]
  
  df_gp1 <- ls_p1 %>%
    rbindlist(idcol = "bin") %>%
    mutate(bin = as.numeric(bin))
  
  # filter: for each term, keep only the bin with the highest p-value
  if(filter){
    
    df_gp1 <- df_gp1 %>%
      group_by(term) %>% 
      slice_min(order_by = pvalue, n = 1)
    
    # filter: only needed if we include many or insignificant variables
    # only keep description labels for n top pvalues per bin
    # these still have to be significant
    top_labels <- df_gp1 %>%
      group_by(bin) %>%
      slice_min(order_by = pvalue, n = 2) %>%
      filter(pvalue < 0.05)
    
    df_gp1 <- df_gp1 %>%
      group_by(bin) %>%
      mutate(description = ifelse(description %in% top_labels$description, description, NA))
    
  } # end filter
  
  # point plot: log10(pvalue) over age bins
  # here we can just plot genetic vs non-genetic as sanity check for the subanalysis
  palette <- colorRampPalette(RColorBrewer::brewer.pal(9, name = 'RdBu'))(length(breaks_mean))
  
  res$plot <- df_gp1 %>%
    ggplot(aes(x = bin, y = -log10(pvalue), fill = factor(bin, levels = breaks_mean))) +
    geom_point() +
    geom_jitter() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_label_repel(aes(label = description), size = 3.0,
                     color = "black", max.overlaps = Inf,
                     force_pull = 0.01,
                     min.segment.length = 0) +
    scale_fill_manual(values = palette, 
                      name = "Mean age (years)",
                      breaks = breaks_mean, 
                      labels = format(round(breaks_mean, 3), nsmall = 1),
                      guide = guide_legend(override.aes = list(label = ""))) +
    scale_x_continuous(expand = expand_scale(mult = c(0.1, 0.1))) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) +
    theme_classic() +
    theme(legend.position = show_legend) +
    coord_cartesian(xlim = c(0, fix_x)) +
    xlab("Age (years)")
  
  # OR plot
  if(odds_plot){
    res$plot_or <- df_gp1 %>%
      ggplot(aes(x = bin, y = odds, fill = factor(bin, levels = breaks_mean))) +
      geom_point() +
      geom_jitter() +
      geom_label_repel(aes(label = description), size = 3.0,
                       color = "black", max.overlaps = Inf,
                       force_pull = 0.01,
                       min.segment.length = 0) +
      scale_fill_manual(values = palette, 
                        name = "Mean age (years)",
                        breaks = breaks_mean, 
                        labels = format(round(breaks_mean, 3), nsmall = 1),
                        guide = guide_legend(override.aes = list(label = ""))) +
      theme_classic() +
      theme(legend.position = show_legend) +
      coord_cartesian(xlim = c(0, fix_x)) +
      xlab("Age (years)") +
      ylab("Odds ratio (censored, log scale)") +
      scale_x_continuous(expand = expand_scale(mult = c(0.1, 0.1))) +
      scale_y_continuous(oob = scales::oob_squish_infinite,
                         trans = scales::log_trans(base = 10),
                         expand = expand_scale(mult = c(0, .1)))
  }
  return(res)
}

#' This function takes a dataframe of {PatientId, status} and a dataframe
#' of {PatientId, ConceptID and ContactAge}, merges them and finds the number of
#' encounters per patient per group; encounters from the majority group are
#' downsampled to correct for encounter frequency imbalance in the dataset
#' @param df_match1 dataframe of matched case-control pairs
#' @param df essentially a lookup table of ContactAge and ConceptID by PatientId
#' @param verbose logical flag; whether to print the group imbalance to the console
#' @returns df_ss dataframe like df_match1, downsampled to equal groups
downsampleMatch <- function(df_match1, df, 
                            verbose = TRUE){
  df_ss <- df_match1 %>%
    left_join(df[ ,c("PatientId", "ContactAge", "ConceptID")], by = "PatientId") 
  
  # force unique ConceptIDs per ContactAge
  df_ss <- df_ss %>%
    distinct(PatientId, ConceptID, ContactAge, .keep_all = TRUE)
  
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
  
  if(verbose){
    print("Ratio of mean number of unique encounters per patient per group:")
    print(ratio_imb)
  }
  
  # sample a fraction of encounters per patient from the majority group
  df_ss_min <- df_ss %>% 
    # ungroup() %>%
    group_by(PatientId) %>%
    filter(status == label_maj) %>%
    slice_sample(prop = ratio_imb)
  
  # rowbind back with the minority group
  df_ss <- df_ss %>%
    ungroup() %>%
    filter(status != label_maj) %>%
    rbind(df_ss_min)
  
  # unnest ConceptIDs again
  df_ss <- df_ss %>% 
    unnest(ConceptID)
  
  return(df_ss)
}

### QQ-Plot and lambda value functions from https://slowkow.com/notes/ggplot2-qqplot/
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

#' Create a quantile-quantile plot with ggplot2.
#'
#' Assumptions:
#'   - Expected P values are uniformly distributed.
#'   - Confidence intervals assume independence between tests.
#'     We expect deviations past the confidence intervals if the tests are
#'     not independent.
#'     For example, in a genome-wide association study, the genotype at any
#'     position is correlated to nearby positions. Tests of nearby genotypes
#'     will result in similar test statistics.
#'
#' @param ps Vector of p-values.
#' @param ci Size of the confidence interval, 95% by default.
#' @return A ggplot2 plot.
#' @examples
#' library(ggplot2)
#' gg_qqplot(runif(1e2)) + theme_grey(base_size = 24)
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}
