# helper functions for HPO analysis

### visualize relative term frequency as enrichment plot
#' @param data data.frame of Surgery, ConceptID (UMLS), term (HPO), and two groups to compare (TRUE; FALSE)
#' @param ontology HPO ontologyIndex object
#' @return res list of data.frame and plot object
enrichmentPlot <- function(data,
                           ontology){
  
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
  
  # data characteristics
  Y_tot <- max(df_group$Y) # total number of observations in group
  N_tot <- max(df_group$N) # total number of observations in !group
  n_tests <- nrow(df_group) # number of multiple tests
  
  # setup plot
  concept_vis_input.df3 <- df_group %>% 
    mutate(Y_out = Y_tot-Y,
           N_out = N_tot-N) %>% 
    mutate(pvalue = fish_test_it(Y,Y_out,N,N_out,"pvalue"),
           odds = fish_test_it(Y,Y_out,N,N_out,"odds"),
           freq1 = Y/Y_tot,
           freq2 = N/N_tot,
           color_sig = ifelse(p.adjust(pvalue, "holm") < 0.001, "<", ">"),
           size_sel = -log10(pvalue)*4) %>%
    filter(freq1 > 0.05 | freq2 > 0.05) # minimum term frequency filter
  
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
    geom_label_repel(aes(label = expcat_text), color = "black", max.overlaps = 8, size = 3, force_pull = 0.4) +
    theme(axis.text = element_text(color = "black"),
          axis.line = element_line(color = "black")) +
    guides(color = "none")
  
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
#' @return res list of longitudinal plot
longitudinalPlot <- function(df_genes, df_match1, show_legend = "none", fix_x = 25){
  res <- list()
  
  # create age bins
  df_g1 <- df_genes %>%
    # subgroup by matched patients
    filter(PatientId %in% df_match1$PatientId) %>%
    # get group label
    left_join(df_match1[ ,c("PatientId", "group")] %>% unique, by = "PatientId") %>%
    ungroup() %>%
    mutate(bin = cut_number(ContactAge, n = 10))
  
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
    
    # check if tbin contains any case observations
    # this is not the case for the CDKL5 cohort, due to sample size limitations
    # if true, then skip this bin
    if(is.null(df_group$Y)){next}
    
    # df_group can also be used for enrichment plots
    df_group <- df_group %>% 
      mutate(Y_out = max(df_group$Y)-Y,
             N_out = max(df_group$N)-N) %>% 
      mutate(pvalue = fish_test_it(Y, Y_out, N, N_out, "pvalue"),
             odds = fish_test_it(Y, Y_out, N, N_out, "odds"),
             freq1 = Y/max(df_group$Y),
             freq2 = N/max(df_group$N),
             color_sig = ifelse(p.adjust(pvalue, "holm") < 0.05, "<", ">"),
             size_sel = -log10(pvalue)*4)
    
    # filter by positive ORs: we only want observations for the cases
    df_group <- df_group %>%
      filter(odds > 1)
    
    # get n best p-values for each bin for longitudinal plot
    df_group <- df_group %>%
      ungroup() %>%
      select(term, description, pvalue) %>% 
      slice_min(order_by = pvalue, n = 8, with_ties = FALSE) # non-trivial to choose
    
    # filter: by ancestors
    min_set <- ontologyIndex::minimal_set(ont_hpo, df_group$term)
    df_group <- df_group[df_group$term %in% min_set, ]
    
    # return
    ls_p1[[i]] <- df_group
  }
  
  # fix bin to correspond to mean age of bin for graph
  seq <- seq(1, length(breaks_binned), 1)
  breaks_mean <- sapply(seq, function(i) {mean(breaks_binned[i:(i+1)])}) %>%
    na.omit
  
  # filter: for each term, keep only the bin with the highest p-value
  names(ls_p1) <- breaks_mean[1:length(ls_p1)]
  
  df_gp1 <- ls_p1 %>%
    rbindlist(idcol = "bin") %>%
    mutate(bin = as.numeric(bin))
  
  df_gp1 <- df_gp1 %>%
    group_by(term) %>% 
    slice_min(order_by = pvalue, n = 1)
  
  # point plot: log10(pvalue) over age bins
  # here we can just plot genetic vs non-genetic as sanity check for the subanalysis
  palette <- colorRampPalette(RColorBrewer::brewer.pal(9, name = 'RdBu'))(length(breaks_mean))
  
  res$plot <- df_gp1 %>%
    ggplot(aes(x = bin, y = -log10(pvalue), fill = factor(bin, levels = breaks_mean))) +
    geom_point() +
    geom_label_repel(aes(label = description), size = 3.0,
                     color = "black", max.overlaps = Inf,
                     force_pull = 0.01) +
    scale_fill_manual(values = palette, 
                      name = "Mean age (years)",
                      breaks = breaks_mean, 
                      labels = format(round(breaks_mean, 3), nsmall = 1),
                      guide = guide_legend(override.aes = list(label = ""))) +
    scale_x_continuous(expand = expand_scale(mult = c(0.1, 0))) +
    scale_y_continuous(expand = expand_scale(mult = c(0.4, 0.4))) +
    theme_classic() +
    theme(legend.position = show_legend) +
    coord_cartesian(xlim = c(0, fix_x)) +
    xlab("Age (years)")
}

