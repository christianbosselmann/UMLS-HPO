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
  
  # plot
  concept_vis_input.df3 <- df_group %>% 
    mutate(Y_out = Y_tot-Y,
           N_out = N_tot-N) %>% 
    mutate(pvalue = fish_test_it(Y,Y_out,N,N_out,"pvalue"),
           odds = fish_test_it(Y,Y_out,N,N_out,"odds"),
           freq1 = Y/Y_tot,
           freq2 = N/N_tot,
           color_sig = ifelse(p.adjust(pvalue, "bonferroni") < 1e-5, "<", ">"),
           size_sel = -log10(pvalue)*4)
  
  res$data <- concept_vis_input.df3
  
  max_freq <- c(concept_vis_input.df3$freq1, concept_vis_input.df3$freq2) %>% max() 
  
  top_sig <- head(sort(concept_vis_input.df3$pvalue, decreasing = FALSE), n = 15) # n of top p-values to label
  
  res$plot <- concept_vis_input.df3 %>% 
    mutate(expcat_text = ifelse(pvalue %in% top_sig, description, NA)) %>% # label only those with 20 lowest p values
    ggplot(aes(x = freq2, y = freq1, color = color_sig)) +
    geom_point(aes(size = size_sel), show.legend = FALSE) +
    theme_classic(base_size = 20) +
    # coord_fixed(xlim = c(0, 0.25), ylim = c(0, 0.25)) +
    coord_cartesian(xlim = c(0, max_freq), ylim = c(0, max_freq)) +
    geom_abline(slope = 1, linetype = "dashed") +
    scale_color_manual(values = c("red", "black")) +
    labs(y = "Case",
         x = "Control") +
    geom_label_repel(aes(label = expcat_text), color = "black", max.overlaps = 12, size = 3, force_pull = 0.5) +
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



