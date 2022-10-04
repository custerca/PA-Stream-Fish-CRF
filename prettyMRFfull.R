prettyMRF = function(data, 
                     MRF_mod, 
                     node_names, 
                     covariate, 
                     obscovquan,
                     cutoff, 
                     plot,
                     noi=NULL,
                     node_color = 'gray',
                     noi_color = 'gold',
                     neg_color = 'red4',
                     pos_color = 'blue',
                     alphanet = 0.5,
                     edge_cutoff = 0.5,
                     vcex=2,
                     vsize = 20){
  
  coef_matrix <- MRF_mod$direct_coef_means
  interaction_coefficients <- coef_matrix[, 2:(nrow(coef_matrix) + 1)]  +
    (Reduce(`+`, MRF_mod$indirect_coef_mean) /
       length(MRF_mod$indirect_coef_mean))
  
  if(missing(node_names)){
    node_names <- rownames(coef_matrix)
  }
  dimnames(interaction_coefficients) <- list(node_names, node_names)
  
  #### Extract indirect effect matrix that matches the covariate name ####
  indirect_coef_names <- names(MRF_mod$indirect_coef_mean)
  which_matrix_keep <- grepl(covariate, indirect_coef_names)
  covariate_matrix <- MRF_mod$indirect_coef_mean[which_matrix_keep][[1]]
  rownames(covariate_matrix) <- node_names
  colnames(covariate_matrix) <- node_names
  baseinteraction_matrix <- interaction_coefficients
  
  #observed_cov_quantiles <- obscovquan
  
  #par(mfrow = c(1, length(observed_cov_quantiles)))
  pred_values <- (covariate_matrix * obscovquan) + baseinteraction_matrix
  cont.cov.mats <- igraph::graph.adjacency(pred_values, weighted = T, mode = "undirected")
  
  mat1 <- cont.cov.mats
  
  
  # Change node color
  if(!is.null(noi)){
    # select the nodes having these names
    selnodes <- V(mat1)[name %in% noi]
    
    # get their network neighborhood 
    selegoV <- ego(mat1, order=1, nodes = selnodes, mode = "all", mindist = 0)
    
    # turn the returned list of igraph.vs objects into a graph
    selegoG <- induced_subgraph(mat1,unlist(selegoV))
    
    igraph::V(selegoG)$color <- node_color
    V(selegoG)$color[V(selegoG)$name==noi] <- noi_color
    
  }else{
    # select the nodes having these names
    selnodes <- V(mat1)
    
    # get their network neighborhood 
    selegoV <- ego(mat1, order=1, nodes = selnodes, mode = "all", mindist = 0)
    
    # turn the returned list of igraph.vs objects into a graph
    selegoG <- induced_subgraph(mat1,unlist(selegoV))
    
    igraph::V(selegoG)$color <- node_color

  }

  # cols <- c(grDevices::adjustcolor(neg_color,alpha.f = alphanet),
  #           grDevices::adjustcolor(pos_color,alpha.f = alphanet))
  
  # Change edge color
  igraph::E(selegoG)$color <- ifelse(igraph::E(selegoG)$weight < 0, neg_color, pos_color)
  
  # Delete edges that represent weak interactions
  selegoG <- igraph::delete.edges(selegoG, which(abs(igraph::E(selegoG)$weight) <= edge_cutoff))
  # plot the subgraph
  V(selegoG)$name <-  str_replace_all(V(selegoG)$name, "[.]", "_")
  V(selegoG)$name <- str_replace_all(V(selegoG)$name, "[_]", " ")
  # Adjust weights for easier visualisation
  # of interaction strengths
  igraph::E(selegoG)$width <- abs(igraph::E(selegoG)$weight) * 6
  # Change vertices labels
  n.v <- length(V(selegoG)$name) # number of vertices
  # # Change labels to letters
  # if(node.letters){
  #   V(selegoG)$name <- LETTERS[1:n.v]
  # }
  # Label font
  V(selegoG)$label.font = 2
  V(selegoG)$label.cex = 0.5
  # Circle layout
  # l1 <- layout_in_circle(selegoG)
  # plot(selegoG,
  #      layout=l1,
  #      vertex.label=V(selegoG)$name,
  #      vertex.label.cex=vcex,
  #      vertex.size=vsize)
  return(selegoG)
  
}
