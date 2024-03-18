#' Get Nodes of Networks
#' @description get all potential nodes of each network.
#' @param mat frequency matrix
#' @param K number of CMs
#' @return data frame of all potential nodes of each network.
#' @export
getNode <- function(mat, K){
  # norm
  mat_norm <- mat
  mat_norm <- mat_norm - apply(mat_norm, MARGIN = 1, FUN = min)
  mat_norm <- mat_norm / apply(mat_norm, MARGIN = 1, FUN = max)

  # nmf
  res <- NMF::nmf(mat_norm, rank = K, method = "nsNMF", seed = rep(77, 6), nrun = 30)

  # basis
  w <- NMF::basis(res)
  colnames(w) <- sprintf("CM%02d", 1:ncol(w))

  # df
  w_df <- reshape2::melt(w)
  colnames(w_df) <- c("subCluster", "cm", "weight")

  return(w_df)
}


#' Get Edges of networks
#' @description Get all specifically correlated cell-subset pairs as potential edges of networks
#' @param mat frequency matrix
#' @importFrom dplyr %>%
#' @return data frame of all potential edges of networks.
#' @export
getEdge <- function(mat){
  # 0) cor test
  cor_ls <- t(mat) %>% psych::corr.test()

  # 1) coefficient
  cor_mat <- cor_ls$r
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  cor_df <- reshape2::melt(cor_mat, na.rm = TRUE)
  colnames(cor_df) <- c("node1", "node2", "correlation")

  # 2) pval and fdr
  pval_mat <- cor_ls$p
  pval_mat[lower.tri(pval_mat, diag = TRUE)] <- NA
  pval_df <- reshape2::melt(pval_mat, na.rm = TRUE)
  colnames(pval_df) <- c("node1", "node2", "pval")
  # fdr
  pval_df$pval_fdr <- stats::p.adjust(pval_df$pval, method = "fdr")

  # 3) spe
  cor_mat <- cor_ls$r
  spe_mat <- matrix(NA, nrow(cor_mat), ncol(cor_mat))
  colnames(spe_mat) <- rownames(spe_mat) <- rownames(cor_mat)
  pmt_mat <- spe_mat
  N <- nrow(spe_mat)
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      tmp <- c(cor_mat[i, -i], cor_mat[-c(i, j), j]) %>% unname()
      # spe
      spe <- sum(tmp <= cor_mat[i, j]) / length(tmp)
      spe_mat[i, j] <- spe
      # pmt
      shuf <- dplyr::cume_dist(tmp)
      pmt <- sum(shuf >= spe) / length(shuf)
      # pmt <- sum(tmp >= cor_mat[i, j]) / length(tmp) # alternative
      pmt_mat[i, j] <- pmt
    }
  }
  spe_df <- reshape2::melt(spe_mat, na.rm = TRUE)
  pmt_df <- reshape2::melt(pmt_mat, na.rm = TRUE)
  colnames(spe_df) <- c("node1", "node2", "spe")
  colnames(pmt_df) <- c("node1", "node2", "pmt")

  # merge
  cor_df <- Reduce(f = merge, x = list(cor_df, pval_df, spe_df, pmt_df))
  table(cor_df$spe + cor_df$pmt)

  # filter 1
  cdt1 <- cor_df$correlation > 0.2
  cdt2 <- cor_df$pval_fdr <= 0.05
  cor_df <- cor_df[cdt1 & cdt2, ]

  # filter 2
  cor_df$pmt_fdr <- stats::p.adjust(cor_df$pmt, method = "fdr")
  cdt3 <- cor_df$pmt_fdr <= 0.1
  cor_df <- cor_df[cdt3, ]

  return(cor_df)
}


#' Create and Visualize Multicellular Networks
#' @description Create and visualize multicellular networks
#' @param node data frame of nodes of each network
#' @param edge data frame of edges of networks
#' @param n number of nodes
#' @import dplyr
#' @importFrom igraph E
#' @importFrom igraph V
#' @importFrom igraph layout_in_circle
#' @importFrom graphics par
#' @importFrom graphics title
#' @return plot networks
#' @export
creNet <- function(node, edge, n = 10){
  # filter node
  node <- node %>%
    group_by(cm) %>%
    arrange(desc(weight)) %>%
    slice(1:n) %>%
    ungroup()
  # graph
  graph <- igraph::graph_from_data_frame(edge, directed = FALSE)
  col_fun <- circlize::colorRamp2(
    breaks = stats::quantile(E(graph)$spe, probs = seq(0, 1, length = 6)),
    colors = ggsci::pal_material(palette = "grey", n = 10)(10)[3:8]
  )
  E(graph)$color <- col_fun(E(graph)$spe)
  E(graph)$width <- 1

  # plot
  par(mfrow = c(2, 6), mar = c(0, 0, 0, 0) + 0.5) # 3，4
  for (cm in unique(node$cm)) {
    sub_graph <- igraph::subgraph(
      graph,
      vids = intersect(V(graph)$name, node$subCluster[node$cm == cm])
    )
    sub_graph <- igraph::delete.vertices(sub_graph, V(sub_graph)[igraph::degree(sub_graph) == 0])
    igraph::plot.igraph(
      sub_graph,
      layout = layout_in_circle, # layout_with_fr
      xlim = c(-1.2, 1.2),
      ylim = c(-1.2, 1.2),
      # vertex
      vertex.size = 50,
      vertex.label = V(sub_graph)$short_name,
      vertex.label.cex = 5 / 8,
      vertex.label.color = "black"
      #vertex.label.family = "Arial"
    )
    title(cm, cex.main = 7 / 8, line = -0.5)
  }

}



