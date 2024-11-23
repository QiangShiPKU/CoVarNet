usethis::use_package(package = "RColorBrewer", type = "Imports")
usethis::use_package(package = "circlize", type = "Imports")
usethis::use_package(package = "ggsci", type = "Imports")
usethis::use_package(package = "ComplexHeatmap", type = "Imports")
usethis::use_package(package = "grid", type = "Imports")
usethis::use_package(package = "viridis", type = "Imports")
usethis::use_package(package = "tidytext", type = "Imports")
usethis::use_package(package = "dendextend", type = "Imports")
usethis::use_package(package = "ggplot2", type = "Depends")



#' Visualize the network of each CM
#'
#' @param each The result of function "cm_network".
#' @param Layout The layout methods in igraph (Default is layout_in_circle).
#' @return A graph of cellular network of each CM.
#' @export
#'
#' @examples
#' gr.igraph_each(each)
#' gr.igraph_each(each,Layout=layout_in_circle)
gr.igraph_each<-function(each,...){
  var_args <- list(...)
  Layout <- if(!is.null(var_args[["Layout"]])) var_args[["Layout"]] else layout_in_circle

  node=each$node
  edge=each$edge
  edge=edge[,c("subCluster1","subCluster2","cm","correlation","pval","pval_fdr","spe","majorCluster1","majorCluster2")]
  graph <- graph_from_data_frame(edge, directed = FALSE)

  meta_sc=data.frame(subCluster=c(edge$subCluster1,edge$subCluster2),majorCluster=c(edge$majorCluster1,edge$majorCluster2))
  meta_sc=meta_sc[!duplicated(meta_sc$subCluster),]
  rownames(meta_sc)=meta_sc$subCluster

  # vertex attr
  V(graph)$majorCluster <- meta_sc[V(graph)$name, "majorCluster"]
  colors <- RColorBrewer::brewer.pal(length(unique(V(graph)$majorCluster)), "Set3")
  names(colors)=unique(V(graph)$majorCluster)
  V(graph)$frame.color <- V(graph)$color <- colors[V(graph)$majorCluster]


  # edge attr
  col_fun <- circlize::colorRamp2(
    # breaks = quantile(-log10(E(graph)$pval_fdr), probs = c(0, 0.7)),
    breaks = quantile(E(graph)$spe, probs = seq(0, 1, length = 6)),
    # colors = RColorBrewer::brewer.pal(9, "Reds")[2:9],
    colors = ggsci::pal_material(palette = "grey", n = 10)(10)[3:8]
  )
  # E(graph)$color <- col_fun(-log10(E(graph)$pval_fdr))
  E(graph)$color <- col_fun(E(graph)$spe)
  E(graph)$width <- 1

  par(mfrow = c(2, 6), mar = c(0, 0, 0, 0) + 0.5) # 3，4
  for (cm in unique(node$cm)) {
    sub_graph <- subgraph(
      graph,
      vids = intersect(V(graph)$name, node$subCluster[node$cm == cm])
    )
    sub_graph <- delete.vertices(sub_graph, V(sub_graph)[igraph::degree(sub_graph) == 0])
    plot.igraph(
      sub_graph,
      layout = Layout, # layout_with_fr
      xlim = c(-1.2, 1.2),
      ylim = c(-1.2, 1.2),
      # vertex
      vertex.size = 50,
      vertex.label.cex = 5 / 8,
      vertex.label.color = "black",
      #vertex.label.family = "Arial"
      edge.curved=FALSE
    )
    title(cm, cex.main = 7 / 8, line = -0.5) #, family = "Arial"
  }
}

#' Visualize the global network
#'
#' @param global The result of function "cm_network".
#' @param Layout The layout methods in igraph (Default is layout_in_circle).
#' @return A graph of the global network.
#' @export
#'
#' @examples
#' gr.igraph_global(global)
#' gr.igraph_global(global,Layout=layout_with_fr)
gr.igraph_global<-function(global,...){
  var_args <- list(...)
  Layout <- if(!is.null(var_args[["Layout"]])) var_args[["Layout"]] else layout_with_fr

  edge=global$edge
  node=global$node
  edge=edge[,c("subCluster1","subCluster2","correlation","pval","pval_fdr","spe","majorCluster1","majorCluster2")]
  graph <- graph_from_data_frame(edge, directed = FALSE)

  # vertex attr
  rownames(node)=node$subCluster
  V(graph)$majorCluster <- node[V(graph)$name, "majorCluster"]

  colors <- RColorBrewer::brewer.pal(length(unique(V(graph)$majorCluster)), "Set3")
  names(colors)=unique(V(graph)$majorCluster)
  V(graph)$frame.color <- V(graph)$color <- colors[V(graph)$majorCluster]

  # edge attr
  col_fun <- circlize::colorRamp2(
    # breaks = quantile(-log10(E(graph)$pval_fdr), probs = c(0, 0.7)),
    breaks = quantile(E(graph)$spe, probs = seq(0, 1, length = 6)),
    # colors = RColorBrewer::brewer.pal(9, "Reds")[2:9],
    colors = ggsci::pal_material(palette = "grey", n = 10)(10)[3:8]
  )
  # E(graph)$color <- col_fun(-log10(E(graph)$pval_fdr))
  E(graph)$color <- col_fun(E(graph)$spe)
  E(graph)$width <- 1

  lgd <- ComplexHeatmap::Legend(
    col_fun = col_fun,
    at = c(min(E(graph)$spe), max(E(graph)$spe)),
    labels = c("Low", "High"),
    title = expression("Specificity"),
    direction = "horizontal",
    grid_height = unit(3, "mm"),
    legend_width = unit(10, "mm"),
    labels_gp = grid::gpar(fontsize = 6),
    title_gp = grid::gpar(fontsize = 6),
    title_position = c("topcenter")
  )
  plot.igraph(
    graph,
    layout = Layout,
    # vertex
    vertex.size = 6,
    vertex.label.cex = 6 / 12,
    vertex.label.color = "black",
    # vertex.label.family = "Arial"
  )
  ComplexHeatmap::draw(lgd, x = unit(0.1, "npc"), y = unit(0.05, "npc"))
}



#' Visualizing the distribution of CMs across samples
#' @param nmf_res The result of function nmf which includes the coefficient matrix.
#' @param meta The groupings of each sample(not required). If provided, this file should contain at least two columns:"sampleID" and the groupings category.
#' @param group The column name for groupings in meta.
#' @return The distribution of CMs across samples.
#' @export
#'
#' @examples
#' gr.distribution(nmf_res)
#' gr.distribution(nmf_res,meta=meta,group="tissue")
gr.distribution<-function(nmf_res,...){
  var_args <- list(...)
  meta <- if(!is.null(var_args[["meta"]])) var_args[["meta"]] else NULL
  group <- if(!is.null(var_args[["group"]])) var_args[["group"]] else NULL

  h=scoef(nmf_res) #得到nmf的coeffcient matrix
  if(is.null(rownames(h))){
    rownames(h) <- sprintf("CM%02d", 1:ncol(w))
  }
  sorted_rownames <- order(rownames(h))
  h <- h[sorted_rownames,]

  id <- apply(h, MARGIN = 2, FUN = which.max)

  if(is.null(meta)){
    meta_sp<-as.data.frame(sprintf("CMT%02d", 1:nrow(h))[id])
    rownames(meta_sp)<-colnames(h)
    colnames(meta_sp)<-"CMT"
    meta_sp["sampleID"]<-rownames(meta_sp)
    meta_sp<-meta_sp[order(meta_sp$CMT),]
  }
  else{
    meta_sp <- meta[!duplicated(meta$sampleID), ]
    rownames(meta_sp) <- meta_sp$sampleID
    meta_sp<-meta_sp[colnames(h),]
    meta_sp$CMT <- sprintf("CMT%02d", 1:nrow(h))[id]
    if(is.null(group)){
      meta_sp<-meta_sp[order(meta_sp$CMT),]
    }
    else{
      meta_sp <- meta_sp[
        order(meta_sp$CMT, meta_sp[,group]),
      ]
    }
  }

  h=h[,rownames(meta_sp)]

  CMTnames<-gsub("CM", "CMT",rownames(h))
  colors <- RColorBrewer::brewer.pal(nrow(h), "Set3")
  colors=colors[1:nrow(h)]
  names(colors)=CMTnames

  if(is.null(group)){
    ann_top <- ComplexHeatmap::HeatmapAnnotation(
      CMT = meta_sp$CMT,
      # color
      col=list(CMT=colors),
      #col = list(CMT = color$CMT,Tissue = color$tissue),
      border = TRUE, # show_legend = c(FALSE, FALSE, TRUE),
      annotation_name_side = "left",
      annotation_name_gp = grid::gpar(fontsize = 6),
      annotation_legend_param = list(
        legend_direction = "horizontal", ncol = 2,
        labels_gp = grid::gpar(fontsize = 6),
        title_gp = grid::gpar(fontsize = 6),
        grid_width = unit(3, "mm"),
        grid_height = unit(2, "mm")
      )
    )
  }
  else{
    ann_top <- ComplexHeatmap::HeatmapAnnotation(
      CMT = meta_sp$CMT,
      Group = meta_sp[,group],
      # color
      col=list(CMT=colors),
      #col = list(CMT = color$CMT,Tissue = color$tissue),
      border = TRUE, # show_legend = c(FALSE, FALSE, TRUE),
      annotation_name_side = "left",
      annotation_name_gp = grid::gpar(fontsize = 6),
      annotation_legend_param = list(
        legend_direction = "horizontal", ncol = 2,
        labels_gp = grid::gpar(fontsize = 6),
        title_gp = grid::gpar(fontsize = 6),
        grid_width = unit(3, "mm"),
        grid_height = unit(2, "mm")
      )
    )
  }
  ComplexHeatmap::Heatmap(h,
                          width = unit(3, "in"), height = unit(1.5, "in"),
                          name = "CM abundance", border = TRUE,
                          show_column_names = FALSE, row_names_side = "left",
                          cluster_rows = FALSE, cluster_columns = FALSE,
                          row_names_gp = grid::gpar(fontsize = 6),
                          column_title = "CM abundance in samples",
                          column_title_gp = grid::gpar(fontsize = 8),
                          col = circlize::colorRamp2(seq(0, 0.5, length.out = 4), viridis::viridis(4)),
                          top_annotation = ann_top,
                          heatmap_legend_param = list(
                            legend_direction = "horizontal", color_bar = "continuous",
                            title_position = "lefttop",
                            labels_gp = grid::gpar(fontsize = 6), title_gp = grid::gpar(fontsize = 6),
                            grid_height = unit(3, "mm")
                          )
  ) -> ht

  ComplexHeatmap::draw(ht,
                       heatmap_legend_side = "bottom", annotation_legend_side = "right",
                       merge_legends = FALSE, # adjust_annotation_extension = TRUE,
  )
  dup <- (which(!duplicated(meta_sp$CMT)) - 1)
  fract <- dup / nrow(meta_sp)
  width <- c(fract[-1], 1) - fract
  ComplexHeatmap::decorate_heatmap_body("CM abundance", {
    grid::grid.rect(
      x = unit(fract, "native"), y = unit(1 - (0:(length(dup) - 1)) / length(dup), "native"),
      width = unit(width, "native"), height = unit(1 / length(unique(meta_sp$CMT)), "native"),
      hjust = 0, vjust = 1, gp = grid::gpar(col = "white", fill = NA, lty = 1, lwd = 2)
    )
  })



}




#' Visualize the weights of all cell subsets in CMs
#' @param nmf_res The result of function nmf, which includes the feature matrix.
#' @return A heatmap displaying the weights of all cell subsets within each cellular module.
#' @export
#'
#' @examples
#' gr.weight_all(nmf_res)
gr.weight_all<-function(nmf_res){

  w <- basis(nmf_res)
  if(is.null(colnames(w))){
    colnames(w) <- sprintf("CM%02d", 1:ncol(w))
  }

  ComplexHeatmap::Heatmap(
    matrix = w, name = "weight",
    # col = viridis(n = 5),
    circlize::colorRamp2(seq(0, 0.1, length.out = 5), viridis::viridis(n = 5)),
    width = unit(5, "cm"), height = unit(10, "cm"),
    row_names_gp = grid::gpar(fontsize = 4.5), column_names_gp = grid::gpar(fontsize = 4.5),
    clustering_method_rows = "ward.D2"
  ) -> ht1

  ComplexHeatmap::draw(ht1)
}



#' Visualize cell subsets with the highest weights in each CM
#' @param nmf_res The result of function nmf, which includes the feature matrix.
#' @param num The number of cell subsets with the highest weights in each CM.
#' @return The top cell subsets sorted by NMF weights for each cellular module.
#' @export
#'
#' @examples
#' gr.weight_top(nmf_res,num=15)
gr.weight_top<-function(nmf_res,...){
  var_args <- list(...)
  num <- if(!is.null(var_args[["num"]])) var_args[["num"]] else 15

  w <- basis(nmf_res)
  if(is.null(colnames(w))){
    colnames(w) <- sprintf("CM%02d", 1:ncol(w))
  }

  sorted_colnames <- order(colnames(w))
  w <- w[, sorted_colnames]

  w_df <- reshape2::melt(w)
  head(w_df)
  colnames(w_df) <- c("subCluster", "cm", "weight")

  w_df <- w_df %>%
    group_by(cm) %>%
    arrange(desc(weight), .by_group = TRUE) %>%
    top_n(num, weight) %>%
    as.data.frame()

  ggplot(w_df, aes(x = tidytext::reorder_within(subCluster, weight, cm), y = weight, fill = cm)) +
    geom_col() +
    # scale_fill_manual(values = color$CM) +
    # facet_grid(cm ~ ., scales = "free") +
    facet_wrap(~cm, scales = "free", nrow = 3) +
    coord_flip() +
    # coord_fixed(ratio = 10 / 1) +
    labs(x = "", y = "Weight") +
    tidytext::scale_x_reordered() +
    # scale_x_discrete(limits = rev(levels(w_df$subCluster))) +
    theme_classic() +
    theme(
      text = element_text(size = 6, color = "black"),
      axis.text = element_text(size = 6, color = "black"),
      legend.key.size = unit(3, "mm"),
      legend.position = "none",
      # plot.title = element_text(size = 8, hjust = 0.5),
      strip.text.x = element_text(size = 7),
      strip.background = element_blank(),
      aspect.ratio = 1.6 / 1
    )
}




#' Calculate the spatial proximity between cell subset pairs within several specified cellular modules.
#'
#' @param coords The coordinates of each spot.
#' @param st_fq The frequencies of each cell subset in each spot.
#' @param ref The "ref$filter" object enables us to obtain the components of each cellular module.
#' @param cm_select A character containing the names of specified cellular modules. For example: c("CM01","CM02").
#'
#' @return A heatmap displaying the spatial proximity between cell subset pairs.
#' @export
#'
#' @examples
#' gr.proximity(coords,st_fq,c("CM01","CM02"))
gr.proximity<-function(coords,st_fq,ref,cm_select){
  nb <- spdep::knearneigh(coords, k = 60) ### K个最近邻居（参与影响细胞丰度的邻居）
  nb <- spdep::knn2nb(nb)
  weight <- spdep::nb2listwdist(
    neighbours = nb,
    x = sp::SpatialPoints(coords),
    type = "idw", style = "W", alpha = 2
  )
  # type:the intended type of distance modelling. idw反比例，次数为alpha=2.
  # style对weight进行标准化，W为 row standardised，对X的所有邻居标准化到和为1

  sc_ls=colnames(st_fq)
  DM <- matrix(nrow = ncol(st_fq), ncol = ncol(st_fq))
  rownames(DM) <- colnames(DM) <- sc_ls
  diag(DM) <- 0
  for (i in sc_ls) {
    for (j in sc_ls) {
      if (i != j) {
        dimoran <- spdep::moran_bv(
          x = st_fq[, i],
          y = st_fq[, j],
          listw = weight, nsim = 2
        )
        DM[i, j] <- dimoran$t0
      }
    }
  }

  node<-ref$filter
  node <- node[node$cm %in% cm_select, ]

  colors <- RColorBrewer::brewer.pal(length(unique(node$cm)), "Set3")
  colors=colors[1:length(unique(node$cm))]
  names(colors)=unique(node$cm)

  ann_top <- ComplexHeatmap::HeatmapAnnotation(
    CM = node$cm,
    col=list(CM=colors),
    simple_anno_size = unit(3, "mm"),
    border = TRUE, annotation_name_gp = grid::gpar(fontsize = 6),
    annotation_legend_param = list(
      labels_gp = grid::gpar(fontsize = 6),
      title_gp = grid::gpar(fontsize = 6),
      grid_width = unit(3, "mm"),
      grid_height = unit(2, "mm")
    )
  )
  ann_left <- ComplexHeatmap::rowAnnotation(
    CM = node$cm,
    col=list(CM=colors),
    show_legend=FALSE,
    simple_anno_size = unit(3, "mm"),
    border = TRUE, annotation_name_gp = grid::gpar(fontsize = 6),
    annotation_legend_param = list(
      labels_gp = grid::gpar(fontsize = 6),
      title_gp = grid::gpar(fontsize = 6),
      grid_width = unit(3, "mm"),
      grid_height = unit(2, "mm")
    )
  )


  ileum <- DM[node$subCluster, node$subCluster]
  diag(ileum) <- 1
  ileum[ileum == 0] <- 1
  ComplexHeatmap::Heatmap(
    matrix = ileum,
    col = circlize::colorRamp2(seq(0, quantile(ileum, 0.7) %>% unname(), length.out = 4), RColorBrewer::brewer.pal(n = 7, name = "RdBu")[4:7]),
    border = TRUE,
    width = unit(6, "cm"), height = unit(6, "cm"),
    column_names_rot = 90, top_annotation = ann_top, left_annotation = ann_left,
    row_names_gp = grid::gpar(fontsize = 6), column_names_gp = grid::gpar(fontsize = 6),
    column_title_gp = grid::gpar(fontsize = 8),
    column_dend_height = unit(5, "mm"), row_dend_width = unit(5, "mm"),
    # clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
    heatmap_legend_param = list(
      # at = seq(0, 3.5, 0.5),
      # labels = c("0", "0.5", "1", "(+)", "2", "(++)", "3", "(+++)"),
      title = "Moran's I",
      # legend_direction = "horizontal",
      grid_width = unit(3, "mm"), # legend_height = unit(1, "cm"),
      labels_gp = grid::gpar(fontsize = 6), title_gp = grid::gpar(fontsize = 6)
      # title_position = "topcenter"
    )
  )->p1
  ComplexHeatmap::draw(p1)
}




