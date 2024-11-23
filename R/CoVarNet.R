

usethis::use_package(package = "NMF", type = "Depends")
usethis::use_package(package = "dplyr", type = "Depends")
usethis::use_package(package = "cluster", type = "Imports")
usethis::use_package(package = "psych", type = "Imports")
usethis::use_package(package = "igraph", type = "Depends")
usethis::use_package(package = "Seurat", type = "Depends")
usethis::use_package(package = "sp", type = "Imports")
usethis::use_package(package = "spdep", type = "Imports")




#' Calculate raw cell subset frequency matrix
#' @description Calculate frequency matrix from cell annotation dataframe.
#' @param meta Cell annotation dataframe (The dataframe should contain at least three columns: "sampleID", "majorCluster","subCluster").
#'
#' @return Raw cell subsets frequency matrix.
#' @export
#'
#' @examples
#' freq_calculate(meta)
freq_calculate<-function(meta){
  meta_sc <- meta[!duplicated(meta$subCluster), c("majorCluster", "subCluster")]
  rownames(meta_sc) <- meta_sc$subCluster
  meta_sc <- meta_sc[order(meta_sc$subCluster), ]

  mat_ct <- table(meta$subCluster, meta$sampleID)
  mat_ct <- matrix(
    data = mat_ct,
    ncol = ncol(mat_ct),
    dimnames = dimnames(mat_ct)
  ) %>% as.data.frame()

  ### add majorCluster
  id <- match(rownames(mat_ct), meta_sc$subCluster)
  mat_ct <- data.frame(
    majorCluster = factor(meta_sc$majorCluster[id]),
    mat_ct,
    check.names = FALSE
  )

  ### frequency matrix (scale)
  mat_fq <- mat_ct %>%
    group_by(majorCluster) %>%
    mutate_if(is.numeric, .funs = list(~ . / sum(.))) %>%
    as.data.frame()
  major=mat_fq$majorCluster
  mat_fq[is.na(mat_fq)] <- 0
  rownames(mat_fq) <- rownames(mat_ct)

  return(mat_fq)
}

#' Normalize cell subset frequency matrix
#'
#' @param mat_fq Raw cell subset frequency matrix.
#' @param normalize Normalize methods:"minmax" represents Min-Max Normalization,"zscore" represents Z-Score Normalization.
#'
#' @return Normalized cell subset frequency matrix.
#' @export
#'
#' @examples
#' freq_normalize(mat_fq)
#' freq_normalize(mat_fq,normalize="minmax")
freq_normalize<-function(mat_fq,...){
  var_args <- list(...)
  normalize <- if(!is.null(var_args[["normalize"]])) var_args[["normalize"]] else "minmax"

  mat_fq<-mat_fq[,-1]
  cluster=rownames(mat_fq)
  mat_fq<-sapply(mat_fq,as.numeric)
  rownames(mat_fq)<-cluster

  if(normalize=="minmax"){
    mat_normalize <- mat_fq
    mat_normalize <- mat_normalize - apply(mat_normalize, MARGIN = 1, FUN = min)
    mat_normalize <- mat_normalize / apply(mat_normalize, MARGIN = 1, FUN = max)
  }
  else if(normalize=="zscore"){
    mat_normalize <- mat_fq
    mat_normalize <- mat_normalize - apply(mat_normalize, MARGIN = 1, FUN = mean)
    mat_normalize <- mat_normalize / apply(mat_normalize, MARGIN = 1, FUN = sd)
    mat_normalize[mat_normalize < 0] <- 0
  }
  return(mat_normalize)
}

#' Normalize cell subset frequency matrix in spatial data
#'
#' @param spe A seurat object(spatial transcriptomics data). The meta data should contain annotation information for each subCluster and be consistent with the reference, representing the abundance of each cell subset in each spot.
#' @param ref The "ref$ann" object contains all the subCluster and corresponding majorCluster information, so that we can compute the frequencies of cell subsets within corresponding broad cell types.
#' @param normalize Normalize methods:"minmax" represents Min-Max Normalization,"zscore" represents Z-Score Normalization.
#'
#' @return A seurat object(normalized spatial transcriptomics data).
#' @export
#'
#' @examples freq_normalize_st(spe,ref,normalize="minmax")
freq_normalize_st<-function(spe,ref,...){
  var_args <- list(...)
  normalize <- if(!is.null(var_args[["normalize"]])) var_args[["normalize"]] else "minmax"

  # add majorCluster
  mat_ct<-t(spe@meta.data[,rownames(ref$ann)])
  id <- match(rownames(mat_ct), ref$ann$subCluster)
  mat_ct <- data.frame(
    majorCluster = factor(ref$ann$majorCluster[id]),
    mat_ct,
    check.names = FALSE
  )

  ### frequency matrix (scale)
  mat_fq <- mat_ct %>%
    group_by(majorCluster) %>%
    mutate_if(is.numeric, .funs = list(~ . / sum(.))) %>%
    as.data.frame()
  major=mat_fq$majorCluster
  mat_fq[is.na(mat_fq)] <- 0
  rownames(mat_fq) <- rownames(mat_ct)

  mat_fq<-mat_fq[,-1]
  cluster=rownames(mat_fq)
  mat_fq<-sapply(mat_fq,as.numeric)
  rownames(mat_fq)<-cluster

  if(normalize=="minmax"){
    mat_normalize <- mat_fq
    mat_normalize <- mat_normalize - apply(mat_normalize, MARGIN = 1, FUN = min)
    mat_normalize <- mat_normalize / apply(mat_normalize, MARGIN = 1, FUN = max)
  }
  else if(normalize=="zscore"){
    mat_normalize <- mat_fq
    mat_normalize <- mat_normalize - apply(mat_normalize, MARGIN = 1, FUN = mean)
    mat_normalize <- mat_normalize / apply(mat_normalize, MARGIN = 1, FUN = sd)
    mat_normalize[mat_normalize < 0] <- 0
  }
  return(mat_normalize)
}


#' Calculate correlation coefficients between cell subsets
#'
#' @param mat_fq Raw cell subset frequency matrix.
#' @param method Methods for calculating correlation:"pearson" or "spearman". Default is "pearson".
#'
#' @return A data frame containing correlation between all pairs of cell subsets, including the following columns: "subCluster1", "subCluster2", "majorCluster1", "majorCluster2", "correlation", "pval", "pval_fdr", "spe". "spe" is the specificity of subset pairs. For each element 𝑟,where i<j in cor_pair, its background (complementary) set is defined as {rik,rkj|k!=i,j}. The specificity is defined as the fraction of elements in the background set that not exceeding rij.
#' @export
#'
#' @examples
#' pair_correlation(mat_fq)
#' pair_correlation(mat_fq,method="pearson")
pair_correlation<-function(mat_fq,...){
  var_args <- list(...)
  method <- if(!is.null(var_args[["method"]])) var_args[["method"]] else "pearson"

  meta_sc <- as.data.frame(mat_fq[,1])
  rownames(meta_sc) <- rownames(mat_fq)
  meta_sc$subCluster=rownames(meta_sc)
  colnames(meta_sc)<-c("majorCluster","subCluster")

  mat_fq<-mat_fq[,-1]
  cluster=rownames(mat_fq)
  mat_fq<-sapply(mat_fq,as.numeric)
  rownames(mat_fq)<-cluster

  # 1) correlation
  cor_ls <- t(mat_fq) %>% psych::corr.test(method = method)
  cor_mat <- cor_ls$r
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  cor_df <- reshape2::melt(cor_mat, na.rm = TRUE)
  colnames(cor_df) <- c("subCluster1", "subCluster2", "correlation")

  # 2) pval and fdr
  pval_mat <- cor_ls$p
  pval_mat[lower.tri(pval_mat, diag = TRUE)] <- NA
  pval_df <- reshape2::melt(pval_mat, na.rm = TRUE)
  colnames(pval_df) <- c("subCluster1", "subCluster2", "pval")
  # fdr
  pval_df$pval_fdr <- p.adjust(pval_df$pval, method = "fdr")

  # 3) specificity
  cor_mat <- cor_ls$r
  spe_mat <- matrix(NA, nrow(cor_mat), ncol(cor_mat))
  colnames(spe_mat) <- rownames(spe_mat) <- rownames(cor_mat)
  N <- nrow(spe_mat)
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      bg <- c(cor_mat[i, -i], cor_mat[-c(i, j), j]) %>% unname()
      # specificity
      spe_mat[i, j] <- mean(bg <= cor_mat[i, j])
    }
  }
  spe_df <- reshape2::melt(spe_mat, na.rm = TRUE)
  colnames(spe_df) <- c("subCluster1", "subCluster2", "spe")

  # merge
  cor_df <- Reduce(f = merge, x = list(cor_df, pval_df, spe_df))
  cor_df$majorCluster1 <- meta_sc$majorCluster[match(cor_df$subCluster1, meta_sc$subCluster)]
  cor_df$majorCluster2 <- meta_sc$majorCluster[match(cor_df$subCluster2, meta_sc$subCluster)]
  return(cor_df)
}



#' Construct cellular network
#' @description Interconnect co-occurring nodes via specifically correlated edges to construct CM networks.
#' @param NMFres Co-occurring nodes in each CM obtained from NMF results.
#' @param Corres Specifically correlated edges obtained from correlation results.
#' @param top_n Select the upper limit of the number of cell subsets in each CM.
#' @param corr Filter specifically correlated edges with correlation greater than corr.
#' @param fdr Filter specifically correlated edges with pval_fdr less than fdr.
#' @return A list contains the following parts."global": the global network, "each": the network of each CM, "raw": NMF result, "filter": the components of each CM, "ann": the correspondence between subClusters and majorClusters. It can be used to draw networks with the functions "gr.igraph_each" and "gr.igraph_global", export the components of CMs, or serve as a reference to recover cellular modules in Tutorial 2.

#' @export
#'
#' @examples
#' cm_network(NMF_K12,cor_pair)
#' cm_network(NMF_K12,cor_pair,top_n=10,corr=0.2,fdr=0.05)
cm_network<-function(NMFres,Corres,...){
  var_args <- list(...)
  top_n <- if(!is.null(var_args[["top_n"]])) var_args[["top_n"]] else 10
  corr <- if(!is.null(var_args[["corr"]])) var_args[["corr"]] else 0.2
  fdr <- if(!is.null(var_args[["pval_fdr"]])) var_args[["pval_fdr"]] else 0.05

  listA=Corres[,c("subCluster1","majorCluster1")]
  listB=Corres[,c("subCluster2","majorCluster2")]
  colnames(listA)=c("subCluster","majorCluster")
  colnames(listB)=c("subCluster","majorCluster")
  ann=rbind(listA,listB)
  ann=ann[!duplicated(ann$subCluster), ]
  rownames(ann)=ann$subCluster


  w <- basis(NMFres)
  if(is.null(colnames(w))){
    colnames(w) <- sprintf("CM%02d", 1:ncol(w))
    rownames(coef(NMF_K12))<- sprintf("CM%02d", 1:ncol(w))
    colnames(basis(NMF_K12))<- sprintf("CM%02d", 1:ncol(w))
  }

  sorted_colnames <- order(colnames(w))
  w <- w[, sorted_colnames]

  weight <- reshape2::melt(w)
  colnames(weight) <- c("subCluster", "cm", "weight")
  # cumsum
  weight <- weight %>%
    group_by(cm) %>%
    arrange(cm, desc(weight)) %>%
    ungroup()
  # top_n subsets of each CM
  weight <- weight %>%
    group_by(cm) %>%
    arrange(desc(weight)) %>%
    slice(1:top_n) %>%
    ungroup()


  meta_sc=data.frame(subCluster=c(Corres$subCluster1,Corres$subCluster2),majorCluster=c(Corres$majorCluster1,Corres$majorCluster2))
  meta_sc <- meta_sc[!duplicated(meta_sc$subCluster),]
  rownames(meta_sc) <- meta_sc$subCluster
  meta_sc <- meta_sc[order(meta_sc$subCluster), ]

  ###########################
  # (1) network
  ########################
  # filter 1
  cdt1 <- Corres$correlation > corr
  cdt2 <- Corres$pval_fdr <= fdr
  pl_df <- Corres[cdt1 & cdt2, ]

  # filter 2
  n_all <- nrow(meta_sc) # subset number
  spe_cutoff <- 1 - ((top_n - 1) * 2 - 1) / ((n_all - 1) * 2 - 1)
  cdt3 <- pl_df$spe >= spe_cutoff
  pl_df <- pl_df[cdt3, ]


  graph <- graph_from_data_frame(pl_df, directed = FALSE)

  # global
  node_global <- V(graph)$name
  node_global <- data.frame(
    subCluster = node_global,
    majorCluster=meta_sc$majorCluster[match(node_global, meta_sc$subCluster)]
  )
  node_global <- node_global[order(node_global$subCluster), ]

  # each
  node_each <- matrix(NA, nrow = 0, ncol = 2) %>% as.data.frame()
  colnames(node_each) <- c("cm", "subCluster")
  edge_each <- matrix(NA, nrow = 0, ncol = ncol(pl_df) + 1) %>% as.data.frame()
  colnames(edge_each) <- c("cm", colnames(pl_df))
  for (cm in unique(weight$cm)) {
    sub_graph <- subgraph(
      graph,
      vids = intersect(V(graph)$name, weight$subCluster[weight$cm == cm])
    )
    sub_graph <- delete.vertices(sub_graph, V(sub_graph)[degree(sub_graph) == 0])
    if(length(V(sub_graph))!=0){
      tmp1 <- data.frame(cm = cm, subCluster = V(sub_graph)$name)
      node_each <- rbind(node_each, tmp1)
      tmp2 <- pl_df[(pl_df$subCluster1 %in% tmp1$subCluster) & (pl_df$subCluster2 %in% tmp1$subCluster), ]
      edge_each <- rbind(edge_each, data.frame(cm = cm, tmp2))
    }
  }
  node_each$majorCluster <- meta_sc$majorCluster[match(node_each$subCluster, meta_sc$subCluster)]
  rownames(edge_each) <- NULL

  weight=merge(node_each,weight)
  weight=weight[order(weight$cm,-weight$weight),]
  rownames(weight)=as.list(1:nrow(weight))

  cm_network<-list(list(node_global,pl_df),list(node_each,edge_each),NMFres,weight,ann)
  names(cm_network)<-c("global","each","raw","filter","ann")
  names(cm_network$global)<-c("node","edge")
  names(cm_network$each)<-c("node","edge")
  return(cm_network)
}



#' Recover cellular modulrs in scRNA_seq data
#' @description Utilize the reference file to recover the sample CMT in scRNA_seq data.
#' @param ref The parameter "ref" is the result of the function cm_network, which serves as the reference file of cellular modules recovery.
#' @param mat_fq Normalized sample-cell frequency matrix.
#'
#' @return A matrix containing the abundance of the CMs in each sample.
#' @export
#'
#' @examples
#' sc_cm_recover(ref,mat_fq)
sc_cm_recover<-function(ref,mat_fq){
  W=data.frame(ref$raw@fit@W)
  NMFpredict <- function(W, new_data){
    to_predict=as.matrix(new_data)
    my_method <- function (i, v, x, copy = FALSE, eps = .Machine$double.eps, ...)
    {
      w <- .basis(x)
      h <- .coef(x)
      nb <- nbterms(x)
      nc <- ncterms(x)
      h <- NMF:::std.divergence.update.h(v, w, h, nbterms = nb, ncterms = nc, copy = copy)
      #w <- NMF:::std.divergence.update.w(v, w, h, nbterms = nb, ncterms = nc, copy = copy)
      if (i%%10 == 0) {
        h <- pmax.inplace(h, eps, icterms(x))
        #w <- pmax.inplace(w, eps, ibterms(x))
      }
      if (copy) {
        #.basis(x) <- w
        .coef(x) <- h
      }
      return(x)
    }

    ws = W
    ws <- ws[apply(to_predict, 1, function(x) var(x) > 0),]
    to_predict = to_predict[apply(to_predict, 1, function(x) var(x) > 0),]
    ws = as.matrix(ws)

    dummy = rnmf(ncol(W), to_predict)

    my.seeding.method <- function(model, target){
      basis(model) <- ws #estim.r@fit@W
      # initialize H randomly
      coef(model) <- dummy@H
      # return updated object
      return(model)
    }

    nmf_method <- NMFStrategy('my-method', 'brunet', Update = my_method, objective = 'KL', Stop='connectivity')

    new_nmf = nmf(to_predict, ncol(W), nrun = 1, method = nmf_method, seed = my.seeding.method, .opt='P1')
    return(new_nmf)
  }
  new_nmf<-NMFpredict(W,mat_fq)
  h<-new_nmf@fit@H

  rownames(h)<-colnames(W)
  sorted_rownames <- order(rownames(h))
  h <- h[sorted_rownames,]

  return(h)
}


#' Recover cellular modules in spatial transcriptomics data
#'
#' @param weight The weights of the cell subsets that make up each cellular module.
#' @param spe A seurat object(normalized spatial transcriptomics data) containing the frequencies of each cell subset in each spot.
#'
#' @return A seurat object containg the abundance of each cellular module in each spot.
#' @export
#'
#' @examples st_cm_recover(weight,spe)
st_cm_recover<-function(weight,spe){
  column_names = unique(weight$cm)

  for (module in column_names){
    wgt <- weight %>% filter(cm == module)
    mat1 <- spe@meta.data[,wgt$subCluster]
    mat2 <- sweep(mat1, 2, wgt$weight, '*')
    spe=AddMetaData(object = spe,
                    metadata = rowSums(mat2),
                    col.name = module)
  }

  # 归一化spe.obs中的指定列
  values <- spe@meta.data[, column_names]
  q99 <- quantile(reshape2::melt(values)[,"value"], 0.99)
  spe@meta.data[,column_names]=spe@meta.data[,column_names]/q99
  return(spe)
}

#' Compute colocalization score of each cellular module based on Spearman correlation coefficients
#'
#' @param st_fq The frequencies of each cell subset in each spot.
#' @param ref The "ref$filter" object enables us to obtain the components of each cellular module.
#' @param method Methods for calculating correlation:"pearson" and "spearman". Default is "spearman".
#'
#' @return A data frame containing the colocalization score of each cellular module.
#' @export
#'
#' @examples
#' colocalization(st_fq,ref,method="spearman")
colocalization<-function(st_fq,ref,...){
  var_args <- list(...)
  method <- if(!is.null(var_args[["method"]])) var_args[["method"]] else "spearman"

  node=ref$filter
  cor_mat <- cor(st_fq, method = method)
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  cor_df <- reshape2::melt(cor_mat, na.rm = TRUE)
  colnames(cor_df) <- c("subCluster1", "subCluster2", "correlation")

  stat <- data.frame(matrix(NA, nrow = 0, ncol = 3))
  for (cm in unique(node$cm)) {
    sc <- node$subCluster[node$cm == cm]
    cdt1 <- cor_df$subCluster1 %in% sc
    cdt2 <- cor_df$subCluster2 %in% sc
    cor_inner <- median(cor_df$correlation[cdt1 & cdt2])
    bg <- cor_df$correlation[cdt1 | cdt2]
    score <- mean(cor_inner >= bg)
    stat <- rbind(stat, data.frame(cm, score))
  }
  colnames(stat) <- c("cm", "score")
  return(stat)
}


#' Compute aggregation score of each cellular module based on Global Bivariate Moran's I
#' @param coords The coordinates of each spot
#' @param st_fq The frequencies of each cell subset in each spot.
#' @param ref The "ref$filter" object enables us to obtain the components of each cellular module.
#'
#' @return A data frame containing the aggregation score of each cellular module.
#' @export
#'
#' @examples
#' bvMoranI(coords,st_fq,ref)
bvMoranI<-function(coords,st_fq,ref){
  node=ref$filter

  neighbors=60
  type="idw"
  style="W"
  alpha=2

  nb <- spdep::knearneigh(coords, k = neighbors) ### K个最近邻居（参与影响细胞丰度的邻居）
  nb <- spdep::knn2nb(nb)

  weight <- spdep::nb2listwdist(
    neighbours = nb,
    x = sp::SpatialPoints(coords),
    type = type, style = style, alpha = alpha
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

  cor_mat <- as.matrix(DM)
  diag(cor_mat) <- NA
  cor_df <- reshape2::melt(cor_mat, na.rm = TRUE)
  colnames(cor_df) <- c("subCluster1", "subCluster2", "MoranI")

  stat <- data.frame(matrix(NA, nrow = 0, ncol = 3))
  for (cm in unique(node$cm)) {
    sc <- node$subCluster[node$cm == cm]
    cdt1 <- cor_df$subCluster1 %in% sc
    cdt2 <- cor_df$subCluster2 %in% sc
    cor_inner <- median(cor_df$MoranI[cdt1 & cdt2])
    bg <- cor_df$MoranI[cdt1 | cdt2]
    score <- mean(cor_inner >= bg)
    stat <- rbind(stat, data.frame(cm, score))
  }
  colnames(stat) <- c("cm", "score")
  return(stat)
}
