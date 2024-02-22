library(Seurat)
library(Libra)
library(rlang)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(scDblFinder)

#functions for RunEdgeRPseudobulk
library(magrittr, include.only = c("%<>%"))
library(tibble, include.only=c("rownames_to_column"))
library(edgeR, include.only=c("DGEList","calcNormFactors", "estimateDisp",
                              "glmQLFit","glmQLFTest","glmFit","glmLRT", 
                              "topTags", "cpm"))

library(purrr, include.only=c("map"))
library(stats, include.only=c("model.matrix"))
library(forcats, include.only=c("fct_recode"))
# Modify it to handle multiome
Load10xData <- function(data_dir, sample_names, sub_rna_dir = "filtered_feature_bc_matrix", with.rna.only=F,with.multiome.rna.only=F,with.multiome=F,sub_tcr_dir = "TCR", human_only = F){
  data_S_list <- list()
  data_dir0 <- data_dir
  if (with.rna.only){
    for(i in sample_names){
      data_dir <- paste(data_dir0, i, sub_rna_dir, sep = "/")
      print(data_dir)
      data_S_list[[i]] <- CreateSeuratObject(Read10X(data.dir = data_dir), project = i)
    }
  }
  if (with.multiome){
    # For output from CellRanger >= 3.0 with multiple data types
    #list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
    data_dir0 <- data_dir
    for(i in sample_names){
      data_dir <- paste(data_dir0, i, sub_rna_dir, sep = "/")
      print(data_dir)
      data <- Read10X(data.dir = data_dir)
      data_S_list[[i]] <- CreateSeuratObject(counts = data$`Gene Expression`, project = i)
      data_S_list[[i]][["Peaks"]] <- CreateAssayObject(counts = data$`Peaks`)
    }
  }
  if (with.multiome.rna.only){
    # For output from CellRanger >= 3.0 with multiple data types
    #list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
    data_dir0 <- data_dir
    for(i in sample_names){
      data_dir <- paste(data_dir0, i, sub_rna_dir, sep = "/")
      print(data_dir)
      data <- Read10X(data.dir = data_dir)
      data_S_list[[i]] <- CreateSeuratObject(counts = data$`Gene Expression`, project = i)
    }
  }
  if (human_only){
    for(i in names(data_S_list)){
      data_S_list[[i]] <- data_S_list[[i]][grep("^hg19", rownames(data_S_list[[i]])), ]
    }
  }
  if (T){
    if (length(grep("^mt",rownames(data_S_list[[1]]))) > 0) pattern_mt <- "^mt"
    if (length(grep("^MT",rownames(data_S_list[[1]]))) > 0) pattern_mt <- "^MT"
    if (length(grep("^hg19-MT",rownames(data_S_list[[1]]))) > 0) pattern_mt <- "^hg19-MT"
    print(pattern_mt)
    for(i in names(data_S_list)){
      data_S_list[[i]][["percent.mt"]] <- PercentageFeatureSet(
        object = data_S_list[[i]], pattern = pattern_mt
      )
    }
  }
  
  data_S_list
}


# 
get_quality_vlnplot <- function(data_S_list, file = "./figure/quality_control.pdf", width=8,height=4,log = F){
  pdf(file, width = width, height=height)
  for(i in names(data_S_list)){
    print(VlnPlot(data_S_list[[i]], idents = data_S_list[[i]]@meta.data$orig.ident, log = log,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05))
  }
  dev.off()
}

get_output_name <- function(main_name, prefix = "", type = ""){
  if (prefix == ""){
    if (type == "") return(main_name)
    if (type == "data") return(paste0("./data/", main_name))
    if (type == "figure") return(paste0("./figure/", main_name))
  }
  if (prefix != ""){
    if (type == "") return(paste0(prefix, "_", main_name))
    if (type == "data") return(paste0("./data/", prefix, "_", main_name))
    if (type == "figure") return(paste0("./figure/", prefix, "_", main_name))
  }
}

get_metric_summary <- function(data_S_list, data_S_updated_list = NULL, ngene_lth = 500, ntranscript_lth = 1000){
  metric_list <- list()
  metric_list[['Cell number']] <- unlist(lapply(data_S_list, function(x){
    ncol(x)
  }))
  metric_list[[paste('% Cells with >', ntranscript_lth, 'Transcripts')]] <- unlist(lapply(data_S_list, function(x){
    percent(mean(Matrix::colSums(x) > ntranscript_lth), accuracy = 0.1)
  }))
  metric_list[[paste('% Cells with >', ngene_lth, 'genes')]] <- unlist(lapply(data_S_list, function(x){
    percent(mean(Matrix::colSums(x@assays$RNA@counts > 0) > ngene_lth), accuracy = 0.1)
  }))
  metric_list[['Median # Transcripts']] <- unlist(lapply(data_S_list, function(x){
    median(Matrix::colSums(x))
  }))
  metric_list[['Median # Genes']] <- unlist(lapply(data_S_list, function(x){
    median(Matrix::colSums(x@assays$RNA@counts > 0))
  }))
  metric_list[['Median mitochondrial ratio']] <- unlist(lapply(data_S_list, function(x){
    percent(median(na.omit(x$percent.mt/100.)), accuracy = 0.1)
  }))
  if (!is.null(data_S_updated_list)){
    metric_list[['Cell number after filtering']] <- unlist(lapply(data_S_updated_list, function(x){
      ncol(x)
    }))
  }
  metric_df <- t(as.data.frame.list(metric_list))
  rownames(metric_df) <- names(metric_list)
  metric_df <- metric_df[c(1, 7, 2:6), ]
  metric_df
}


FilterCells <- function(data_S_list, ngene_lth = NULL, ngene_hth = NULL, mt_hth, nRNA_lth = NULL, nRNA_hth = NULL, ngene_lqth = NULL, ngene_hqth = NULL){
  if (is.null(ngene_lth)){
    for(i in names(data_S_list)){
      data_S_list[[i]]$ngene_lth <- quantile(data_S_list[[i]]$nFeature_RNA,ngene_lqth)
      data_S_list[[i]]$ngene_hth <- quantile(data_S_list[[i]]$nFeature_RNA,ngene_hqth)
      data_S_list[[i]]$mt_hth <- mt_hth
      data_S_list[[i]] <- subset(
        x = data_S_list[[i]], 
        subset = nFeature_RNA > ngene_lth & nFeature_RNA < ngene_hth & percent.mt <= mt_hth
      )
      data_S_list[[i]] <- RenameCells(data_S_list[[i]], add.cell.id = i)
    }
  }
  if (is.null(ngene_lqth)){
    for(i in names(data_S_list)){
      data_S_list[[i]]$ngene_lth <- ngene_lth
      data_S_list[[i]]$ngene_hth <- ngene_hth
      data_S_list[[i]]$mt_hth <- mt_hth
      #data_S_list[[i]]$nRNA_lth <- nRNA_lth
      #data_S_list[[i]]$nRNA_hth <- nRNA_hth
      data_S_list[[i]] <- subset(
        x = data_S_list[[i]], 
        #subset = nFeature_RNA > ngene_lth & nFeature_RNA <= ngene_hth & percent.mt < mt_hth & nCount_RNA > nRNA_lth & nCount_RNA < nRNA_hth
        subset = nFeature_RNA > ngene_lth & nFeature_RNA <= ngene_hth & percent.mt <= mt_hth
      )
      data_S_list[[i]] <- RenameCells(data_S_list[[i]], add.cell.id = i)
    }
  }
  data_S_list
}

MergeData <- function(data_S_list){
  sample_names <- names(data_S_list)
  n_sample <- length(sample_names)
  data_S <- merge(
    x = data_S_list[[1]], data_S_list[sample_names[2:n_sample]]
  )
  data_S
}

FindDEGs <- function(object, group.by = "all", compare.by = "orig.ident", assay = "RNA", features = NULL, min.cell = 25, 
                     min.pct = 0.01, logfc.threshold = log(2), pseudocount.use = 0.01, ...){
  #object <- sub_data_S
  #compare.by  <- "combined_da_tmp"
  DefaultAssay(object) <- assay
  DEGs_list <- list()
  Idents(object) <- compare.by
  #print(Idents(object))
  if (group.by == "all"){
    object$tmp_all_for_DEG <- "all"
    group.by <- "tmp_all_for_DEG"
  }
  object_all <- object
  for (group_i in unique(na.omit(object_all[[group.by]][,1]))){
    #print(group_i)
    object <- object_all[, which(object_all[[group.by]] == group_i)]
    print(paste("Processing", group_i, "..."))
    if(min(table(object[[compare.by]])) < min.cell | length(table(object[[compare.by]])) == 1){
      print(paste("Cell number for some condition is below", min.cell))
      next
    }
    cluster_markers <- FindAllMarkers(
      object, assay = assay, features = features, only.pos = T, min.pct = min.pct, logfc.threshold = logfc.threshold, 
      pseudocount.use = pseudocount.use, ...
    )
    cluster_markers$pct.diff <- cluster_markers$pct.1 - cluster_markers$pct.2
    #cluster_markers <- cbind(cluster_markers, AverageExpression(object, features = rownames(cluster_markers)), assays = assay)
    
    #for (compare_i in unique(object[[compare.by]][,1])){
    #sub_dat <- slot(object[, which(object[[compare.by]] == compare_i)]@assays[[assay]], name = "data")[rownames(cluster_markers),]
    #cluster_markers[, paste0("avg_", compare_i)] <- apply(sub_dat, 1, mean)
    #}
    DEGs_list[[group_i]] <- cluster_markers
  }
  DEGs_list
}

ProcessSingleSample <- function(data_S_list){
  for(i in names(data_S_list)){
    data_S_list[[i]] <- NormalizeData(data_S_list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
    data_S_list[[i]] <- FindVariableFeatures(data_S_list[[i]], selection.method = "vst", nfeatures = 2000)
    data_S_list[[i]] <- ScaleData(data_S_list[[i]])
    data_S_list[[i]] <- RunPCA(data_S_list[[i]], features = VariableFeatures(object = data_S_list[[i]]))
  }
  data_S_list
}

DoubletDetection <- function(data_S_list){
  for(i in names(data_S_list)){
    sweep.res.list <- paramSweep_v3(data_S_list[[i]], PCs = 1:30, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    pred_doub_percent <- nrow(data_S_list[[i]]@meta.data)*0.001*0.008
    nExp_poi <- round(pred_doub_percent*nrow(data_S_list[[i]]@meta.data))
    data_S_list[[i]] <- doubletFinder_v3(data_S_list[[i]], PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  }
  data_S_list
}

scDoubletDetection <- function(data_S_list){
  for(i in names(data_S_list)){
    #convert to sce
    data_S_list[[i]][["RNA"]] <- as(object = data_S_list[[i]][["RNA"]], Class = "Assay")
    data_sce <- as.SingleCellExperiment(data_S_list[[i]])
    data_sce <- scDblFinder(data_sce)
    data_S_list[[i]]$scDblFinder.class <- data_sce$scDblFinder.class
    data_S_list[[i]]$scDblFinder.score <- data_sce$scDblFinder.score
  }
  data_S_list
}


ClusterMetricTable <- function(seu_obj, cluster_col, metric_cols){
  for(col in metric_cols){
    metric_table <- table(seu_obj@meta.data[[cluster_col]], seu_obj@meta.data[[col]])
    write.table(metric_table, file = paste0("cluster_", col, "_metrics.txt"))
    metric_df <- as.data.frame.matrix(metric_table) 
    metric_df_long <- metric_df
    #fixing error where numerical rownames are out of order
    if(suppressWarnings(is.na(as.integer(rownames(metric_df)[1])))){
      metric_df_long$cluster <- rownames(metric_df)
    }
    else{
      metric_df_long$cluster <- as.integer(rownames(metric_df))
    }
    metric_df_long <- melt(metric_df_long, id.vars = "cluster")
    colnames(metric_df_long)[colnames(metric_df_long) == "value"] ="cell_counts"
    colnames(metric_df_long)[colnames(metric_df_long) == "variable"] ="sample_id"
    return(metric_df_long)
  }
}



StackedVlnMarkerPlot3 <- function(features_df, seu_obj, plot_name, col_name,feat_level){
  #get the total list of genes and make sure features selected
  #are in the list
  #features_df$grouping <- as.factor(features_df$grouping)
  all.genes <- row.names(seu_obj)
  features <- features_df$Feat[(features_df$Feat %in% all.genes)]
  
  #get the metadata of object for cluster information
  cell_meta <- seu_obj@meta.data
  cluster_info <- cell_meta[c(col_name)]
  counts <- as.data.frame(as.matrix(seu_obj@assays$RNA@data[features,]))
  counts <- t(counts)
  count_cluster_df <- merge(counts, cluster_info,  by = 'row.names', all = TRUE)
  #return(count_cluster_df)
  count_cluster_df$Cell <- count_cluster_df$Row.names
  count_cluster_df$Idents <- count_cluster_df[[col_name]]
  
  count_cluster_df <- reshape2::melt(count_cluster_df, id.vars = c("Cell","Idents"), measure.vars = features,
                                     variable.name = "Feat", value.name = "Expr")
  
  count_cluster_df <- left_join(count_cluster_df,features_df, by="Feat")
  #return(count_cluster_df)
  
  #need to refactor if you have subsetted
  #count_cluster_df$Idents <- factor(count_cluster_df$Idents)
  
  #avg <- sapply(X = split(x = count_cluster_df, f = count_cluster_df$Idents),
  #              FUN = function(df) { return(tapply(X = df$Expr, INDEX = df$Feat, FUN = mean)) })
  
  
  #L2Norm <- function(mat, MARGIN){
  #  normalized <- sweep(x = mat, MARGIN = MARGIN,
  #                      STATS = apply(X = mat, MARGIN = MARGIN,
  #                                    FUN = function(x){ sqrt(x = sum(x ^ 2)) }), FUN = "/")
  #  normalized[!is.finite(x = normalized)] <- 0
  #  return(normalized)
  #}
  
  # Performs hierarchical clustering
  #idents.order <- hclust(d = dist(t(L2Norm(mat = avg, MARGIN = 2))))$order
  #avg <- avg[,idents.order]
  #avg <- L2Norm(mat = avg, MARGIN = 1)
  #mat <- hclust(d = dist(avg))$merge
  
  # Order feature clusters by position of their "rank-1 idents"
  #position <- apply(X = avg, MARGIN = 1, FUN = which.max)
  #orderings <- list()
  #for (i in 1:nrow(mat)) {
  #  x <- if (mat[i,1] < 0) -mat[i,1] else orderings[[mat[i,1]]]
  #  y <- if (mat[i,2] < 0) -mat[i,2] else orderings[[mat[i,2]]]
  #  x.pos <- min(x = position[x])
  #  y.pos <- min(x = position[y])
  #  orderings[[i]] <- if (x.pos < y.pos) { c(x, y) } else { c(y, x) }
  #}
  #features.order <- orderings[[length(orderings)]]
  
  # Update Feature and Identity factor orders
  #count_cluster_df$Idents <- factor(count_cluster_df$Idents, levels = ident_level)
  count_cluster_df$Feat <- factor(count_cluster_df$Feat, levels = feat_level)
  
  #feature_labeller <- as_labeller(variable,value)
  
  # Plot stacked violin plot with reordered identity classes and features - horizontal
  f <- ggplot(count_cluster_df, aes(Expr, Idents, fill = Feat)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    scale_x_continuous(expand = c(0, 0), labels = function(x)
      c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    facet_grid(cols = vars(Feat), scales = "free")  +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "none", panel.spacing = unit(0, "lines"),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = NA, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
    ggtitle(plot_name) + xlab("Expression Level") + ylab("Identity")
  
  
  
  return(f)
}

StackedVlnMarkerPlot4 <- function(features_df, seu_obj, plot_name, col_name,sample_cols,cluster_level){
  #get the total list of genes and make sure features selected
  #are in the list
  features_df$grouping <- as.factor(features_df$grouping)
  all.genes <- row.names(seu_obj)
  features <- features_df$Feat[(features_df$Feat %in% all.genes)]
  
  #get the metadata of object for cluster information
  cell_meta <- seu_obj@meta.data
  cluster_info <- cell_meta[c(col_name)]
  counts <- as.data.frame(as.matrix(seu_obj@assays$RNA@data[features,]))
  counts <- t(counts)
  count_cluster_df <- merge(counts, cluster_info,  by = 'row.names', all = TRUE)
  #return(count_cluster_df)
  count_cluster_df$Cell <- count_cluster_df$Row.names
  count_cluster_df$Idents <- count_cluster_df[[col_name]]
  
  count_cluster_df <- reshape2::melt(count_cluster_df, id.vars = c("Cell","Idents"), measure.vars = features,
                                     variable.name = "Feat", value.name = "Expr")
  
  count_cluster_df <- left_join(count_cluster_df,features_df, by="Feat")
  #return(count_cluster_df)
  
  #need to refactor if you have subsetted
  count_cluster_df$Idents <- factor(count_cluster_df$Idents)
  
  avg <- sapply(X = split(x = count_cluster_df, f = count_cluster_df$Idents),
                FUN = function(df) { return(tapply(X = df$Expr, INDEX = df$Feat, FUN = mean)) })
  
  
  L2Norm <- function(mat, MARGIN){
    normalized <- sweep(x = mat, MARGIN = MARGIN,
                        STATS = apply(X = mat, MARGIN = MARGIN,
                                      FUN = function(x){ sqrt(x = sum(x ^ 2)) }), FUN = "/")
    normalized[!is.finite(x = normalized)] <- 0
    return(normalized)
  }
  
  # Performs hierarchical clustering
  idents.order <- hclust(d = dist(t(L2Norm(mat = avg, MARGIN = 2))))$order
  avg <- avg[,idents.order]
  avg <- L2Norm(mat = avg, MARGIN = 1)
  mat <- hclust(d = dist(avg))$merge
  
  # Order feature clusters by position of their "rank-1 idents"
  position <- apply(X = avg, MARGIN = 1, FUN = which.max)
  orderings <- list()
  for (i in 1:nrow(mat)) {
    x <- if (mat[i,1] < 0) -mat[i,1] else orderings[[mat[i,1]]]
    y <- if (mat[i,2] < 0) -mat[i,2] else orderings[[mat[i,2]]]
    x.pos <- min(x = position[x])
    y.pos <- min(x = position[y])
    orderings[[i]] <- if (x.pos < y.pos) { c(x, y) } else { c(y, x) }
  }
  #features.order <- orderings[[length(orderings)]]
  
  # Update Feature and Identity factor orders
  count_cluster_df$Idents <- factor(count_cluster_df$Idents, levels = levels(count_cluster_df$Idents)[idents.order])
  #count_cluster_df$Feat <- factor(count_cluster_df$Feat, levels = levels(count_cluster_df$Feat)[features.order])
  
  #feature_labeller <- as_labeller(variable,value)
  
  # Plot stacked violin plot with reordered identity classes and features - horizontal
  f <- ggplot(count_cluster_df, aes(Expr, factor(Idents), fill = Feat)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    scale_x_continuous(expand = c(0, 0), labels = function(x)
      c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    facet_grid(cols = vars(grouping, Feat), scales = "free")  +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "none", panel.spacing = unit(0, "lines"),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = NA, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
    ggtitle(plot_name) + xlab("Expression Level") + ylab("Identity")
  
  metric_table <- ClusterMetricTable(seu_obj, col_name, c("condition", "label"))
  #return(metric_table)
  metric_table$cluster <- factor(metric_table$cluster,levels=cluster_level)
  #metric_table$cluster <- factor(metric_table$cluster, levels = levels(metric_table$cluster)[idents.order])
  metric_table$condition <- factor(sapply(as.character(metric_table$sample_id),function(name) strsplit(name,split = "_")[[1]][1]))
  
  #add total counts
  totals <- metric_table %>%
    group_by(cluster) %>%
    summarize(total = sum(cell_counts))
  
  t <- ggplot(metric_table) +
    aes(x = cluster, y = cell_counts, fill = sample_id, label = sample_id,col=condition) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=sample_cols)+ 
    labs(x = "cluster", title = "Cell Counts (Percentage)") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    ) +
    geom_text(data=totals, aes(x=cluster, label=total, y=total, fill=NULL), position = position_fill(vjust = -0.05)) +
    coord_flip()
  
  return(metric_table)
}


# Differntially Expressed Genes using Libra
RunLibra <- function(SeuratObj, celltype_col_name, id_col_name, label_col_name, 
                     num_replicates){
  meta = subset(SeuratObj@meta.data, select=c(({{celltype_col_name}}), 
                                              ({{id_col_name}}), ({{label_col_name}})))
  colnames(meta)[which(colnames(meta)== ({{celltype_col_name}}))] = 'cell_type'
  colnames(meta)[colnames(meta)==({{id_col_name}})] = 'replicate'
  X = SeuratObj@assays$RNA@counts
  #missing_cell_type = !(is.na(meta$cell_type))
  #X = X[,missing_cell_type]
  #meta = meta[missing_cell_type,]
  DE = Libra::run_de(X,meta=meta,min_cells=50,min_reps={{num_replicates}},min_features=0)
  return(DE)
}

RunEdgeRPseudobulkBatchEffects <- function(SeuratObj, celltype_col_name, id_col_name, 
                                           label_col_name, num_replicates,
                                           batch_effect=TRUE){
  #EdgeR LRT DE w/ batch effects, but Libra Formatting
  meta = subset(SeuratObj@meta.data, select=c(({{celltype_col_name}}), 
                                              ({{id_col_name}}), ({{label_col_name}})))
  colnames(meta)[which(colnames(meta)== ({{celltype_col_name}}))] = 'cell_type'
  colnames(meta)[colnames(meta)==({{id_col_name}})] = 'replicate'
  X = SeuratObj@assays$RNA@counts
  
  pseudobulks = Libra::to_pseudobulk(
    input = X,
    meta = meta,
    replicate_col = 'replicate',
    cell_type_col = 'cell_type',
    label_col = label_col_name,
    min_cells = 50,
    min_reps = {{num_replicates}},
    min_features = 0,
  )
  
  de_type='pseudobulk'
  de_method='edgeR'
  de_family='LRT'
  
  results = map(pseudobulks, function(x) {
    # create targets matrix
    
    
    # create design
    
    if (batch_effect){
      targets = data.frame(group_sample = colnames(x)) %>%
        mutate(group = gsub(".*\\:", "", group_sample)) %>% 
        mutate(batch = gsub("\\D","", group_sample))
      ## optionally, carry over factor levels from entire dataset
      if (is.factor(meta$label)) {
        targets$group %<>% factor(levels = levels(meta$label))
      }
      if (n_distinct(targets$group) > 2)
        return(NULL)
      design = model.matrix(~batch+group, data = targets)
    }
    else{
      targets = data.frame(group_sample = colnames(x)) %>%
        mutate(group = gsub(".*\\:", "", group_sample))
      ## optionally, carry over factor levels from entire dataset
      if (is.factor(meta$label)) {
        targets$group %<>% factor(levels = levels(meta$label))
      }
      if (n_distinct(targets$group) > 2)
        return(NULL)
      design = model.matrix(~group, data = targets)
    }
    DE = tryCatch({
      y = DGEList(counts = x, group = targets$group) %>%
        calcNormFactors(method = 'TMM') %>%
        estimateDisp(design)
      test = {
        fit = glmFit(y, design = design)
        test = glmLRT(fit)
      }
      res = topTags(test, n = Inf) %>%
        as.data.frame() %>%
        rownames_to_column('gene') %>%
        # flag metrics in results
        mutate(de_family = 'pseudobulk',
               de_method = de_method,
               de_type = de_type)
    }, error = function(e) {
      message(e)
      data.frame()
    })
    
    
  })
  results %<>% bind_rows(.id = 'cell_type')
  DE <- results
  suppressWarnings(
    colnames(DE) %<>%
      fct_recode('p_val' = 'p.value',  ## DESeq2
                 'p_val' = 'pvalue',  ## DESeq2
                 'p_val' = 'p.value',  ## t/wilcox
                 'p_val' = 'P.Value',  ## limma
                 'p_val' = 'PValue'  , ## edgeR
                 'p_val_adj' = 'padj', ## DESeq2/t/wilcox
                 'p_val_adj' = 'adj.P.Val',      ## limma
                 'p_val_adj' = 'FDR',            ## edgeER
                 'avg_logFC' = 'log2FoldChange', ## DESEeq2
                 'avg_logFC' = 'logFC', ## limma/edgeR
                 'avg_logFC' = 'avg_log2FC' # Seurat V4
      )
  ) %>%
    as.character()
  
  DE %<>%
    # calculate adjusted p values
    group_by(cell_type) %>%
    mutate(p_val_adj = p.adjust(p_val, method = 'BH')) %>%
    # make sure gene is a character not a factor
    mutate(gene = as.character(gene)) %>%
    # invert logFC to match Seurat level coding
    mutate(avg_logFC = avg_logFC * -1) %>%
    dplyr::select(cell_type,
                  gene,
                  avg_logFC,
                  p_val,
                  p_val_adj,
                  de_family,
                  de_method,
                  de_type
    ) %>%
    ungroup() %>%
    arrange(cell_type, gene)
  return(DE)
}
