`%notin%` <- Negate(`%in%`)

tabla <- function(...) {
  call <- as.list(match.call())[-1] # first position is the function_name
  
  custom_args <- list(useNA = "no") # could extend this list for more customization
  
  overlap_args <-  names(call) %in% names(custom_args) # handle overlapping args
  if (!any(overlap_args)) call <- c(call, custom_args)
  
  do.call(table, call) # exectue table() with the custom settings
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

get_quality_vlnplot_v0 <- function(data_S_list, file = "./figure/quality_control.pdf"){
  pdf(file, width = 8, height=4)
  for(i in names(data_S_list)){
    print(VlnPlot(data_S_list[[i]], idents = data_S_list[[i]]@meta.data$orig.ident, log =T,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.25))
  }
  dev.off()
}

get_quality_vlnplot <- function(data_S_list, file = "./figure/quality_control.pdf", log = F){
  pdf(file, width = 8, height=4)
  for(i in names(data_S_list)){
    print(VlnPlot(data_S_list[[i]], idents = data_S_list[[i]]@meta.data$orig.ident, log = log,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05))
  }
  dev.off()
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

FilterCells <- function(data_S_list, ngene_lth = NULL, ngene_hth = NULL, mt_hth, ngene_lqth = NULL, ngene_hqth = NULL){
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
      data_S_list[[i]] <- subset(
        x = data_S_list[[i]], 
        subset = nFeature_RNA > ngene_lth & nFeature_RNA <= ngene_hth & percent.mt <= mt_hth
      )
      data_S_list[[i]] <- RenameCells(data_S_list[[i]], add.cell.id = i)
    }
  }
  
  data_S_list
}

Load10xData <- function(data_dir, sample_names, sub_rna_dir = "filtered_feature_bc_matrix", with.protein = F, with.tcr = F, sub_tcr_dir = "TCR", human_only = F){
  data_S_list <- list()
  data_dir0 <- data_dir
  if (!with.protein){
    for(i in sample_names){
      data_dir <- paste(data_dir0, i, sub_rna_dir, sep = "/")
      print(data_dir)
      data_S_list[[i]] <- CreateSeuratObject(Read10X(data.dir = data_dir), project = i)
    }
  }
  if (with.protein){
    # For output from CellRanger >= 3.0 with multiple data types
    #list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
    data_dir0 <- data_dir
    for(i in sample_names){
      data_dir <- paste(data_dir0, i, sub_rna_dir, sep = "/")
      print(data_dir)
      data <- Read10X(data.dir = data_dir)
      data_S_list[[i]] <- CreateSeuratObject(counts = data$`Gene Expression`, project = i)
      data_S_list[[i]][["ADT"]] <- CreateAssayObject(counts = data$`Antibody Capture`)
    }
  }
  if (with.tcr){
    for(i in sample_names){
      data_dir <- paste(data_dir0, i, sub_tcr_dir, sep = "/")
      print(data_dir)
      data_S_list[[i]] <- add_clonotype(data_S_list[[i]], data_dir)
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

MergeData <- function(data_S_list){
  sample_names <- names(data_S_list)
  n_sample <- length(sample_names)
  data_S <- merge(
    x = data_S_list[[1]], data_S_list[sample_names[2:n_sample]]
  )
  data_S
}

add_clonotype_v0 <- function(seurat_obj, tcr_folder,  prefix = ""){
  tcr_orig <- read.csv(paste(tcr_folder,"filtered_contig_annotations.csv", sep="/"))
  
  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  tcr_orig$barcode <- gsub("-1", "", tcr_orig$barcode)
  tcr <- tcr_orig[!duplicated(tcr_orig$barcode), ]
  
  # Only keep the barcode and clonotype columns. 
  # We'll get additional clonotype info from the clonotype table.
  tcr <- tcr[,c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_folder,"clonotypes.csv", sep="/"))
  
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  
  # Reorder so barcodes are first column and set them as rownames.
  extract_TRB <- function(barcode, trb_dat){
    cdr3 <- trb_dat[trb_dat$barcode == barcode,]$cdr3
    if (all(cdr3 == "None") | is.null(cdr3)) return("")
    else paste(sort(unique(cdr3[cdr3 != "None"])), collapse = "_")
  }
  rownames(tcr) <- paste0(prefix,tcr$barcode)
  tcr <- tcr[, c("clonotype_id", "cdr3s_aa")]
  tcr_orig_trb <- tcr_orig[tcr_orig$chain == "TRB",]
  tcr$TRB_aa <- as.factor(sapply(rownames(tcr), FUN = function(x){extract_TRB(x, tcr_orig_trb)}))

  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}

add_clonotype <- function(seurat_obj, tcr_folder,  prefix = ""){
  #tcr_folder <- data_dir
  tcr_orig <- read.csv(paste(tcr_folder,"filtered_contig_annotations.csv", sep="/"), stringsAsFactors = F)
  names(tcr_orig)[names(tcr_orig) == "raw_clonotype_id"] <- "clonotype_id"
  clono <- read.csv(paste(tcr_folder,"clonotypes.csv", sep="/"), stringsAsFactors = F)
  tcr_orig <- merge(tcr_orig, clono)
  tcr_orig_trb <- tcr_orig[tcr_orig$chain == "TRB",]
  tcr_orig_tra <- tcr_orig[tcr_orig$chain == "TRA",]
  valid_barcodes <- names(which(table(tcr_orig_trb$barcode) == 1))
  print(paste(length(valid_barcodes), "valid barcodes with unique TRB are found."))
  #length(names(which(table(tcr_orig_tra$barcode) == 1)))
  #length(intersect(names(which(table(tcr_orig_trb$barcode) == 1)),
  #                 names(which(table(tcr_orig_tra$barcode) == 1))))
  selected_tcr_orig_trb <- tcr_orig_trb[which(tcr_orig_trb$barcode %in% valid_barcodes),]
  selected_tcr_orig_trb_to_incorp <- data.frame(clonotype_id = selected_tcr_orig_trb$clonotype_id,
                                                contig_id = selected_tcr_orig_trb$contig_id,
                                                v_gene_TRB = selected_tcr_orig_trb$v_gene,
                                                d_gene_TRB = selected_tcr_orig_trb$d_gene,
                                                j_gene_TRB = selected_tcr_orig_trb$j_gene,
                                                c_gene_TRB = selected_tcr_orig_trb$c_gene,
                                                cdr3_TRB = selected_tcr_orig_trb$cdr3,
                                                cdr3_nt_TRB = selected_tcr_orig_trb$cdr3_nt,
                                                cdr3s_aa = selected_tcr_orig_trb$cdr3s_aa,
                                                cdr3s_nt = selected_tcr_orig_trb$cdr3s_nt)
  rownames(selected_tcr_orig_trb_to_incorp) <- selected_tcr_orig_trb$barcode
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=selected_tcr_orig_trb_to_incorp)
  
  selected_tcr_orig_tra_to_incorp <- data.frame(clonotype_id = selected_tcr_orig_trb$clonotype_id,
                                                v_gene_TRA = NA,
                                                d_gene_TRA = NA,
                                                j_gene_TRA = NA,
                                                c_gene_TRA = NA,
                                                cdr3_TRA = NA,
                                                cdr3_nt_TRA = NA)
  rownames(selected_tcr_orig_tra_to_incorp) <- rownames(selected_tcr_orig_trb_to_incorp)
  for (i in rownames(selected_tcr_orig_tra_to_incorp)){
    selected_tcr_orig_tra_to_incorp[i, "v_gene_TRA"] = paste(tcr_orig_tra[tcr_orig_tra$barcode == i, "v_gene"], collapse = ";")
    selected_tcr_orig_tra_to_incorp[i, "d_gene_TRA"] = paste(tcr_orig_tra[tcr_orig_tra$barcode == i, "d_gene"], collapse = ";")
    selected_tcr_orig_tra_to_incorp[i, "j_gene_TRA"] = paste(tcr_orig_tra[tcr_orig_tra$barcode == i, "j_gene"], collapse = ";")
    selected_tcr_orig_tra_to_incorp[i, "c_gene_TRA"] = paste(tcr_orig_tra[tcr_orig_tra$barcode == i, "c_gene"], collapse = ";")
    selected_tcr_orig_tra_to_incorp[i, "cdr3_TRA"] = paste(tcr_orig_tra[tcr_orig_tra$barcode == i, "cdr3"], collapse = ";")
    selected_tcr_orig_tra_to_incorp[i, "cdr3_nt_TRA"] = paste(tcr_orig_tra[tcr_orig_tra$barcode == i, "cdr3_nt"], collapse = ";")
  }
  selected_tcr_orig_tra_to_incorp[] <- lapply(selected_tcr_orig_tra_to_incorp, factor)
  clono_seurat <- AddMetaData(object=clono_seurat, metadata=selected_tcr_orig_tra_to_incorp)
  #clono_seurat$v_gene_TRB[1:10]
  #clono_seurat$v_gene_TRA[1:10]
  return(clono_seurat)
}

############
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

SaveDEGs <- function(DEGs_list, path = "./data/differential_analysis/", prefix = "", suffix = "_DEGs", format = "xlsx", ...){
  if (!dir.exists(path)) dir.create(path)
  if (format == "xlsx"){
    require(openxlsx)
    for(i in names(DEGs_list)){
      DEGs_list[[i]] <- cbind(gene_name = rownames(DEGs_list[[i]]), DEGs_list[[i]])
      write.xlsx(DEGs_list[[i]], paste0(path, prefix, i, suffix, ".", format))
    }
  }
}

CreateIdentByGene <- function(object, gene_id, pos_lab = "+", neg_lab = "-", assay = "RNA", slot = "data"){
  object[[gene_id]] <- ifelse(slot(object@assays[[assay]], slot)[gene_id, ] > 0, paste0(gene_id, pos_lab), paste0(gene_id, neg_lab))
  object
}


suppressPackageStartupMessages({
  library(rlang)
})

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by)]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    group.use <- groups.use[, c(i, additional.group.by), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(foo=rep(x = levels(x = group.use[[i]]), times = lines.width))
      placeholder.groups[additional.group.by] = NA
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.levels <- levels(x = group.use[[i]])
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    #group.use = group.use[order(group.use[[i]]), , drop=F]
    group.use <- group.use[with(group.use, eval(parse(text=paste('order(', paste(c(i, additional.group.by), collapse=', '), ')', sep='')))), , drop=F]
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))){
        if (colname == group.by){
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))), "#FFFFFF")
        } else {
          cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = grid::gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        #temp <- as.data.frame(cols[[colname]][levels(group.use[[colname]])])
        #colnames(temp) <- 'color'
        #temp$x <- temp$y <- 1
        #temp[['name']] <- as.factor(rownames(temp))
        
        #temp <- ggplot(temp, aes(x=x, y=y, fill=name)) + geom_point(shape=21, size=5) + labs(fill=colname) + theme(legend.position = "bottom")
        #legend <- get_legend(temp)
        #multiplot(plot, legend, heights=3,1)
        
        if ((colname == group.by) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}




suppressPackageStartupMessages({
  library(rlang)
})

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          #x.divs <- pbuild$layout$panel_params[[1]]$x.major
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}



FeaturePlotNew <- function(object, features, gradient = NULL, combine = T, ...){
  if(is.null(gradient)){
    gradient <- rev(brewer_pal(palette = "RdYlBu")(5))
  }
  if(combine){
    plot_grid(plotlist = lapply(
      FeaturePlot(object = object, features = features, combine = F, ...),
      function(x) x + scale_color_gradientn(colors = gradient)
    ))
  }else{
    lapply(
      FeaturePlot(object = object, features = features, combine = F, ...),
      function(x) x + scale_color_gradientn(colors = gradient)
    )
  }
}
