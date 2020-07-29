#' Some Useful Utilization function for Seurat object
#' Theses function is needed to packaged into a pack
#' @author  wlt
#' @description Some Useful Utilization function for Seurat object
#' @param seuratv3 object
#' @importFrom Seurat
#' @importFrom tidyverse
#' @importFrom ComplexHeamp
#' @importFrom RColorBrewer
###----myScHeatmap-----
#' @description drawheatmap, rows are interested gene, columns are interested cells;
#' @description values can be choose raw log-normlized expression or scale expresion
myScHeatmap <- function (object=NULL, 
                         features = NULL, 
                         cells = NULL,
                         slot = "data", 
                         assay = "RNA",
                         base.sample=NULL,
                         idents.subset = NULL,
                         sample.subset = NULL,
                         idents.order=tmp.cell.order,
                         scale.min=-4,
                         scale.max=4, 
                         show_rownames = T,
                         show_columnnames = F,
                         show_row_dend = T,
                         cluster_rows=F,
                         cluster_columns=F,
                         color.plot=NULL,
                         raster = TRUE,
                         raster.device="CairoPNG") {
  ### check input
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were assigned to zero as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  ### prepare data to plot
  data <- as.data.frame(GetAssayData(object = object, slot = slot)[features, cells, drop = FALSE])
  if(exists("bad.features")){
    tmp.df <- data.frame(matrix(data = rep(0,length(bad.features)*ncol(data)),
                                nrow = length(bad.features)),
                         row.names = bad.features,
                         stringsAsFactors = F)
    colnames(tmp.df) <- colnames(data)
    data <- rbind(tmp.df,data)
    features <- rownames(data)
  }
  
  if (any(!c("cell.type","sample") %in% colnames(object[[]]))) {
    stop("Please add cell.type and sample to seurat metadata")
  }
  data.meta <- object[[]] %>%
    select(cell.type,sample) 
  data.meta$cell.type <- factor(data.meta$cell.type,
                                levels = idents.order)
  data.meta <- data.meta %>%
    rownames_to_column("cell.id") %>%
    arrange(cell.type,sample) %>%
    column_to_rownames("cell.id")
  
  if(is.null(color.plot)){
    cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
    warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
    warmcold <- colorRampPalette(c(rev(cold(21)), warm(20)))
    color.plot <- warmcold(100)
  }
  
  set.seed(42)
  ha <- HeatmapAnnotation(df = data.meta)
  tmp.color.maplist <- get_color_mapping_list(ha)
  
  
  if(is.null(idents.subset)){
    idents.subset <- unique(as.character(data.meta$cell.type))
  }
  if(is.null(sample.subset)){
    sample.subset <- unique(as.character(data.meta$sample))
  }
  
  if(is.null(base.sample)){
    tmp.meta <- data.meta %>%
      rownames_to_column("cell.id") %>%
      filter(cell.type %in% idents.subset) %>%
      filter(sample %in% sample.subset) %>%
      arrange(cell.type,sample) %>%
      column_to_rownames("cell.id")
    
    tmp.plot <- data[features,rownames(tmp.meta)]
    if(slot=="scale.data"){
      tmp.plot[tmp.plot > scale.max] = scale.max
      tmp.plot[tmp.plot < scale.min] = scale.min
    }
    
    
    ha <- HeatmapAnnotation(df = tmp.meta,
                            col = list(cell.type=tmp.color.maplist[[1]]@colors[idents.subset],
                                       sample=tmp.color.maplist[[2]]@colors[sample.subset]))
    col_funs <- circlize::colorRamp2(seq(min(tmp.plot),max(tmp.plot),length.out = length(color.plot)),
                                     colors = color.plot)
    
    ht <- Heatmap(matrix = tmp.plot,
                  name = "Exp", 
                  show_column_names = show_columnnames,
                  show_row_names = show_rownames,
                  cluster_columns = cluster_columns,
                  cluster_rows = cluster_rows,
                  show_row_dend = show_row_dend,
                  col = col_funs,
                  use_raster = raster,
                  heatmap_legend_param = list(at = pretty(c(min(tmp.plot)+0.5, 0, max(tmp.plot)-0.5),n=5)),
                  row_names_gp = gpar(fontsize = 8), 
                  top_annotation = ha,
                  column_split = tmp.meta$cell.type,
                  column_title = NULL,
                  raster_device = raster.device)
  }else{
    
    if (slot!="data") {
      stop("Please using data slot")
    }
    
    tmp.res.list <- sapply(
      idents.subset,
      FUN = function(ii){
        tmp.meta <- data.meta %>%
          rownames_to_column("cell.id") %>%
          filter(cell.type %in% ii) %>%
          filter(sample %in% sample.subset) %>%
          arrange(cell.type,sample) %>%
          column_to_rownames("cell.id")
        
        tmp.base.index <- data.meta %>%
          rownames_to_column("cell.id") %>%
          filter(sample==base.sample) %>%
          filter(cell.type==ii) %>%
          column_to_rownames("cell.id")
        tmp.base.index <- rownames(tmp.base.index)
        
        tmp.base <- rowMeans(data[,tmp.base.index])
        tmp.base.sd <- ifelse(tmp.base==0,1,apply(data[,tmp.base.index],1,sd))
        
        tmp.plot <- data[,rownames(tmp.meta)] %>%
          mutate_all(., ~(.x-tmp.base)/tmp.base.sd)
        tmp.plot[is.na(tmp.plot)] <- 0
        rownames(tmp.plot) <- rownames(data)
        return(tmp.plot)
      },
      simplify = F,
      USE.NAMES = T
    )
    
    tmp.plot <-  bind_cols(tmp.res.list) %>%
      as.data.frame %>%
      `row.names<-`(., row.names(tmp.res.list[[1]]))
    
    ###test distribution
    ###summary(unlist(tmp.plot))
    
    tmp.meta <- data.meta %>%
      rownames_to_column("cell.id") %>%
      filter(cell.type %in% idents.subset) %>%
      filter(sample %in% sample.subset) %>%
      arrange(cell.type,sample) %>%
      column_to_rownames("cell.id")
    
    tmp.plot <- tmp.plot[features,rownames(tmp.meta)]
    tmp.plot[tmp.plot > scale.max] = scale.max
    tmp.plot[tmp.plot < scale.min] = scale.min
    
    ha <- HeatmapAnnotation(df = tmp.meta,
                            col = list(cell.type=tmp.color.maplist[[1]]@colors[idents.subset],
                                       sample=tmp.color.maplist[[2]]@colors[sample.subset]))
    col_funs <- circlize::colorRamp2(seq(scale.min,scale.max,
                                         length.out = length(color.plot)),
                                     colors = color.plot)
    ht <- Heatmap(matrix = tmp.plot,
                  name = "z-score", 
                  show_column_names = show_columnnames,
                  show_row_names = show_rownames,
                  cluster_columns = cluster_columns,
                  cluster_rows = cluster_rows,
                  show_row_dend = show_row_dend,
                  col = col_funs,
                  use_raster = raster,
                  heatmap_legend_param = list(at = pretty(c(min(tmp.plot)+0.5, 0, max(tmp.plot)-0.5),n=5)),
                  row_names_gp = gpar(fontsize = 8), 
                  top_annotation = ha,
                  column_split = tmp.meta$cell.type,
                  column_title = NULL,
                  raster_device = raster.device)
    
  }
  
  return(ht)
}


###----myPseudoBulkHeatmap----
#' @description draw heatmap,rows are intersted gene,cols are interested cells
#' @description values are row scale expression or based_sample scale expression
myPseudoBulkHeatmap <- function(object=pbmc, 
                                features = levels(p$data$Feature), 
                                cells = NULL,
                                slot = "data", 
                                assay = "RNA",
                                seed = 42,
                                pseudo.rep = 3,
                                idents.order=NULL,
                                sample.order=NULL,
                                idents.subset = NULL,
                                sample.subset = NULL,
                                base.sample = NULL,
                                useExp.mean = T,
                                scale.min=-2.5,
                                scale.max=2.5, 
                                show_rownames=T,
                                show_columnnames=T,
                                cluster_rows=F,
                                cluster_columns=F,
                                show_row_dend=T,
                                color.plot=mypal,
                                raster = TRUE,
                                raster.device="CairoPNG",
                                return.data=T) {
  
  ### check input
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were assigned to zero as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  ### prepare data to plot
  data <- as.data.frame(GetAssayData(object = object, slot = slot)[features, cells, drop = FALSE])
  if(exists("bad.features")){
    tmp.df <- data.frame(matrix(data = rep(0,length(bad.features)*ncol(data)),
                                nrow = length(bad.features)),
                         row.names = bad.features,
                         stringsAsFactors = F)
    colnames(tmp.df) <- colnames(data)
    data <- rbind(tmp.df,data)
    features <- rownames(data)
  }
  
  if (any(!c("cell.type","sample") %in% colnames(object[[]]))) {
    stop("Please add cell.type and sample to seurat metadata")
  }
  data.meta <- object[[]] %>%
    select(cell.type,sample) 
  
  data.meta$cell.type <- factor(data.meta$cell.type,
                                levels = idents.order)
  data.meta$sample <- factor(data.meta$sample,
                             levels = sample.order)
  
  data.meta <- data.meta %>%
    rownames_to_column("cell.id") %>%
    arrange(cell.type,sample) %>%
    column_to_rownames("cell.id")
  
  
  if(is.null(color.plot)){
    cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
    warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
    warmcold <- colorRampPalette(c(rev(cold(21)), warm(20)))
    color.plot <- warmcold(100)
  }
  
  
  
  if(is.null(idents.subset)){
    idents.subset <- levels(data.meta$cell.type)
  }
  if(is.null(sample.subset)){
    sample.subset <- levels(data.meta$sample)
  }
  
  ###PseudoRep is strategy random parrtation of the the sample
  set.seed(seed)
  tmp.meta <-  data.meta %>%
    rownames_to_column("cell.id") %>%
    filter(cell.type %in% idents.subset) %>%
    filter(sample %in% sample.subset) %>%
    group_by(cell.type,sample) %>%
    mutate(pseudo.replicate=sample(paste0("PseudoRep",1:pseudo.rep), 
                                   size = length(sample), 
                                   replace = TRUE)) %>%
    arrange(cell.type,sample,pseudo.replicate) %>%
    ungroup()
  
  # test partition reseult
  # tmp.meta %>%
  #   count(cell.type,sample,pseudo.replicate) %>%
  #   summarise(res=sum(n)) %>%
  #   unlist()
  
  tmp.plot <- data[,tmp.meta$cell.id] %>%
    rownames_to_column("gene.name") %>%
    gather(key = cell.id,
           value = gene.exp, 
           -gene.name,
           factor_key=F) %>%
    left_join(tmp.meta) %>% 
    group_by(cell.type,sample,pseudo.replicate,gene.name) %>%
    summarise(mean.gene.exp=ifelse(useExp.mean,mean(expm1(gene.exp)),mean(gene.exp))) %>%
    group_by(gene.name) %>%
    mutate(MetaCell=paste0("MetaCell",1:length(gene.name))) %>%
    ###relase grouping infor
    ungroup()
  
  if(pseudo.rep>1){
    tmp.meta <- tmp.plot %>%
      select(-mean.gene.exp,-gene.name) %>%
      unique() %>%
      column_to_rownames("MetaCell")
  }else{
    tmp.meta <- tmp.plot %>%
      select(-mean.gene.exp,-gene.name,-pseudo.replicate) %>%
      unique() %>%
      column_to_rownames("MetaCell")
  }
  
  if(is.null(base.sample)){
    tmp.plot <- tmp.plot %>%
      select(gene.name,mean.gene.exp,MetaCell) %>%
      spread(key = MetaCell,value = mean.gene.exp) %>%
      column_to_rownames("gene.name") %>%
      pheatmap:::scale_rows()
    tmp.plot[is.na(tmp.plot)] <- 0
  }else{
    ### I don't have any good idea by using sapply with dplyr
    tmp.plot.1 <- tmp.plot %>%
      select(gene.name,mean.gene.exp,MetaCell) %>%
      spread(key = MetaCell,value = mean.gene.exp) %>%
      column_to_rownames("gene.name")
    
    tmp.res.list <- sapply(
      idents.subset,
      FUN = function(ii){
        tmp.meta.1 <- tmp.meta %>%
          rownames_to_column("cell.id") %>%
          filter(cell.type %in% ii) %>%
          filter(sample %in% sample.subset) %>%
          arrange(cell.type,sample) %>%
          column_to_rownames("cell.id")
        
        tmp.base.index <- tmp.meta %>%
          rownames_to_column("cell.id") %>%
          filter(sample==base.sample) %>%
          filter(cell.type==ii) %>%
          column_to_rownames("cell.id")
        tmp.base.index <- rownames(tmp.base.index)
        
        if(pseudo.rep > 1){
          tmp.base <- rowMeans(tmp.plot.1[,tmp.base.index])
          tmp.base.sd <- ifelse(tmp.base==0,1,apply(tmp.plot.1[,tmp.base.index],1,sd))
        }else{
          tmp.base <- tmp.plot.1[,tmp.base.index]
          tmp.base.sd <- 1
        }
        
        gene.name <- rownames(tmp.plot.1)
        tmp.plot.1 <- tmp.plot.1[,rownames(tmp.meta.1)] %>%
          mutate_all(., ~(.x-tmp.base)/tmp.base.sd)
        # ###test distribution
        # summary(unlist(tmp.plot.1))
        tmp.plot.1[is.na(tmp.plot.1)] <- 0
        rownames(tmp.plot.1) <- gene.name
        return(tmp.plot.1)
      },
      simplify = F,
      USE.NAMES = T
    )
    
    tmp.plot <-  bind_cols(tmp.res.list) %>%
      as.data.frame %>%
      `row.names<-`(., row.names(tmp.res.list[[1]]))
  }
  
  
  
  
  ###test distribution
  ###summary(unlist(tmp.plot))
  ###using this to assure the right order
  tmp.plot <- tmp.plot[features,rownames(tmp.meta)]
  tmp.plot[tmp.plot > scale.max] = scale.max
  tmp.plot[tmp.plot < scale.min] = scale.min
  
  set.seed(42)
  ha <- HeatmapAnnotation(df = data.meta)
  tmp.color.maplist <- get_color_mapping_list(ha)
  set.seed(42)
  ha <- HeatmapAnnotation(df = tmp.meta,
                          col = list(cell.type=tmp.color.maplist[[1]]@colors[idents.subset],
                                     sample=tmp.color.maplist[[2]]@colors[sample.subset]))
  col_funs <- circlize::colorRamp2(seq(scale.min,
                                       scale.max,
                                       length.out = length(color.plot)),
                                   colors = color.plot)
  ht <- Heatmap(matrix = tmp.plot,
                name = ifelse(is.null(base.sample),"z-score","FC"), 
                show_column_names = show_columnnames,
                show_row_names = show_rownames,
                cluster_columns = cluster_columns,
                cluster_rows = cluster_rows,
                col = col_funs,
                use_raster = raster,
                show_row_dend = show_row_dend,
                heatmap_legend_param = list(at = pretty(c(min(tmp.plot)+0.5, 0, max(tmp.plot)-0.5),n=5)),
                row_names_gp = gpar(fontsize = 8), 
                top_annotation = ha,
                column_split = tmp.meta$cell.type,
                column_title = NULL,
                raster_device = raster.device)
  if(return.data){
    data.PseudoBulkHeatmap <- list(mat=tmp.plot,meta.mat=tmp.meta)
    return(list(ht=ht,data.PseudoBulkHeatmap=data.PseudoBulkHeatmap))
  }
  return(ht)
}


###--------MyCellPercent--------
#'@description  draw cell percents
MyCellPercent <- function(object=NULL,cells=NULL,
                          Idents.use="cell.type",
                          Idents.order=NULL,
                          MyTitle=NULL,
                          color=NULL){
  if(!is.null(cells)){
    object <- subset(object,cells = cells)
  }
  data.plot <- object[[]]
  test <- plyr::ddply(data.plot,
                      c("sample",Idents.use),
                      summarise,n=length(cell.type))
  
  test[[Idents.use]] <- factor(test[[Idents.use]],
                               levels = Idents.order)
  if(is.null(color)){
    ### get color palette
    data.meta <- object[[]] %>%
      select(cell.type,sample) 
    data.meta$cell.type <- factor(data.meta$cell.type,
                                  levels = levels(object))
    data.meta <- data.meta %>%
      arrange(cell.type,sample) 
    set.seed(42)
    ha <- HeatmapAnnotation(df = data.meta)
    tmp.color.maplist <- get_color_mapping_list(ha)
    color <- tmp.color.maplist[[1]]@colors
  }
  p <- ggplot(test,aes_string(x="sample", y="n",
                              fill=Idents.use))+
    geom_bar(stat = "identity",position = "fill")+
    theme_light()+
    scale_y_continuous(breaks = c(0,0.20,0.40,0.60,0.80,1.00),
                       labels = c("0","20","40","60","80","100"))+
    xlab("Sample ID")+
    ylab("Cell type percentage")+
    scale_fill_manual(values = color)+
    ggtitle(MyTitle)
  return(p)
}
