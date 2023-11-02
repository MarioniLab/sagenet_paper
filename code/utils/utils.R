# library('leiden')
# library('biomaRt')
library('spatstat')
library('igraph')
library('viridis')


celltype_colours = c("Epiblast" = "#635547",
                     "Primitive Streak" = "#DABE99",
                     "Caudal epiblast" = "#9e6762",
                     
                     "PGC" = "#FACB12",
                     
                     "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Definitive endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     "Gut tube" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#C9EBFB",
                     "Splanchnic mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Lateral plate mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     "Mixed mesenchymal mesoderm" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "Erythroid" = "#f79083",
                     "Blood progenitors" = "#f9decf",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     
                     "Visceral endoderm" = "#F6BFCB",
                     "ExE endoderm" = "#7F6874",
                     "ExE ectoderm" = "#989898",
                     "Parietal endoderm" = "#1A1A1A",
                     
                     "Unknown" = "#FFFFFF",
                     "Low quality" = "#e6e6e6",
                     
                     # somitic and paraxial types
                     # colour from T chimera paper Guibentif et al Developmental Cell 2021
                     "Cranial mesoderm" = "#77441B",
                     "Anterior somitic tissues" = "#F90026",
                     "Sclerotome" = "#A10037",
                     "Dermomyotome" = "#DA5921",
                     "Posterior somitic tissues" = "#E1C239",
                     "Presomitic mesoderm" = "#9DD84A"
)



.palette1   = c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", 
    "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A", 
    "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4", 
    "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", 
    "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
    )


.palette2   = c(
    "#FF7F00", "#FDB462", "#BEAED4",
    "#7BAFDE", "#E78AC3", "#B2DF8A",  
    "#B17BA6", "#E7298A", "#33A02C",
    "#55A1B1", "#A6761D", "#E6AB02", 
    "#8DD3C7", "#7570B3", "#8600bf")

.palette3 = c(
    "#666666", "#999999", "#aa8282", "#d4b7b7", 
    "#8600bf", "#ba5ce3", "#808000", "#aeae5c", 
    "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")


method_ord  = c('True Space', 'sagenet', 'tangram_cells', 'novosparc', 'sagenet_markers', 'exprs')
method_cols = c("#E1C239", '#33A02C', '#FF7F00', '#882E72', '#B2DF8A', '#666666')
names(method_cols) = method_ord

# method_ord  = c('True Space', 'sagenet', 'tangram_cells', 'novosparc', 'sagenet_markers', 'exprs')
# method_cols = c("#E1C239", '#1e90ff', '#f79083', '#B2DF8A', '#8DD3C7', '#BEAED4')
# names(method_cols) = method_ord

embryo_ord  = c('scRNAseq', 'embryo1_2', 'embryo1_5', 'embryo2_2', 'embryo2_5', 'embryo3_2', 'embryo3_5')
embryo_cols = c("#DC050C", "#8600bf", "#ba5ce3", "#E7298A", "#E78AC3", "#1965B0", "#7BAFDE")
names(embryo_cols) = embryo_ord

embryo_code = c(
    'scRNAseq'  = 'scRNAseq',
    'embryo1_2' = 'E1_L1', 
    'embryo1_5' = 'E1_L2', 
    'embryo2_2' = 'E2_L1', 
    'embryo2_5' = 'E2_L2',
    'embryo3_2' = 'E3_L1',
    'embryo3_5' = 'E3_L2')

method_code = c(
    'True Space'      = 'True Space',
    'sagenet'         = 'SageNet', 
    'tangram_cells'   = 'Tangram', 
    'novosparc'       = 'novoSpaRc', 
    'sagenet_markers' = "SageNet-SIG", 
    'exprs'           = 'Expression')

class__code = c(
    'Compact ventricular myocardium',
    'Trabecular ventricular myocardium (1)',
    'Trabecular ventricular myocardium (2)',
    'Trabecular ventricular myocardium (3)',
    'Atrial myocardium',
    'Outflow tract / large vessels',
    'Atrioventricular mesenchyme & valves',
    'Mediastinal mesenchyme & vessels',
    'Cavities with blood & immune cells',
    'Epicardium')

class__cols = c(
    'Atrial myocardium' = '#A10037',
    'Compact ventricular myocardium' = '#77441B',
    'Trabecular ventricular myocardium (1)' = '#aa8888',
    'Trabecular ventricular myocardium (2)' = '#826366',
    'Trabecular ventricular myocardium (3)' = '#5B3F45',
    'Outflow tract / large vessels' = '#DC050C',
    'Atrioventricular mesenchyme & valves' = '#cc7818',
    'Mediastinal mesenchyme & vessels' = '#EF4E22',
    'Cavities with blood & immune cells' = '#808000',
    'Epicardium' = '#0F4A9C')
class__ord = names(class__cols)

cell_type_cols = c(
 )


cell_type_ord = names(cell_type_cols)
# mart = useMart('ensembl', dataset='hsapiens_gene_ensembl', host='www.ensembl.org')

# go_cellcycle = getBM(
#     attributes = c('ensembl_gene_id','external_gene_name'), 
#     filters    = 'go', 
#     values     = 'GO:0007049', 
#     mart       = mart)

# go_translation = getBM(
#     attributes = c('ensembl_gene_id','external_gene_name'), 
#     filters    = 'go', 
#     values     = 'GO:0006412', 
#     mart       = mart)

# go_ribosome1 = getBM(
#     attributes = c('ensembl_gene_id','external_gene_name'), 
#     filters    = 'go', 
#     values     = 'GO:0005840', 
#     mart       = mart)

# go_ribosome2 = getBM(
#     attributes = c('ensembl_gene_id','external_gene_name'), 
#     filters    = 'go', 
#     values     = 'GO:0042254', 
#     mart       = mart)

# ex_genes = unique(c(
#     go_cellcycle$external_gene_name, 
#     go_translation$external_gene_name, 
#     go_ribosome1$external_gene_name,
#     go_ribosome2$external_gene_name)) 


nngraph_comm <- function(nngraph, min_cells = 100, res_param=0.1){
    comm_vec = leiden(nngraph, resolution_parameter=res_param)
    comm_dt  = table(comm_vec)
    to.keep  = max(as.numeric(names(comm_dt)[comm_dt > min_cells]))
    comm_vec = ifelse(comm_vec <= to.keep, comm_vec, NA)
    names(comm_vec) = names(V(nngraph))
    comm_vec
}

get_network <- function(data, rho = .1, threshold = .1){
    data     = data.table(data)
    S        = stats::cov.wt(data, method='ML')$cov
    C        = stats::cov2cor(S)
    res      = glasso::glasso(C, rho=rho)
    AM       = abs(res$wi) > threshold
    diag(AM) = F
    rownames(AM) = colnames(AM) = colnames(data)
    g.lasso  = graph_from_adjacency_matrix(AM)
    # names(N(g.lasso)) = colnames(data)
    as(g.lasso, 'igraph')
}

analyze_network <- function(graph_obj, res = 1){
    adj_mat = igraph::as_adj(graph_obj)
    comm  = leiden(adj_mat, resolution_parameter=res)
    isolated_comm = min(which(table(comm) == 1))
    names(comm) = colnames(adj_mat)
    comm = ifelse(comm >= isolated_comm, isolated_comm, comm)

    comm_dt = data.table(ID=names(comm), community=comm) %>%
        setkey(community) %>% 
        .[, color := c(.palette1, .palette2)[community]] %>%
        setkey(ID) %>%
        .[names(V(graph_obj))]

    comm_dt
    # markers = unique(comm) %>% 
    #     map(~names(base::sort(colSums(adj_mat[names(comm)[comm==.x], names(comm)[comm==.x]]), decreasing=T))[1:2]) %>%
    #     unlist

}

graph_col_comm <- function(graph, lay, grp, title=NULL, labels){
    igraph::V(graph)$color <- grp
    v <-  igraph::V(graph)
    # sprintf(comm_out, title) %>% pdf()
    p = plot.igraph(
        graph,
        vertex.size = 6,
        layout = lay,
        vertex.label = labels,
        vertex.frame.color = igraph::V(graph)$color,
        vertex.label.family = 'Helvetica',
        vertex.label.dist = 0,
        vertex.label.cex = .5,
        vertex.label.font = 0.1,
        vertex.label.color = '#695c59',
        main=title,
        edge.arrow.mode='-',
        directed=FALSE,
        edge.arrow.size=0)
    # dev.off()
    p
}

graph_col_act <- function(graph, values, contrast=1, lay, title){
    grp_range = c(min(values)^(contrast)/sum(values^(contrast)), max(values)^(contrast)/sum(values^(contrast)))
    grp_vals  = seq(grp_range[1],grp_range[2],length.out=9)
    grp_cols  = circlize::colorRamp2(grp_vals, viridis::viridis(9))
    igraph::V(graph)$color = grp_cols(values^(contrast)/sum(values^(contrast)))
    p = plot.igraph(graph,
        vertex.size = 5,
        layout = lay,
        vertex.frame.color = igraph::V(graph)$color,
        vertex.label = "",
        main=title)
    p
}

plot_spatial <- function(dim_df, labels, label_cols=c(.palette1, .palette2, .palette3), title='', label_title='label', hide_legend=TRUE, sz=0.5){
    dim_dt = data.table(label=labels,
                         dim1=unlist(dim_df[,1]), dim2=unlist(dim_df[,2])) 
    dim_plot = dim_dt %>%
        ggplot +
        aes(dim1, dim2, color=label) +
        geom_point(size=sz) +
        theme_minimal() + 
        theme(axis.text= element_blank(), 
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=label_cols, na.value='#e6e6e6' , drop = TRUE, limits = unique(labels)) +
        labs(title=title, x='', y='', color=label_title) 
    if(hide_legend)
        dim_plot = dim_plot + theme(legend.position='none')
    dim_plot
}

plot_2d <- function(dim_df, labels, label_cols=c(.palette1, .palette2, .palette3), title='', label_title='label', hide_legend=TRUE, sz=1, alpha=1, shape=16){
    dim_dt = data.table(label=labels,
                         dim1=unlist(dim_df[,1]), dim2=unlist(dim_df[,2]), alpha=alpha) 
    dim_plot = dim_dt %>%
        ggplot +
        aes(dim1, dim2, color=label) +
        geom_point(size=sz, alpha=alpha, shape=shape) +
        theme_minimal() + 
        theme(axis.text= element_blank(), 
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=label_cols, na.value='#e6e6e6' , drop = TRUE, limits = unique(labels)) +
        labs(title=title, x='', y='', color=label_title, alpha='confidence score')
    if(hide_legend)
        dim_plot = dim_plot + theme(legend.position='none')
    dim_plot
}


plot_2d_cont <- function(dim_df, labels, label_cols=nice_cols, title='', label_title='label', sz=3, hide_legend=TRUE, shape=16){
    dim_dt = data.table(label=labels,
                         dim1=unlist(dim_df[,1]), dim2=unlist(dim_df[,2])) 
    dim_plot = dim_dt %>%
        ggplot +
        aes(dim1, dim2, color=label) +
        # geom_hex(bins = 30) + 
        geom_point(size=sz, shape=shape) +
        # coord_fixed() +
        scale_color_viridis(na.value='#e6e6e6') +
        theme_bw() + 
        theme(
            axis.text= element_blank(), 
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
        labs(title=title, x='', y='', color=label_title) 
        if(hide_legend)
            dim_plot = dim_plot + theme(legend.position='none')
    dim_plot
}


read_preds_cells <- function(preds_f, func='mean'){
    preds     = preds_f %>% readH5AD
    ref_locs  = rowData(preds)[,c('x', 'y')] %>% as.matrix
    preds_p   = t(assays(preds)[['X']])
    if(func == 'max')
        preds_p %<>% 
            apply(1, function(x) ifelse(x == max(x), 1, 0)) %>%
            t 
    preds_p %*% ref_locs
}

read_preds_novosparc <- function(preds_f, func='mean'){
    preds     = zellkonverter::readH5AD(preds_f)
    ref_locs  = rowData(preds)[,c('x', 'y')] %>% as.matrix
    preds_p   = t(assays(preds)[['X']])
    if(func == 'max')
        preds_p %<>% 
            apply(1, function(x) ifelse(x == max(x), 1, 0)) %>%
            t 
    preds_p %*% ref_locs
}

read_preds_tangram <- function(preds_f, func='mean'){
    preds     = anndata::read_h5ad(preds_f)
    ref_locs  = preds$obs[,c('x', 'y')] %>% as.matrix
    preds_p   = t(preds$X)
    if(func == 'max')
        preds_p %<>% 
            apply(1, function(x) ifelse(x == max(x), 1, 0)) %>%
            t 
    preds_p %*% ref_locs
}


read_preds_cluster <- function(preds_f){
    print(preds_f)
    preds     = preds_f %>% readH5AD
    preds_p   = t(assays(preds)[['X']])
    preds_p
}

read_preds <- function(preds_f){
  preds = preds_f %>% fread %>% .[-1,-1] %>% as.matrix %>%
    apply(1, function(x) exp(x)/sum(exp(x))) %>%
    t
  colnames(preds) = paste('label', 1:dim(preds)[2], sep='_')

  preds_class = apply(preds, 1, function(x) ifelse(max(x) <0.5, NA, colnames(preds)[which.max(x)]))

  list(preds=preds, preds_class=preds_class)
}

# read_preds_tangram <- function(preds_f){
#   preds = preds_f %>% fread %>% .[-1,-1] %>% as.matrix 
#   colnames(preds) = paste('label', 1:dim(preds)[2], sep='_')

#   preds_class = apply(preds, 1, function(x) ifelse(max(x) <0, NA, colnames(preds)[which.max(x)]))

#   list(preds=preds, preds_class=preds_class)
# }


read_imp <- function(imp_f){
  imp = imp_f %>% fread %>% setnames('V1', 'GENE')
  genes = imp$GENE
  imp %<>% .[,-1] %>% as.matrix
  rownames(imp) = genes
  imp_abs = abs(imp)
  # imp = apply(imp, 2, function(x) (x)/(max(x)))
  imp_dt = apply(imp, 2, function(x) rownames(imp)[order(-x)[1:n_markers]]) %>%
      as.data.table
  imp_dt

}


get_confidence <- function(sce){
    meta_dt = colData(sce) %>% as.data.table 
    meta_dt %<>% .[, .SD, .SDcols = names(meta_dt) %like% '^ent']
    rowSums(meta_dt)
}


