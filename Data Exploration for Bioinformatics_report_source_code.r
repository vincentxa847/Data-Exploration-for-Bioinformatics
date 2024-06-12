library(ggplot2)
library(dplyr)
library(org.Hs.eg.db) # Annotation for human, data was collected from IMR-90 cell line, so using database of human 


#### File handle ####
File_Handle = function(file_path)
{
  table = read.table(file_path,header=TRUE, row.names=1, sep= "\t")
  return(table)
}

DE_Senes_MtD_vs_Prolif = File_Handle("./Data/DE_Senes_MtD_vs_Prolif.csv")
DE_Senes_MtD_vs_Senes = File_Handle("./Data/DE_Senes_MtD_vs_Senes.csv")
DE_Senes_vs_Prolif = File_Handle("./Data/DE_Senes_vs_Prolif.csv")
Expression_table = File_Handle("./Data/EM.csv")
Gene_Annotation = File_Handle("./Data/Human_Background_GRCh38.p13.csv")
Sample_sheet = File_Handle("./Data/sample_sheet.csv")

#### Theme ####
# Making custom theme
my_theme = ggplot2::theme(
  plot.title = element_text(size=rel(1.1),face="bold"),
  panel.background = element_rect(fill = "white"),
  axis.line = element_line(size = rel(1), linetype=1),
  legend.position = "right",
  legend.key.size = unit(0.5, 'cm')
)
color_for_upregulate = "#701820" # the color for sig_up and down in MA plot
color_for_downregulate = "#1948A4"
color_for_no_change = "#D8D7DF"
color_for_Prolif = "#283870"
color_for_Senes = "#6E3F2B"
color_for_Senes_MtD = "#41553A" 

# Making "master" table for plot (select significant differential genes, add gene information, add mean, deal with row name) 
Select = function(de_table,p.adj_thresholds,log2fold_threshold)
{
  # merge table and select significant differential genes
  de_table = merge(de_table,Gene_Annotation, by.x = 0, by.y = 0, header = TRUE ) # Adding annotation table
  de_table = merge(Expression_table,de_table, by.x = 0, by.y = 1,header = TRUE ) # Adding expression table
  row.names(de_table) = de_table[,1]
  de_table = de_table[,-1]
  de_table$ mean = rowMeans(de_table[,1:9])
  de_table_sig = subset(de_table, p.adj < p.adj_thresholds & abs(log2fold) > log2fold_threshold)
  
  # sort by p.adj value
  sort_order = order(de_table_sig[,"p.adj"],decreasing = FALSE)
  de_table_sig = de_table_sig[sort_order,]
  
  return(de_table_sig)
  
}

Sig_Senes_MtD_vs_Prolif = Select(DE_Senes_MtD_vs_Prolif,0.001,2)
Sig_Senes_MtD_vs_Prolif_UP = subset(Sig_Senes_MtD_vs_Prolif,log2fold > 0)
Sig_Senes_MtD_vs_Prolif_DOWN = subset(Sig_Senes_MtD_vs_Prolif,log2fold < 0)

Sig_Senes_MtD_vs_Senes = Select(DE_Senes_MtD_vs_Senes,0.001,2)
Sig_Senes_MtD_vs_Senes_UP = subset(Sig_Senes_MtD_vs_Senes,log2fold > 0)
Sig_Senes_MtD_vs_Senes_DOWN = subset(Sig_Senes_MtD_vs_Senes,log2fold < 0)

Sig_Senes_vs_Prolif = Select(DE_Senes_vs_Prolif,0.001,2)
Sig_Senes_vs_Prolif_UP = subset(Sig_Senes_vs_Prolif,log2fold > 0)
Sig_Senes_vs_Prolif_DOWN = subset(Sig_Senes_vs_Prolif,log2fold < 0)

#### Volcano plot ####
plot_volcano = function(de_table,de_table_sig_UP,de_table_sig_DOWN,plot_title,p_threshold,fold_threshold)
{
  # remove rows contain NA
  de_table = na.omit(de_table) 
  
  # sort by p.adj value
  sort_order__for_up = order(de_table_sig_UP[,"p.adj"],decreasing = FALSE)
  de_table_sig_UP = de_table_sig_UP[sort_order__for_up,]
  sort_order__for_down = order(de_table_sig_DOWN[,"p.adj"],decreasing = FALSE)
  de_table_sig_DOWN = de_table_sig_DOWN[sort_order__for_down,]
  
  # pick out the top 5 
  de_table_sig_up_top5 = de_table_sig_UP[1:5,]
  print(de_table_sig_up_top5)
  de_table_sig_down_top5 = de_table_sig_DOWN[1:5,]
  print(de_table_sig_down_top5)
  # making plot
  result = ggplot(de_table, aes(x=log2fold, y= -log10(p.adj))) + 
    geom_point(aes(colour="a")) + 
    geom_point(data= de_table_sig_UP, aes(colour="b")) + 
    geom_point(data= de_table_sig_DOWN, aes(colour="c")) +
    xlim(-10,10) + # set scale for x axis to make sure consistency between 3 groups when plotting
    my_theme +
    geom_vline(xintercept=-2,linetype="dashed") +
    geom_vline(xintercept=2,linetype="dashed") +
    geom_hline(yintercept=-log10(0.001),linetype="dashed") + 
    labs(title=plot_title, x="Log2 fold change", y="-Log10 p") + 
    scale_color_manual(values=c(color_for_no_change , color_for_upregulate, color_for_downregulate), labels = c("No change", "Up", "Down"), name="") +
    ggrepel::geom_text_repel(data= de_table_sig_up_top5, aes(label=SYMBOL,colour = "b"),size = 2.5,fontface = 2, show.legend=FALSE)+ # changing label size to avoid OVERLAP, if label comes from different group, then repel cannot work functionally.
    ggrepel::geom_text_repel(data= de_table_sig_down_top5, aes(label=SYMBOL,colour = "c"),size = 2.5,fontface = 2, show.legend=FALSE) # changing fontface to bold(2) to make it clearly
  
  
  return(result)
}
Volcano_Plot_For_Senes_MtD_vs_Prolif = plot_volcano(DE_Senes_MtD_vs_Prolif,Sig_Senes_MtD_vs_Prolif_UP,Sig_Senes_MtD_vs_Prolif_DOWN,"Senes_MtD_vs_Prolif",0.001,2) # total sig:1798, upregulate: 907, downregulate: 891
Volcano_Plot_For_Senes_MtD_vs_Prolif
Volcano_Plot_For_Senes_MtD_vs_Senes = plot_volcano(DE_Senes_MtD_vs_Senes,Sig_Senes_MtD_vs_Senes_UP,Sig_Senes_MtD_vs_Senes_DOWN,"Senes_MtD_vs_Senes",0.001,2) # total sig:1572, upregulate: 466, downregulate: 1106
Volcano_Plot_For_Senes_MtD_vs_Senes
Volcano_Plot_For_Senes_vs_Prolif = plot_volcano(DE_Senes_vs_Prolif,Sig_Senes_vs_Prolif_UP,Sig_Senes_vs_Prolif_DOWN,"Senes_vs_Prolif",0.001,2) # total sig:1291, upregulate: 862, downregulate: 429
Volcano_Plot_For_Senes_vs_Prolif 

#### MA plot ####
plot_MA = function(de_table,de_table_sig,plot_title,p_threshold,fold_threshold)
{
  # add mean for de_table
  de_table = na.omit(de_table) # remove rows contain NA
  de_table = merge(Expression_table,de_table, by.x = 0, by.y = 0,header = TRUE ) # Adding expression table
  de_table$ mean = rowMeans(de_table[,2:10])
  
  # select significant differential genes (up-regulate and down-regulate)
  de_table_sig_up = subset(de_table_sig,log2fold > 0)
  de_table_sig_down = subset(de_table_sig,log2fold < 0)
  
  # sort by p.adj value
  sort_order__for_up = order(de_table_sig_up[,"p.adj"],decreasing = FALSE)
  de_table_sig_up = de_table_sig_up[sort_order__for_up,]
  sort_order__for_down = order(de_table_sig_down[,"p.adj"],decreasing = FALSE)
  de_table_sig_down = de_table_sig_down[sort_order__for_down,]
  
  # pick out the top 5 
  de_table_sig_up_top5 = de_table_sig_up[1:5,]
  de_table_sig_down_top5 = de_table_sig_down[1:5,]
  
  # making plot
  result = ggplot(de_table, aes(x=log10(mean), y= log2fold)) + 
    geom_point(colour=color_for_no_change) + 
    geom_point(data= de_table_sig_up, colour=color_for_upregulate) + 
    geom_point(data= de_table_sig_down, colour=color_for_upregulate) +
    my_theme +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept=-2,linetype="dashed") +
    geom_hline(yintercept=2,linetype="dashed") +
    labs(title=plot_title, x="mean expression (log10)", y="log2fold") +
    ggrepel::geom_text_repel(data= de_table_sig_up_top5, aes(label=SYMBOL),colour = color_for_upregulate,size = 3.2,fontface = 2) + # changing label size to avoid ggrepel warning : unlabeled data points (too many overlaps).
    ggrepel::geom_text_repel(data= de_table_sig_down_top5, aes(label=SYMBOL),colour = color_for_upregulate,size = 3.2,fontface = 2)
}

MA_Plot_For_Senes_MtD_vs_Prolif = plot_MA(DE_Senes_MtD_vs_Prolif,Sig_Senes_MtD_vs_Prolif,"Senes_MtD_vs_Prolif",0.001,2) # total sig 1798
MA_Plot_For_Senes_MtD_vs_Prolif
MA_Plot_For_Senes_MtD_vs_Senes = plot_MA(DE_Senes_MtD_vs_Senes,Sig_Senes_MtD_vs_Senes,"Senes_MtD_vs_Senes",0.001,2) #  total sig 1106
MA_Plot_For_Senes_MtD_vs_Senes
MA_Plot_For_Senes_vs_Prolif = plot_MA(DE_Senes_vs_Prolif,Sig_Senes_vs_Prolif,"Senes_vs_Prolif",0.001,2) # total sig 1291
MA_Plot_For_Senes_vs_Prolif

#### PCA plot ####
plot_PCA = function(group_expression)
{
  # just take expression table for PCA
  if ("SYMBOL" %in% colnames(group_expression)){
    group_expression = group_expression[,1:9]
  }
  em_scaled = as.matrix(sapply(group_expression,as.numeric))
  pca = prcomp(t(em_scaled))
  pca_coordinates = data.frame(pca$x)
  pca_coordinates$group = Sample_sheet$SAMPLE_GROUP # Adding group for colour
  
  # adding the % variance
  vars = apply(pca$x, 2, var) # var function apply to column
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
  y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")
  
  ggp = ggplot(pca_coordinates,aes(x=PC1, y=PC2, colour = group)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label=row.names(pca_coordinates)),show.legend = FALSE,size = 2) +
    scale_colour_manual(values = c(color_for_Prolif, color_for_Senes,color_for_Senes_MtD)) +
    xlab(x_axis_label) +
    ylab(y_axis_label) 
  
  
  return(ggp)
}

PCA_for_whole_geneset = plot_PCA(Expression_table)
PCA_for_whole_geneset
PCA_for_Sig_Senes_MtD_vs_Prolif = plot_PCA(Sig_Senes_MtD_vs_Prolif)
PCA_for_Sig_Senes_MtD_vs_Prolif
PCA_for_Sig_Senes_MtD_vs_Senes = plot_PCA(Sig_Senes_MtD_vs_Senes)
PCA_for_Sig_Senes_MtD_vs_Senes
PCA_for_Sig_Senes_vs_Prolif =  plot_PCA(Sig_Senes_vs_Prolif)
PCA_for_Sig_Senes_vs_Prolif # the distance between group in Senes is much more apart from each other 

#### Expression density plot ####
plot_density = function(group_expression)
{
  group_expression.m = reshape2::melt(Expression_table)
  
  ggp = ggplot(data = group_expression.m, aes(x=log10(value),fill = variable)) +
    geom_density(show.legend = FALSE) +
    facet_wrap(~variable,nrow=3) +
    # using same colour for data in same comparing group
    scale_fill_manual(values=c(color_for_Prolif,color_for_Prolif,color_for_Prolif, color_for_Senes, color_for_Senes, color_for_Senes , color_for_Senes_MtD, color_for_Senes_MtD, color_for_Senes_MtD)) 
  
  return(ggp)
}

Density_plot = plot_density(Expression_table)
Density_plot

#### Multiple boxplot  ####
make_multi_boxplot = function(candidate_genes, sample_groups,plot_title)
{
  gene_data = Expression_table[candidate_genes,] # input a list of gene (candidate_genes)
  row.names(gene_data) = Gene_Annotation[row.names(gene_data),"SYMBOL"]
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = sample_groups
  gene_data.m = reshape2::melt(gene_data, id.vars="sample_group")
  
  # change the setting of title to better fit in combine plot
  my_theme = theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_blank(),
                   axis.title.y = element_blank(),plot.title = element_text(size = rel(0.8),hjust = 0.5))
  
  ggp = ggplot(gene_data.m,aes(x=variable,y=log10(value),colour = sample_group)) +
    geom_boxplot() +
    scale_colour_manual(values = c(color_for_Prolif, color_for_Senes,color_for_Senes_MtD)) +
    ylim(0,6) +
    labs(title = plot_title) +
    my_theme
  
  return(ggp)
}

candidate_genes = row.names(Sig_Senes_MtD_vs_Prolif[1:10,])
Multiple_boxplot_Sig_Senes_MtD_vs_Prolif_TOP10 = make_multi_boxplot(candidate_genes, Sample_sheet$SAMPLE_GROUP,"Sig_Senes_MtD_vs_Prolif_TOP10")
Multiple_boxplot_Sig_Senes_MtD_vs_Prolif_TOP10

candidate_genes = row.names(Sig_Senes_MtD_vs_Senes[1:10,])
Multiple_boxplot_Sig_Senes_MtD_vs_Senes_TOP10 = make_multi_boxplot(candidate_genes, Sample_sheet$SAMPLE_GROUP,"Sig_Senes_MtD_vs_Senes_TOP10")
Multiple_boxplot_Sig_Senes_MtD_vs_Senes_TOP10

candidate_genes = row.names(Sig_Senes_vs_Prolif[1:10,])
Multiple_boxplot_Sig_Senes_vs_Prolif_TOP10 = make_multi_boxplot(candidate_genes, Sample_sheet$SAMPLE_GROUP,"Sig_Senes_vs_Prolif_TOP10")
Multiple_boxplot_Sig_Senes_vs_Prolif_TOP10

#### Heatmap ####
# Making heatmap for sig genes in different groups
plot_heatmap = function(expression_table_sig,plot_title)
{
  em_sig = expression_table_sig[,1:9]
  em_sig_scaled = data.frame(t(scale(t(em_sig))))
  hm.matrix = as.matrix(em_sig_scaled)
  
  # set gradient colour
  colours = c("#1948A4","#F0D888","#701820")
  colorRampPalette(colours)(100)
  
  # clustering
  y.dist = amap::Dist(hm.matrix, method="spearman")
  
  y.cluster = hclust(y.dist, method="average")
  
  y.dd = as.dendrogram(y.cluster)
  
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  
  y.order = order.dendrogram(y.dd.reorder)
  
  hm.matrix_clustered = hm.matrix[y.order,]
  
  hm.matrix_clustered = as.matrix(hm.matrix_clustered)
  hm.matrix_clustered = reshape2::melt(hm.matrix_clustered) 
  
  ggp = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + 
    scale_fill_gradientn(colours = colorRampPalette(colours)(100)) +
    labs(title=plot_title, x="", y="") + 
    
    # remove redundant y text, also create some space for x text
    theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45,vjust = 0.6),axis.ticks=element_blank()) 
  
  
  return(ggp)
}
Heat_map_Sig_Senes_MtD_vs_Prolif = plot_heatmap(Sig_Senes_MtD_vs_Prolif,"Sig_Senes_MtD_vs_Prolif")
Heat_map_Sig_Senes_MtD_vs_Prolif
Heat_map_Sig_Senes_MtD_vs_Senes = plot_heatmap(Sig_Senes_MtD_vs_Senes,"Sig_Senes_MtD_vs_Senes")
Heat_map_Sig_Senes_MtD_vs_Senes
Heat_map_Sig_Senes_vs_Prolif = plot_heatmap(Sig_Senes_vs_Prolif,"Sig_Senes_vs_Prolif")
Heat_map_Sig_Senes_vs_Prolif

Combined_Plot = ggpubr::ggarrange(Heat_map_Sig_Senes_MtD_vs_Prolif,Heat_map_Sig_Senes_MtD_vs_Senes,Heat_map_Sig_Senes_vs_Prolif,
                          ncol=3,labels = "AUTO")
Combined_Plot

#### Pathway analysis ####
pathway_analysis = function(expression_table_sig,Title)
{
  sig_genes = expression_table_sig[,"SYMBOL"]
  # Biological Id TRanslator to ENTREZID
  sig_genes_entrez = clusterProfiler::bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  pathway_results = clusterProfiler::enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, 
                             ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  # using clusterProfiler in-built plotting tool
  ggp = clusterProfiler::dotplot(pathway_results, showCategory= 10,title= Title)

  return(ggp)
}
Pathway_analysis_for_Sig_Senes_MtD_vs_Prolif = pathway_analysis(Sig_Senes_MtD_vs_Prolif,"Sig_Senes_MtD_vs_Prolif")
Pathway_analysis_for_Sig_Senes_MtD_vs_Prolif
Pathway_analysis_for_Sig_Senes_MtD_vs_Prolif_UP = pathway_analysis(Sig_Senes_MtD_vs_Prolif_UP,"Sig_Senes_MtD_vs_Prolif_UP")
Pathway_analysis_for_Sig_Senes_MtD_vs_Prolif_UP
Pathway_analysis_for_Sig_Senes_MtD_vs_Prolif_DOWN = pathway_analysis(Sig_Senes_MtD_vs_Prolif_DOWN,"Sig_Senes_MtD_vs_Prolif_DOWN")
Pathway_analysis_for_Sig_Senes_MtD_vs_Prolif_DOWN

Pathway_analysis_for_Sig_Senes_MtD_vs_Senes = pathway_analysis(Sig_Senes_MtD_vs_Senes,"Sig_Senes_MtD_vs_Senes")
Pathway_analysis_for_Sig_Senes_MtD_vs_Senes
Pathway_analysis_for_Sig_Senes_MtD_vs_Senes_UP = pathway_analysis(Sig_Senes_MtD_vs_Senes_UP,"Sig_Senes_MtD_vs_Senes_UP")
Pathway_analysis_for_Sig_Senes_MtD_vs_Senes_UP
Pathway_analysis_for_Sig_Senes_MtD_vs_Senes_DOWN = pathway_analysis(Sig_Senes_MtD_vs_Senes_DOWN,"Sig_Senes_MtD_vs_Senes_DOWN")
Pathway_analysis_for_Sig_Senes_MtD_vs_Senes_DOWN


Pathway_analysis_for_Sig_Senes_vs_Prolif = pathway_analysis(Sig_Senes_vs_Prolif,"Sig_Senes_vs_Prolif")
Pathway_analysis_for_Sig_Senes_vs_Prolif
Pathway_analysis_for_Sig_Senes_vs_Prolif_UP = pathway_analysis(Sig_Senes_vs_Prolif_UP,"Sig_Senes_vs_Prolif_UP")
Pathway_analysis_for_Sig_Senes_vs_Prolif_UP
Pathway_analysis_for_Sig_Senes_vs_Prolif_DOWN = pathway_analysis(Sig_Senes_vs_Prolif_DOWN,"Sig_Senes_vs_Prolif_DOWN")
Pathway_analysis_for_Sig_Senes_vs_Prolif_DOWN

Combined_Plot = ggpubr::ggarrange(Pathway_analysis_for_Sig_Senes_MtD_vs_Prolif,Pathway_analysis_for_Sig_Senes_MtD_vs_Senes,Pathway_analysis_for_Sig_Senes_vs_Prolif,
                          ncol=3,labels = "AUTO")
Combined_Plot

#### Signature ####
Add_mean_for_three_groups_to_find_signature = function(table)
{
  table$Prolif_mean = rowMeans(table[,1:3])
  table$Senes_mean = rowMeans(table[,4:6])
  table$Senes_MtD_mean = rowMeans(table[,7:9])
  return (table)
}

metagene_boxplot = function(signature)
{
  signature_expression = data.frame(t(scale(t(signature[,1:9])))) # scale the expression
  metagene = data.frame(colMeans(signature_expression)) 
  names(metagene) = "value"
  metagene$group = Sample_sheet$SAMPLE_GROUP
  
  ggp = ggplot(metagene, aes(x=group, y= value, colour= group)) + geom_boxplot() +
    scale_colour_manual(values = c(color_for_Prolif,color_for_Senes,color_for_Senes_MtD))
  
  return(ggp)
}

Up_in_Senes_Down_in_Senes_MtD_and_Prolif = Add_mean_for_three_groups_to_find_signature(Sig_Senes_MtD_vs_Prolif)
Up_in_Senes_Down_in_Senes_MtD_and_Prolif = subset(Up_in_Senes_Down_in_Senes_MtD_and_Prolif,Senes_mean > Senes_MtD_mean & Senes_mean > Prolif_mean)
metagene_boxplot_Up_in_Senes_Down_in_Senes_MtD_and_Prolif = metagene_boxplot(Up_in_Senes_Down_in_Senes_MtD_and_Prolif)
metagene_boxplot_Up_in_Senes_Down_in_Senes_MtD_and_Prolif
Pathway_analysis_for_Up_in_Senes_Down_in_Senes_MtD_and_Prolif = pathway_analysis(Up_in_Senes_Down_in_Senes_MtD_and_Prolif,"Up_in_Senes_Down_in_Senes_MtD_and_Prolif") + ggtitle("") # remove the title (to avoid overlap when combining plot)
Pathway_analysis_for_Up_in_Senes_Down_in_Senes_MtD_and_Prolif # extracellular matrix organization

#### Venn Diagram ####
genelist_for_MtD_UP = list(Prolif = Sig_Senes_MtD_vs_Prolif_UP[,"SYMBOL"], Senes = Sig_Senes_MtD_vs_Senes_UP[,"SYMBOL"])
VennDiagram_MtD = ggVennDiagram::ggVennDiagram (genelist_for_MtD_UP) +
  scale_fill_gradient(high = color_for_upregulate) 
VennDiagram_MtD # THE NAMES ON THE DIAGRAM SHOW THE "DIFFERENT" COMPARSION BETWEEN THE GROUP

genelist_for_Senes = list(Senes_MtD = Sig_Senes_MtD_vs_Prolif[,"SYMBOL"], Senes = Sig_Senes_vs_Prolif[,"SYMBOL"])
VennDiagram_Senes = ggVennDiagram::ggVennDiagram(genelist_for_Senes) +
  scale_fill_gradient(high = color_for_upregulate) 
VennDiagram_Senes
