library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(amap)
library(clusterProfiler)
library(org.Hs.eg.db) # Annotation for human, data was collected from IMR-90 cell line, so using database of human 
library(ggVennDiagram)

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
