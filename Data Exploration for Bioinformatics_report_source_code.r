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
