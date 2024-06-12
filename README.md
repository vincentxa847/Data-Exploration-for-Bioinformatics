# Using bulk RNA-seq dataset to compare the transcriptomes of proliferating, replicative senescent, and mitochondria-depleted senescent fibroblasts

## Result
Genes that show significant differential expression (p.adj < 0.001, absolute log2fold > 2) in three comparisons (Senes_MtD_vs_Prolif, Senes_MtD_vs_Senes and Senes_vs_Prolif) were selected separately
(with prefix Sig-) for differential gene expression analysis. 
```
#  Function to generate a "master" table that selects differential genes.
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

# The thresholds are set at an adjusted p-value cutoff of 0.001 and an absolute log2 fold change of 2
Sig_Senes_MtD_vs_Prolif = Select(DE_Senes_MtD_vs_Prolif,0.001,2)
Sig_Senes_MtD_vs_Senes = Select(DE_Senes_MtD_vs_Senes,0.001,2)
Sig_Senes_vs_Prolif = Select(DE_Senes_vs_Prolif,0.001,2)
```
MA plots indicate differential genes evenly distribute relative to y axis and tightens as the x axis increase. This means the normalization process of RNA-Seq data assigned is appropriate to avoid
deviate from y=0, also the sequencing depth is appropriate to prevent generally low x value. The result of MA plots shows a good quality of sequencing and processing before analysis.
![MA_plot](https://github.com/vincentxa847/Data-Exploration-for-Bioinformatics/assets/118545004/fec85a33-16b4-4732-bbe3-ca8fd72a6f9f)
*MA plot for quality control of dataset*

Principal component analysis and expression density plot show a consistency for replicates within group, which can generally exclude the possibility of technical issues when preparing samples.
```
#  Function to perform PCA and display the result
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
```
![Density plot and PCA plot](https://github.com/vincentxa847/Data-Exploration-for-Bioinformatics/assets/118545004/fb50b784-6a25-4dea-b304-15a6383d1321)
*PCA and density plot for quality control of dataset*

In *Senes_MtD_vs_Prolif group, PGK1, SLC2A3, ENO2, NDRG1 and PFKFB4 with most significant up-regulation*. 
[PGK1 plays a key role in glycolysis by transferring the phosphoryl-group from 1,3-biophosphoglycerate to ADP, which creates ATP and 3-phosphoglycerate](https://www.nature.com/articles/385275a0).
[SLC2A3 belongs to solute carrier family (SLCs) that localized on membrane to transports glucose and nutrients across membrane](https://doi.org/10.1016/j.cell.2015.07.022). 
ENO2 is an enolase that takes part in the ninth step of glycolysis by reversible converting 2-phosphoglycerate to phosphoenolpyruvate. 
[NDRG1 is a regulator of proliferation, which performs by Fe chelation](https://doi.org/10.1182/blood-2004-05-1866). 
[PFKFB4 is a bifunctional enzyme that transfers the phosphoryl-group from ATP to fructose-6-phosphste and produce fructose-2,6-biphosphate](https://doi.org/10.1038/s41586-018-0018-1).\
*Most significant up-regulation genes in Senes_MtD_vs_Senes are same as in Senes_MtD_vs_Prolif*, it can be inferred that depletion of mitochondria leads to the up-regulation of PGK1,
SLC2A3, ENO2 and PFKFB4 and therefore the stimulation of glycolysis.
```
# Function to generate volcano plot displaying relevant information 
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
```
![Volcano Plot](https://github.com/vincentxa847/Data-Exploration-for-Bioinformatics/assets/118545004/6636ee2f-0cde-4f80-89b9-8a7e9f218e40)
*Volcano plot of differential expression genes in three comparsions*

Pathway analysis was performed to investigate the function of significant differential genes. It can be seen that the genes with up-regulation in
Sig_Senes_MtD_vs_Senes and Sig_Senes_MtD_vs_Prolif are relevant to oxygen levels.
![Sig_Senes_MtD_UP_Pathway](https://github.com/vincentxa847/Data-Exploration-for-Bioinformatics/assets/118545004/215e4265-d991-4dfa-a316-5f67c1591c65)
*Over Representation Analysis (ORA) of genes with up-regulation in comparsions Sig_Senes_MtD_vs_Senes and Sig_Senes_MtD_vs_Prolif*

Genes with down-regulation in Sig_Senes_MtD_vs_Prolif are relevant to the function of cell division, which is same as in Sig_Senes_vs_Prolif. It can be inferred that the rate of cell division decreases in replicative senescent cells, and the depletion of mitochondria is unable to rescue. KIF gene family that shows most significant down-regulation
in both Senes_MtD_vs_Prolif and Senes_vs_Prolif is [for intracellular transport, which is important to cell division](https://doi.org/10.1073/pnas.111145398).
![Sene_DOWN_pathway](https://github.com/vincentxa847/Data-Exploration-for-Bioinformatics/assets/118545004/a70cda61-986d-42be-8971-bb5840c9eb58)
*Over Representation Analysis (ORA) of genes with udown-regulation in comparsions Sig_Senes_vs_Prolif and Sig_Senes_MtD_vs_Prolif*

Clustered heatmaps were made to investigate the expression pattern of significant differential genes in three sample groups. It is worth noticed that part of the genes in Sig_Senes_MtD_vs_Prolif show up-regulation in Senes in relative to Senes_MtD and Prolif. From the result of pathway analysis, it can be inferred that mitochondria depletion
attenuates the function of replicative senescent cells for extracellular matrix organization.
```
# Function for gene Clustering and making heatmap 
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
```
![Cluster_Heatmap](https://github.com/vincentxa847/Data-Exploration-for-Bioinformatics/assets/118545004/07dccb31-1d72-4282-94b5-37e1afb6e9b2)\
*Expression profiles of differentially expressed RNA in Senes_MtD_vs_Prolif, Senes_MtD_vs_Senes and Senes_vs_Prolif*

Venn Diagrams were plotted to investigate the number of genes that overlap in different comparisons.
308 genes overlap in upregulation part of Sig_Senes_MtD_vs_Prolif and Sig_Senes_MtD_vs_Senes, which means that mitochondria depletion makes IMR-90 cells up-regulate some genes to adapt to the change. This change seems to be in a replicative senescent independent manner.
551 genes overlap in Sig_Senes_MtD_vs_Prolif and Sig_Senes_vs_Prolif. It can be inferred that some features of replicative senescent are in a mitochondria independent manner, and
this is consistent with the result of pathway analysis.
```
# Venn Diagram 
genelist_for_MtD_UP = list(Prolif = Sig_Senes_MtD_vs_Prolif_UP[,"SYMBOL"], Senes = Sig_Senes_MtD_vs_Senes_UP[,"SYMBOL"])
VennDiagram_MtD = ggVennDiagram::ggVennDiagram (genelist_for_MtD_UP) + scale_fill_gradient(high = color_for_upregulate) 
VennDiagram_MtD # THE NAMES ON THE DIAGRAM SHOW THE "DIFFERENT" COMPARSION BETWEEN THE GROUP

genelist_for_Senes = list(Senes_MtD = Sig_Senes_MtD_vs_Prolif[,"SYMBOL"], Senes = Sig_Senes_vs_Prolif[,"SYMBOL"])
VennDiagram_Senes = ggVennDiagram::ggVennDiagram(genelist_for_Senes) +
  scale_fill_gradient(high = color_for_upregulate) 
VennDiagram_Senes
```
![Venn Diagram](https://github.com/vincentxa847/Data-Exploration-for-Bioinformatics/assets/118545004/ebabfeb3-97e9-4af2-b73c-29fe3e4b245c)\
*Analysis of differential expression signatures*

## Conclusion
Replicative senescent cells show a down regulation of genes related to cell division, and up-regulation of genes related to extracellular matrix organization. Mitochondria depletion results in a down regulation of genes related to extracellular matrix organization in replicative senescent cells. Analysis of Senes_MtD indicates the up
regulation of glycolysis related genes in the situation of mitochondria depletion.

## Discussion
In this report, the influence of mitochondria depletion was determined, which could result in the stimulation of glycolysis. 
Mitochondria respiration and glycolysis both utilize glucose to generate ATP, with 36 ATP generated in mitochondria respiration and 2 ATP generated in
glycolysis per glucose molecule. Mitochondria respiration is used for most cells when oxygen is sufficient.
However, [proliferating cells like fibroblasts use aerobic glycolysis even when oxygen is abundant](https://doi.org/10.1016/j.molcel.2020.12.012). 
The depletion of mitochondria may result in increasing of this phenomenon, which is said to [promote cell proliferation](https://doi.org/10.1038/ncb3124). 
Further analysis is needed to investigate the impact of mitochondria depletion on replicative senescent state.
