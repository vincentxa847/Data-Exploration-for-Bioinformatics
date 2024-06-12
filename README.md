# Using bulk RNA-seq dataset to compare the transcriptomes of proliferating, replicative senescent, and mitochondria-depleted senescent fibroblasts

## Result
Genes that show significant differential expression (p.adj < 0.001, absolute log2fold > 2) in three comparisons (Senes_MtD_vs_Prolif, Senes_MtD_vs_Senes and Senes_vs_Prolif) were selected separately
(with prefix Sig-) for differential gene expression analysis. 
```
#  A function is created to generate a "master" table that selects differential genes.
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
![Volcano Plot](https://github.com/vincentxa847/Data-Exploration-for-Bioinformatics/assets/118545004/6636ee2f-0cde-4f80-89b9-8a7e9f218e40)
*Volcano plot of differential expression genes in three comparsions*
