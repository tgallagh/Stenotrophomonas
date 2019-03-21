# STENO Phylogenomics, Metabolomics, Transcriptomics, Metatranscriptomics


## Growth data from in vitro experiments
### &nbsp;&nbsp; Scripts
* Main growth data in THB [/code/scripts/R](https://github.com/tgallagh/Stenotrophomonas/blob/master/code/scripts/R/Steno_GrowthCurve_fig2.R)
* Growth data in ASM [/code/scripts/R](https://github.com/tgallagh/Stenotrophomonas/blob/master/code/scripts/R/Steno_GrowthCurveASM_figS2.R)
### &nbsp;&nbsp; Data
* Growth curve data [/GrowthData](https://github.com/tgallagh/Stenotrophomonas/tree/master/GrowthData) 


## Phylogenomics
### &nbsp;&nbsp; Scripts
*  Phylogenomics by A. Chase and A. Oliver
*  contact for info about scripts for phylogenomics
*  compare accessory gene content: [/code/scripts/R](https://github.com/tgallagh/Stenotrophomonas/blob/master/code/scripts/R/accessorygenes_pathwayanalysis_fig1.R)
* look at distrubution of genes across Steno strains [/code/scripts/R](https://github.com/tgallagh/Stenotrophomonas/blob/master/code/scripts/R/gene.distrubution.phyloFigS1.R)
### &nbsp;&nbsp; Data
* 153 strain tree using 21 marker genes: [/data/processed/genome/phylogenetic/Old_153strains](https://github.com/tgallagh/Stenotrophomonas/tree/master/data/processed/genome/phylogenetic/Old_153strains)
* core genome tree (74 strains, 95% AAI clusters): [/data/processed/genome/phylogenetic/newtree_alex](https://github.com/tgallagh/Stenotrophomonas/tree/master/data/processed/genome/phylogenetic/newtree_alex)


## Metabolomics
### &nbsp;&nbsp; Scripts
*  PCOA and volcano plots of metabolomics data in R: [/code/scripts/R](https://github.com/tgallagh/Stenotrophomonas/blob/master/code/scripts/R/metabolites_pcoa_volcanoplot_fig4_figs3.R) 
### &nbsp;&nbsp; Data
* metabolomics data from WCMC: [/data/metabolomics](https://github.com/tgallagh/Stenotrophomonas/tree/master/data/metabolomics)

## Transcriptomics
### &nbsp;&nbsp; Scripts
* Pipeline for data from raw to alignments and HT-Seq Count hits [/code/scripts/Transcriptome](https://github.com/tgallagh/Stenotrophomonas/tree/master/code/scripts/Transcriptome)
* transcriptome analyses in R, including DGE [/code/scripts/R/trancriptome_R_stuff](https://github.com/tgallagh/Stenotrophomonas/tree/master/code/scripts/R/transcriptome_R_stuff)
### &nbsp;&nbsp; Data
* raw data = https://www.ncbi.nlm.nih.gov/bioproject/?term=GSE121347 
* BAM alignments to FLR19 genome: [/data/processed/bowtie2](https://github.com/tgallagh/Stenotrophomonas/tree/master/data/processed/bowtie2)
* HT-Seq output: [/data/processed/HTSeq](https://github.com/tgallagh/Stenotrophomonas/blob/master/data/processed/HTSeq/steno_htseq.xls) 

## Metranscriptomics
### &nbsp;&nbsp; Scripts
* quality filter data in linux environment: [/code/scripts/Metatranscriptome](https://github.com/tgallagh/Stenotrophomonas/blob/master/code/scripts/Metatranscriptome/qualityfilter_Whiteley.sh)
* align to non-Steno CF database: [/code/scripts/Metatranscriptome](https://github.com/tgallagh/Stenotrophomonas/blob/master/code/scripts/Metatranscriptome/bowtie1_CFdatabase.sh)
* align left-overs to Steno pan-genome:[/code/scripts/Metatranscriptome](https://github.com/tgallagh/Stenotrophomonas/blob/master/code/scripts/Metatranscriptome/bowtie2_invitro_pangenome.sh)
* look at gene alignments:   [/code/scripts/Metatranscriptome](https://github.com/tgallagh/Stenotrophomonas/blob/master/code/scripts/Metatranscriptome/pangenomealignments.R)
* look at gene hits in sputum E and F: [/code/scripts/Metatranscriptome](https://github.com/tgallagh/Stenotrophomonas/blob/master/code/scripts/Metatranscriptome/calculate_pangenome_hits.R)

### &nbsp;&nbsp; Data
* meta transcriptomes from Marvin Whiteley's Lab: https://www.pnas.org/content/115/22/E5125.short
* alignment to Steno pan-genome: [/data/processed/CFdata/Whiteley/new_aligned](https://github.com/tgallagh/Stenotrophomonas/tree/master/data/processed/CFdata/Whiteley/new_aligned)
 

