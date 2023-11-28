# Lab_DNA_Methylation

With climate change rapidly impacting ecosystems, phenotypic plasticity that operates on timescales quicker than evolutionary adaptation is crucial for species persistence. Epigenetics, the study of changes to gene expression that do not come from changes in the DNA sequence, is emerging as a promising avenue to understand mechanisms behind phenotypic plasticity. The addition of methyl (CH<sub>3</sub>) groups to cytosine bases adjacent to guanine bases (CpGs), or DNA methylation, is the most studied epigenetic mechanism due to increasing accessibility and affordability of sample preparation and sequencing technologies. Additionally, DNA methylation is environmentally-sensitive and potentially heritable by offspring, suggesting an important role of methylation in intra- and intergenerational plasticity.

In this lab, you'll go through several steps of a typical DNA methylation data analysis workflow using data from [Venkataraman et al. (2022)](http://dx.doi.org/10.1186/s12864-022-08781-5):

1. Use the [Integrative Genomics Viewer](https://igv.org/doc/desktop/#DownloadPage) to visualized cleaned methylation data aligned to the Pacific oyster genome and various differentially methylated loci (DML) tracks
2. Understand parameters used to identify environmentally-sensitive DML
3. Interpret [`bedtools`](https://bedtools.readthedocs.io/en/latest/index.html) output to characterize genomic location of DML
4. Complete a gene enrichment analysis with [`topGO`](https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf) and visualize results with [REVIGO](http://revigo.irb.hr)

## Data visualization with IGV

[IGV](https://igv.org/) is a GUI used to visualize various forms of genomic data. Today, we'll use it to get an understanding of what DNA methylation data "looks" like!

- Download and install IGV on your local machine. You can use the IGV web browser, but I've found it's much slower and glitchier than using it on your local machine.
- Open IGV, which opens a new session. Upload the [Pacific oyster genome](https://gannet.fish.washington.edu/panopea/Cg-roslin/cgigas_uk_roslin_v1_genomic-mito.fa) and [index file](https://gannet.fish.washington.edu/panopea/Cg-roslin/cgigas_uk_roslin_v1_genomic-mito.fa.fai) using the "Load Genome from URL" option.

<img width="568" alt="Screenshot 2023-11-24 at 11 40 24 AM" src="https://github.com/2023-environmental-bioinformatics/Lab_BLAST/assets/22335838/8adbb26e-1409-4849-8448-ef948eee761a">

- Add one methylation data file using the "Load from File" option: [sample file](https://gannet.fish.washington.edu/spartina/project-gigas-oa-meth/output/bismark-roslin/zr3616_1_R1_val_1_val_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz)

- Familiarize yourself with the different ways to visualize the data, such as navigating to a different part of the genome, zooming in and out, changing the color of track, etc. Is the genome methylated throughout, or are there concentrated pockets that are more methylated?

The R package [`methylKit`](https://github.com/al2na/methylKit) is commonly used to identify differentially methylated loci (DML). Users need to define a methylation difference. A 50% difference is commonly used to identify DML, but some studies have used other.

- Add the [BEDfiles of DML](https://gannet.fish.washington.edu/spartina/project-gigas-oa-meth/output/DML-characterization/) with various methylation thresholds. Zoom into the individual loci. Which parameter (25% difference, 50% difference, or 75% difference) is the "best"?
  - [25% difference](https://gannet.fish.washington.edu/spartina/project-gigas-oa-meth/output/DML-characterization/DML-pH-25-Cov5-Fem.csv.bed)
  - [50% difference](https://gannet.fish.washington.edu/spartina/project-gigas-oa-meth/output/DML-characterization/DML-pH-50-Cov5-Fem.csv.bed)
  - [75% difference](https://gannet.fish.washington.edu/spartina/project-gigas-oa-meth/output/DML-characterization/DML-pH-75-Cov5-Fem.csv.bed)

## Characterizing the genomic location of DML

Now that you have a list of DML, the next step is to figure out where they are in the genome!

- [Install BEDtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
  - For Mac OS users, I prefer using the `homebrew` installation method
- Download the Pacific oyster gene track and one of the other tracks.
  - [Genes](http://owl.fish.washington.edu/halfshell/genomic-databank/cgigas_uk_roslin_v1_gene.gff)
  - [Upstream flanking regions](http://owl.fish.washington.edu/halfshell/genomic-databank/cgigas_uk_roslin_v1_upstream.gff) and [downstream flanking regions](http://owl.fish.washington.edu/halfshell/genomic-databank/cgigas_uk_roslin_v1_downstream.gff):
  - [Transposable elements](http://owl.fish.washington.edu/halfshell/genomic-databank/cgigas_uk_roslin_v1_rm.te.bed): These are elements of the genome that can move around to regulate other genes
  - [Long non-coding RNA](http://owl.fish.washington.edu/halfshell/genomic-databank/cgigas_uk_roslin_v1_lncRNA.gff): When transcribed, they can regulate transcription of other genes
  - [Intergenic regions](http://owl.fish.washington.edu/halfshell/genomic-databank/cgigas_uk_roslin_v1_intergenic.bed):
- Use the [`bedtools intersect`](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) to examine the overlap between DML and the various genome features.

```
# Specify path for intersectBed
# -u: Write original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B.
# -a: BEDfile “A”. This is the DML list. Each feature in A is compared to B in search of overlaps.
# -b: BEDfile "B". This is the genome feature file
# Specify path for output file that lists overlaps
# head: View first few lines of output file
# wc -l: Count number of overlaps

{path-to-BEDtools-directory}/intersectBed \
-u \
-a {path-to-BEDfile-with-DML}.bed \
-b {path-to-genome-feature-file} \
> {path-to-output-file}.bed
!head {path-to-output-file}.bed
!wc -l {path-to-output-file}.bed
```

- If you have time, figure out what set of commands you need to create the upstream and downstream flanking region tracks.
  - Hint: look at `bedtools flank`
  - Another hint: the flanking regions shouldn't overlap with any other genome feature

## Gene enrichment analysis

The final step is to understand what biological processes are associated with genes that contain DML. This helps us interpret the function of methylation in regulating various pathways. To do this, we will use `topGO` to conduct a gene enrichment test, then visualize the results with [REVIGO](http://revigo.irb.hr). A lot of gene enrichment tools were created for gene expression analysis, but we can use them for methylation analysis too!

- Install the R package `topGO`

To install this package, start R (version "4.3") and enter:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("topGO")
```

- Download the input files for enrichment analysis
  - [gene2go file](https://gannet.fish.washington.edu/spartina/project-oyster-oa/Haws/GO-annotation/geneid2go-union1x.tab)
  - [DML, transcript IDs, and GO term annotations](https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/output/11-GOterm-annotation/DML-pH-50-Cov5-All.GeneIDs.geneOverlap.transcriptIDs.GOAnnot)
- Conduct an enrichment test for biological process GO terms with [`topGO`](https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf) manual

```
#Create "gene universe"

geneID2GO <- readMappings(file = "{path-to-gene2go-file}") #Loading the GO annotations and GeneIDs. Each line has one transcript ID and all associated GOterms
geneNames <- names(geneID2GO) #Extract names to use as gene universe
```

```
#Create gene list

allDMLGOtermsFiltered <- read.delim(file = "{path-to-DML-annotation-file}", sep = "\t", header = FALSE, col.names(“transcript”, “ID”, “car”, “start”, “stop”, “meth.diff”, “GOterm”, “name”, “GOSlim”, “cat”))
allGenes <- allDMLGOtermsFiltered$transcript #Extract transcript name
allGeneList <- factor(as.integer(geneNames %in% allGenes))
names(allGeneList) <- allGenes
str(allGeneList)
```

```
#Prepare for statistical testing

allGOdataBP <- new("topGOdata", ontology = "BP", allGenes = allGeneList,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO) #Create biological process topGO object
allGOdataBP #Get summary of object
```

Available genes have annotations, but feasible ones are linked to the GO hierarchy that topGO uses (I think...see [this link](https://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html)).

```
#Conduct statistical testing

test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher.allBP <- getSigGroups(allGOdataBP, test.stat)
resultFisher.allBP
```

```
#Extract p-values

pvalFis.allBP <- score(resultFisher.allBP) #Extract p-values
head(pvalFis.allBP)
hist(pvalFis.allBP, 50, xlab = "p-values") #Plot histogram of p-values
```

```
#Summary table with results

allRes.allBP <- GenTable(allGOdataBP, classic = resultFisher.allBP, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis.allBP)) #Create a statistical results table with statistical test results. Order by p-value (classic), and include all results (topNodes)
head(allRes.allBP)
```

```
write.csv(allRes.allBP, "{path-to-output-file}", quote = FALSE, row.names = FALSE) #Save dataframe
```

- Upload the significantly enriched GO terms and p-values ("classic" in the `topGO` output nomenclature) to [REVIGO](http://revigo.irb.hr) using default settings. You can either copy and paste this list or upload a document.

<img width="766" alt="Screenshot 2023-11-25 at 10 31 42 AM" src="https://github.com/yaaminiv/apalm-hypoxia-omics/assets/22335838/730f8f52-2025-4733-9194-32e2b1a0bdc8">

- Click on "scatterplot". This will take you to a semantic similarity plot, where each dot is a different GO term. The size of the dot corresponds to how often it was found in the dataset. Dots that are closer together have more similar biological function. The color of the dot roughly correlates to p-values (or whatever values you include with the GO terms when uploading your list to REVIGO).

<img width="825" alt="Screenshot 2023-11-25 at 10 32 16 AM" src="https://github.com/yaaminiv/apalm-hypoxia-omics/assets/22335838/ffa4f726-093f-48a8-8d08-af42fdd7fd66">

- Mess around with the settings! If you click on a dot, you can see the name of the GO term. What does a "useful" plot look like?
