# Background

DNA replication and gene transcription occur when the higher-order structure of DNA - heterochromatin - is unpacked creating regions of open chromatin. Identifying these open regions is crucial because they often hold the keys to how genes are controlled and aid understanding the genome's regulatory landscape.

Assay for Transposase-Accessible Chromatin using sequencing (ATAC-seq) is a powerful technique that has revolutionized our understanding of the genome's regulatory landscape. It uses sequencing adapters and a cutting enzyme to fragment and tag DNA in the open chromatin regions, after purification, the library can be prepared (including PCR amplification) and sequenced using next-generation sequencing techniques.

# Overview

Here we present a computational pipeline and related techniques for analyzing ATAC-seq data.  The ultimate goal of this pipeline is identify genomic regions whose chromatin accessibility changes when a biological system is perturbed (using a compound or a genetic alteration i.e. CRISPR).  Here we will walk through an example of the outputs of this pipeline using publicly available data to illustrate the main outputs and how we use them.


# Data

The standard pipeline was run on publicly available data from paper "[Chromatin accessibility underlies synthetic lethality of SWI/SNF subunits in ARID1A-mutant cancers](https://elifesciences.org/articles/30506#content)" looking for potential PD markers as well as what an ATAC-seq profile looks like. This paper has ATACseq results of ARID1A-/- cancer cell lines (native or CRISPR knockout) with ARID1B knockdown. 

**Data from GEO series**:  [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101975)

Overview of experiment:
* Biological context (N=2):  TOV21G, HCT116
* wild type and modified with stable ARID1A knockout
* Perturbagens (N=1):  shRNA knockdown of ARID1B
* Negative control (N=1):  wild type / untreated
* Replicates:  N=2



### HCT116 (ACH-000971)
* WT: SRR5876158 & SRR5876159
* ARID1B knockdown: SRR5876160 & SRR5876161
* ARID1A knockout: SRR5876162 & SRR5876163
* ARID1A knockout ARID1B knockdown:SRR5876164 & SRR5876165

### TOV21G (ACH-000885)
* WT: SRR5876661 & SRR5876662
* ARID1B knockdown: SRR5876663 & SRR5876664


# Description of the pipeline

The overall ATACseq architecture is split between two parts: the alignment and the downstream analysis.  Alignment focuses on alignment of reads  the genome and identification of regions of higher read counts - "peaks".  Downstream analysis focuses on comparisons across/between samples.

## Alignment

The alignment includes steps to remove adapters, align to the genome, remove duplicates and filter to reads less than 120 bp in length (to retain fragments that represent open regions). We obtain two major output files that are the BAM files (aligned reads) and BED/narrowPeak files (identifying regions of open chromatin, "peaks").

<img src="images/data_flow/alignment_data_flow.JPG" alt="image" style="width:600px;height:auto;">

You can find a more detailed version of the data flow diagram of the aligmnent [here](https://github.com/FoghornTherapeutics/FHT-ATACseq-pipeline/blob/main/code/alignment/README.md). 

## Analysis

As mentioned one goal of the ATACseq pipeline is to identifiy regions where the chromatin accessibility is changing.  An example of a region like this is shown below in a genome browser view.  The upper three traces correspond to DMSO (vehicle/untreated) samples, and the lower three traces correspond to samples treated with a compound.  The coverage of reads in the upper three traces is consistently higher than in the lower three, indicating that chromatin accessibility has decreased as a result of this treatment.

![](images/data_flow/IGView.JPG)

It is useful to contrast RNA-seq with ATAC-seq to illustrate a challenge of analyzing ATAC-seq data.  In RNA-seq, generally, one uses the predefined exons and then counts the reads that overlap with exons to generate gene expression read counts.  For ATAC-seq there are no pre-defined regions analogous to exons that can be used for generating a count matrix, so an early step of the analysis is to identify these regions (aka peaks).  These peaks are identified for individual samples, and then need to be reconciled across samples in order to allow subsequent comparisons / aggregate analysis across samples.

Once we have identified these common regions aka peaks we can generate a count matrix, and then we can look at PCA (Principal Component Analysis) and at sample-to-sample correlations.  One important aspect of these analyses is for quality control to look for individual samples that are clear outliers from the other replicates.

![](images/data_flow/QC_data_flow.JPG)


Before we dive in the second part, we define some terminology. We call each sample a replicate, and replicates with the same conditions (treatment, cell lines, time, dose, etc.) are a group.  The comparison of two groups (i.e., test vs. negative control; treated vs. untreated etc.) is called a contrast.  Regions of interest are also called "peaks" for short and are genomic locations where open chromatin has been observed for 1 or more samples.

For a given contrast we define the peaks as the union of the peaks for the groups within the contrast.  The peaks in the groups are determined using nIDR, see below for details.  Using the peaks for a contrast we calculate the peak area as the number of reads in each peak for each group, and then compute the differential peak area for each contrast to get a logFC and a pvalue for each peak. The peaks are categorized by the observed change of chromatin accssibility: descrease (down), increase (up) or unchanged. The three groups of peaks are then analyzed to determine genomic locations (promoter vs enhancer), presence of motifs, tornado plots, and footprinting.


![](images/data_flow/DownstreamAnalysis.JPG)



# Results


1) QC check:

The first part of the analysis starts with some QC measures:

   A) **Insert size**

   The insert size plot is an histogram of the distribution of the length of DNA betwen the two sequencing reads (paired-end sequencing). The x-axis represents the insert size (length of DNA between the paired reads). This size can vary depending on how the DNA was fragmented and prepared for sequencing. The y-axis shows the frequency or count of read pairs for each insert size. Each line corresponds to a unique sample. <br/> 
   
Generally, a well-conducted ATAC-seq experiment is expected to produce a plot showing the distribution of fragment sizes with a pattern of decreasing, periodic peaks. These peaks correspond to regions without nucleosomes, also called Nucleosome Free Regions (NFR), which are typically less than 120 base pairs, and regions with mono-, di-, and tri-nucleosomes, approximately at 200, 400, and 600 base pairs, respectively. Deviations from the expected pattern, like multiple peaks or a very broad distribution, can indicate issues such as DNA fragmentation problems or contamination. In this example, it is the classical ATAC-seq pattern with nucleosome free peak, dinucleosome peak, trinucleosome peak, etc.

   NB: The oscillation of the insert size is due to the DNA helix shape (double stranded) including a major and minor group wrapped around each other. The enzyme has a preference for one of the groups.

   <img src="/images/output_results/picard_insert_size.png" alt="image" style="width:700px;height:auto;">
   

   
   B) **Duplication statistics** 

   Duplicate reads are identical or nearly identical sequences that appear multiple times in the sequencing data. They can arise due to biological reasons (e.g., highly repetitive regions in the genome) or technical reasons (e.g., PCR amplification during library preparation). They can skew the data analysis and lead to inaccurate representation of gene or variant frequencies. 

   Once removed, we compute the deduplication statistics to provide insights into the number and proportion of duplicate reads found in the sequencing data. Common metrics include the total number of reads, the number of unique reads, the number of duplicate reads, and the percentage of reads that are duplicates. high levels of duplication may indicate over-amplification during PCR or other library preparation issues. Conversely, very low duplication levels might suggest a diverse library or under-sampling. In this example, about 40% of the reads uniquely mapped to human genome which is a typical number for ATAC seq experiments.
      	
   <img src="/images/output_results/picard_deduplication.png" alt="image" style="width:700px;height:auto;">
   
   C) **Fraction of Reads in Peaks (FRiP) scores**

The FRiP scores are a quantitative measure of how effectively the DNA fragments of interest have been enriched and captured in the experiment. For well-characterized transcription factors or histone modifications with strong and specific binding, FRIP scores are often expected to be higher. A score in the range of 1% to 10% is commonly seen as acceptable, with some experiments even reaching 20% or higher. In this case, the FRiP scores are around 3%, i.e. 3% of the total peaks are found within the peaks. It is low compared to some other results we could see in the past but still acceptable.

   <img src="/images/output_results/FRiP_Table.JPG" alt="image" style="width:550px;height:auto;">



   D) **n Irreproducible Discovery Rate (nIDR)**
   
   IDR is used to measure the consistency of results (such as the identification of transcription factor binding sites, histone marks, or gene expression levels) across different experimental replicates. The standard method involves comparing the rank of results across different replicates to estimate the proportion of findings that are reproducible (consistent across replicates) versus those that are irreproducible (inconsistent or likely to be noise). The nIDR method is technique we have developed from the orignial and detailed [here](). 

One of the nIDR outputs is the Empirical Cumulative Distribution Function (ECDF) plots of the consistency across replicates of all peaks found in a group of replicates. The x-axis is the "min percent rank" which indicates the consistency of a peak across the replicates, a higher value corresponds to a higher consistency of peaks accross replicates. The y-axis is the fraction peaks that have at least that value. The red curve is plotting the actual data and the blue curve is the simulated null distribution. 

Looking at the first plot, there are two vertical lines on the left side, the highest point of the red part (real data) of the vertical line is at y=0.75, x=0.22.  This means that 75% of the real data has a min percent rank (consistency score) of 0.21 or less. At the higher consistency scores around x > 0.6, the red curve of real data is below the blue curve of the null data, indicating that real data has higher consistency scores than the null.  This indicates the replicates have peaks that are consistent. The green dashed line indicates a p-value of 0.1 based on the blue null curve and determines the consistency score threshold to use for keeping the real peaks. In this case, the green dashed line indicates the threshold where 90% of the null peaks have a consistency score below 0.65. Therefore this sets a threshold for choosing peaks with a p value < 0.1. 0.65 is considered to be high, with a null distribution being above the actual one, replicates of the same group show consistency.


   <img src="/images/output_results/HCT116_nIDR.JPG" alt="image" style="width:700px;height:auto;">
   <img src="/images/output_results/TOV21G_nIDR.JPG" alt="image" style="width:700px;height:auto;">


   
   E) **Principal Component Analysis (PCA) plot**
   
   PCA, or Principal Component Analysis, is a statistical technique used to simplify the complexity in high-dimensional data while retaining trends and patterns in a two-dimensional space.

The PCA shows a clear separation between the cell lines with the samples from TOV21G on the upper left and HCT11 cells in the lower right. <br/> 


   <img src="/images/output_results/PCA_all_samples.JPG" alt="image" style="width:760px;height:auto;">
        
Among HCT116 cells only, we see a clear separation of the double treatment (ARID1A knockout and ARID1B knockdown) with the blue points being in the far left. ARID1A knockout samples are also clustering together in the middle in purple points. There are no clear separation between the HCT116 WT cells and ARID1B knockdown. It could be an indication that the double treatment has the strongest effect on WT cells, then ARID1A knockout alone has less strong effect and ARID1B knockdown does not have a strong effect in chromatin accessibility in WT HCT116 cells. 
   
   <img src="/images/output_results/HCT116_PCA.JPG" alt="image" style="width:760px;height:auto;">

There is no clear interpretation for the PCA with TOV21G cells. We also predict a weak effect from ARID1B knockdown in ARID1A-mutant TOV21G cells.

   <img src="/images/output_results/TOV21G_PCA.JPG" alt="image" style="width:700px;height:auto;">
      
  
   
   F) **Sample to sample correlation**

  Another QC measure is to compute the pearson correlation between samples of their feature counts.

  The first visualization is to look at the heatmap of the sample-to-sample correlation, clustered by similarities. Each column and row is a sample.
  
  . Just like in the PCA plot, samples cluste

HCT116 shows higher correlation
say it is consistent w/ PCA
explain box plot of sample-sample correlation - these are the same correlation values from the heatmap but plotted for each sample against the other samples of the same cell line   
   
   * HCT116 shows more overlap between samples than TOV21G.
   * Again, WT HCT116 cells and ARID1B KD HCT116 cells cluster together and ARID1A-mutant TOV21G cells and ARID1B KD TOV21G cluster together which predicts that these two contrasts will present weak results.

     <img src="/images/output_results/Heatmap_sample_to_sample_corr.JPG" alt="image" style="width:550px;height:auto;">
     
        ![](images/output_results/boxplot_sample_to_sample_corr.JPG)
   
   Another example with one of the samples being an outlier is described [here](https://github.com/FoghornTherapeutics/FHT-ATACseq-pipeline/blob/main/QC_example_with_outlier.md).



   DEfine the contrasts!!!
   
3) Peak distribution:

As a reminder, we divide peaks into three categories: Up (p-value < 0.05 &  logFC > 0.5), Down (p-value < 0.05 &  logFC < - 0.5) and Unchanged (p-value < 0.05 or |logFC| < 0.5).

   * ARID1A KO in WT HCT116 cells dramatically altered overall chromatin accessibility resulting in thousands of increased and decreased sites.
   * In contrast, ARID1B KD in WT HCT116 cells had modest effect on chromatin accessibility.
   * When combining ARID1B KD and ARID1A-/-, HCT116 cells resulted in thousands of additional sites that lost accessibility.
   * ARID1A-mutant TOV21G cell line infected with shRNAs to ARID1B (ARID1B KD) showed little effect of chromatin accessibility.
   ![](images/output_results/peak_dist_boxplot.JPG)

   * The upset plot shows overlap between ARID1A KO and ARID1B KD in ARID1A-/- HCT116 cells:
     
   ![](images/output_results/Upset_plot_down.JPG)
   ![](images/output_results/Upset_plot_up.JPG)


     
4) Genome location:


Majority of loss of accessibility occur at distal intergenic regions. Decreased sites are enriched at intronic regions regions while increased sites are enriched at promoters.

   ![](images/output_results/Genomic_location.JPG)

  
     
5) Motif analysis: 
   

For the Motif analysis, we use [HOMER](http://homer.ucsd.edu/homer/motif/). It looks at the peaks in the Target Sequence. for example the peaks loosing chromatin accessibility and compare them to the Total Background Dequences, i.e. the peaks with unchanged chromatin accessibility. It then compares sequences that have known motifs with the reference that are in the Target Sequence.

   * Sites losing chromatin accessibility are strongly enriched in the AP-1 family in ARID1A-/- HCT116 cells and relatively highly enriched in ARID1A-mutant TOV21G cell lines over the total number of peaks.
   * However, motifs in ARID1B KD in WT HCT116 cells have low p-values and are mostly in the TEAD family.


  
  Motifs losing peak accessibility: 
  
   ![](images/output_results/Motifs_down.JPG)
  
  Motifs gaining peak accessibility: 
  
   ![](images/output_results/Motifs_up.JPG)

Heatmap of motifs accross all contrasts, clustered by motifs:

   <img src="/images/output_results/motif_heatmap.JPG" alt="image" width="20%" height="auto">
     

Zoom in of certain clusters in the heatmap:

describe the axis/ color (blue is loss and red is gain) and intensity is the log(log(pValue of Motifs))


   <img src="/images/output_results/motif_heatmap_zoom.JPG" alt="image" style="width:550px;height:auto;">
     



   + ADD SCATTER PLOT OF MOTIFS
     
6) Tornado plot: It is an overall look at the whole data at once looking at a very condensed view, where each row is a peak and the intensity of the color represents the read count.

   We first order peaks within contrast for WT or DMSO depending on the read counts. Then we concatenate peak accessibility going down/unchanched/up. We keep the same order of the peaks for the treatment. We can observe on the top, peaks losing accessibility and show less intensity with the compound/treatment i.e. there are less reads). The opposite is true on the bottom with more intensity in the compound where peaks gain chromatin accessibility.

   peaks losing their color itensity ant in the bottom gaining 




   <img src="/images/output_results/TP_1.JPG" alt="image" style="width:900px;height:auto;">
   <img src="/images/output_results/TP_1.JPG" alt="image" style="width:900px;height:auto;">
     

     



7) Footprinting

[Rgt-hint](https://reg-gen.readthedocs.io/en/latest/hint/tutorial-dendritic-cell.html) generates new bed files that consider peak regions for footprinting. It then finds motifs overlapping with predicted footprints and  generates average ATAC-seq profiles around binding sites of particular TFs. 

We can see in this example that when combining ARID1A knockout and ARID1B knockdown, we have a higher logFC between the two groups (WT vs compound) than we the two other first contrasts. The volcano plots highlights a lot of TF motifs that are downregulated. The heatmap (ordered by desending absolute value of logFC and filtered for only significant logFC) shows a clear difference in footprint scores between WT HCTT16 and ARID1A knockout and ARID1B knockdown HCT116 cells.
 <img src="/images/output_results/footprint.JPG" alt="image" style="width:1000px;height:auto;">

Rgt-hint also outputs profile plots. The x-axis id the base pair +/- 100bp either side of the FOSL1:JUNB motif footprint. The y-axis is the coverage of the BAM file reads. It correlates with the open chromatin and assumes to be TF activity. 
 
<img src="/images/output_results/footprint_profile.JPG" alt="image" style="width:500px;height:auto;">







