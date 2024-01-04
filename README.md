# Overview

Chromatin DNA replication and gene transcription occur when the higher-order structure of DNA becomes loose, which is what we call open chromatin. Identifying these open regions is crucial because they often hold the keys to how genes are controlled and understande the genome's regulatory landscape.

Assay for Transposase-Accessible Chromatin using sequencing (ATACseq) is a powerful technique that has revolutionized our understanding of the genome's regulatory landscape. It uses sequencing adapters, like cutting enzyme, to fragment the open chromatin and, after purification, the library can be amplified by PCR. Then it can be analyzed by NGS to obtain fastqs.


# Data

The standard pipeline is ran on publicly available data from paper "[Chromatin accessibility underlies synthetic lethality of SWI/SNF subunits in ARID1A-mutant cancers](https://elifesciences.org/articles/30506#content)" looking for potential PD markers as well as what an ATAC-seq profile looks like. This paper has ATACseq results of ARID1A-/- cancers with ARID1B KD. 

**Data from GEO series**:  [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101975) <br/>
Biological context (N=2):  TOV21G, HCT116 <br/>
wild type and modified with stable ARID1A KO <br/>
Perturbagens (N=1):  shRNA KD of ARID1B <br/>
Doses (N/A):  just the shRNA no relevant dose <br/>
Negative control (N=1):  wild type / untreated <br/>
Replicates:  N=2



### HCT116 (ACH-000971)
* WT: SRR5876158 & SRR5876159
* ARID1B knockdown: SRR5876160 & SRR5876161
* ARID1A knockout: SRR5876162 & SRR5876163
* ARID1A knockout ARID1B knockdown:SRR5876164 & SRR5876165

### TOV21G (ACH-000885)
* WT: SRR5876661 & SRR5876662
* ARID1B knockdown: SRR5876663 & SRR5876664


# Data Flow

The overall ATACseq architecture is split between two parts: the alignment and the downstream analysis.

## Alignment

The alignment includes several steps to remove adapters with a cutting enzyme, align to the genome, remove duplicates and filter reads, less than 120 to retain fragments that are not wrapped around nucleosomes, only keep reads in open chromatin. We obtain two major output files that are the BAM files (read counts) and narrowPeak (denser regions).

<img src="images/data_flow/alignment_data_flow.JPG" alt="image" style="width:600px;height:auto;">

You can find a more detailed version of the data flow diagram of the aligmnent [here](https://github.com/FoghornTherapeutics/FHT-ATACseq-pipeline/blob/main/code/alignment/README.md). 


The goal of the ATACseq pipeline is to identifiy regions where the chromatin accessibility is changing, just like in the one below. We cansee that the overall number of peaks is reduced when comparing DMSO with "treatment". 


![](images/data_flow/IGView.JPG)


## Analysis


The first part of the downstream analysis is to look at QC measures. We look at the overall samples in the experiment to check for any outliers. The challenge of ATAC is that we cannot compare, like RNA, gene expression across exons, but we need to find a set of peaks that is common to all samples. Basically, we look for common location to look at across the genome. With that new set of peaks, we look at PCA (Principal Component Analysis) and at sample-to-sample correlations.

![](images/data_flow/QC_data_flow.JPG)


Before we dive in the second part, we define some terminology. We call each sample a replicate. We define all replicates with the same conditions (treatment, cell lines, time, dose, etc.) as a group. And the comparison of two groups (i.e., test vs. negative control; treated vs. untreated etc.) is called a contrast.

This part of the downstream analysis compares reads at regions that are meaningful but also in common across samples of the same group, to be compared accross contrasts. Once this new set of peaks is defined, we quantify the number of reads in peaks and compute the differential peak area in each contrast to get a logFC and a pvalue by peaks. The peaks are split by change of chroamtin accssibility: descrease (down), increase (up) or unchanged. The three groups of peaks are used for the standard outputs of the downstream analysis (Genomic locations, motif, tornado plots, footprint).


![](images/data_flow/DownstreamAnalysis.JPG)



# Results


1) QC check:
   
   A) Insert size: classical ATAC-seq pattern with nucleosome free peak, dinucleosome peak, trinucleosome peak, etc.
   NB: The oscillation of the insert size is due to the DNA helix shape (double stranded) including a major and minor group wrapped around each other. The enzyme has a preference for one of the groups.

   <img src="/images/output_results/picard_insert_size.png" alt="image" style="width:700px;height:auto;">
   

   
   B) Duplication statistics: ~40% of the reads uniquely mapped to human genome
      	
   <img src="/images/output_results/picard_deduplication.png" alt="image" style="width:700px;height:auto;">
   
   C) FRiP scores: about 30% of the total peaks are found within the peaks


   <img src="/images/output_results/FRiP_Table.JPG" alt="image" style="width:550px;height:auto;">



   D) nIDR: Replicates of the same group show consistency.
   <img src="/images/output_results/HCT116_nIDR.JPG" alt="image" style="width:700px;height:auto;">
   <img src="/images/output_results/TOV21G_nIDR.JPG" alt="image" style="width:700px;height:auto;">


   
   E) PCA plot: clear separation between cell lines with the Principal Component 1 and 2. PC2 clusters HCT116 by treatment.
   *   WT HCT116 cells and ARID1B KD HCT116 cells cluster together. It could be an indication that ARID1B KD does not have a strong effect in chromatin accessibility in WT HCT116 cells.
   *   We predict similar results between ARID1A-mutant TOV21G cells and ARID1B KD TOV21G cells that cluster all together.
     

        <img src="/images/output_results/PCA_all_samples.JPG" alt="image" style="width:760px;height:auto;">
        <img src="/images/output_results/HCT116_PCA.JPG" alt="image" style="width:760px;height:auto;">
        <img src="/images/output_results/TOV21G_PCA.JPG" alt="image" style="width:700px;height:auto;">
      
  
   
   F) Sample to sample correlation:
   * HCT116 shows more overlap between samples than TOV21G.
   * Again, WT HCT116 cells and ARID1B KD HCT116 cells cluster together and ARID1A-mutant TOV21G cells and ARID1B KD TOV21G cluster together which predicts that these two contrasts will present weak results.

     <img src="/images/output_results/Heatmap_sample_to_sample_corr.JPG" alt="image" style="width:550px;height:auto;">
     
        ![](images/output_results/boxplot_sample_to_sample_corr.JPG)
   
   Another example with one of the samples being an outlier is described [here](https://github.com/FoghornTherapeutics/FHT-ATACseq-pipeline/blob/main/QC_example_with_outlier.md).
   
3) Peak distribution:

As a reminder, we divide peaks into three categories: Up (p-value < 0.05 &  logFC > 0.5), Down (p-value < 0.05 &  logFC < - 0.5) and Unchanged (p-value < 0.05 or |logFC| < 0.5).

   * ARID1A KO in WT HCT116 cells dramatically altered overall chromatin accessibility resulting in thousands of increased and decreased sites.
   * ARID1B KD in WT HCT116 cells had little effect on chromatin accessibility.
   * In contrast, ARID1B KD in ARID1A-/- HCT116 cells resulted in hundreds of changed sites, primarily at regions where accessibility was lost.
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


   <img src="/images/output_results/motif_heatmap_zoom.JPG" alt="image" style="width:550px;height:auto;">
     



   + ADD SCATTER PLOT OF MOTIFS
     
6) Tornado plot: It is an overall look at the whole data at once looking at a very condensed view, where each row is a peak and the intensity of the color represents the read count.

   We first order peaks within contrast for WT or DMSO depending on the read counts. Then we concatenate peak accessibility going down/unchancged/up. We keep the same order of the peaks for the treatment. We can observe on the top, peaks losing accessibility and show less intensity with the compound/treatment i.e. there are less reads). The opposite is true on the bottom with more intensity in the compound where peaks gain chromatin accessibility.

   peaks losing their color itensity ant in the bottom gaining 




   <img src="/images/output_results/TP_1.JPG" alt="image" style="width:900px;height:auto;">
   <img src="/images/output_results/TP_1.JPG" alt="image" style="width:900px;height:auto;">
     

     



7) Footprinting

[Rgt-hint](https://reg-gen.readthedocs.io/en/latest/hint/tutorial-dendritic-cell.html) generates new bed files that consider peak regions for footprinting. It then finds motifs overlapping with predicted footprints and  generates average ATAC-seq profiles around binding sites of particular TFs. 

We can see in this example that when combining ARID1A knockout and ARID1B knockdown, we have a higher logFC between the two groups (WT vs compound) than we the two other first contrasts. The volcano plots highlights a lot of TF motifs that are downregulated. The heatmap (ordered by desending absolute value of logFC and filtered for only significant logFC) shows a clear difference in footprint scores between WT HCTT16 and ARID1A knockout and ARID1B knockdown HCT116 cells.
 <img src="/images/output_results/footprint.JPG" alt="image" style="width:1000px;height:auto;">

Rgt-hint also outputs profile plots. The x-axis id the base pair +/- 100bp either side of the FOSL1:JUNB motif footprint. The y-axis is the coverage of the BAM file reads. It correlates with the open chromatin and assumes to be TF activity. 
 
<img src="/images/output_results/footprint_profile.JPG" alt="image" style="width:500px;height:auto;">







