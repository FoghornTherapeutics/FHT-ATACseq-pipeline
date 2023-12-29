# Overview

Chromatin DNA replication and gene transcription occur when the higher-order structure of DNA becomes loose, which is what we call open chromatin. Identifying these open regions is crucial because they often hold the keys to how genes are controlled and understande the genome's regulatory landscape.

Assay for Transposase-Accessible Chromatin using sequencing (ATACseq) is a powerful technique that has revolutionized our understanding of the genome's regulatory landscape. It uses sequencing adapters, like cutting enzyme, to fragment the open chromatin and, after purification, the library can be amplified by PCR. Then it can be analyzed by NGS to obtain fastqs.


# Data

The standard pipeline is ran on publicly available data from paper "[Chromatin accessibility underlies synthetic lethality of SWI/SNF subunits in ARID1A-mutant cancers](https://elifesciences.org/articles/30506#content)" looking for potential PD markers as well as what an ATAC-seq profile looks like. This paper has ATACseq results of ARID1A-/- cancers with ARID1B KD. 

data from GEO series:  [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101975)

Biological context (N=2):  TOV21G, HCT116
wild type and modified with stable ARID1A KO
Perturbagens (N=1):  shRNA KD of ARID1B
Doses (N/A):  just the shRNA no relevant dose
Negative control (N=1):  wild type / untreated
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



## Alignment

## Analysis




# Results


1) QC check:
   
   A) Insert size: classical ATAC-seq pattern with nucleosome free peak, dinucleosome peak, trinucleosome peak, etc.
   
   B) Duplication statistics: ~40% of the reads uniquely mapped to human genome
   
   C) FRiP scores: about 30% of the total peaks are found within the peaks


   D) nIDR: Replicates of the same group show consistency.

   
   E) PCA (Principal Component Analysis) plot: clear separation between cell lines with the Principal Component 1 and 2. PC2 clusters HCT116 by treatment.


   F) Sample to sample correlation: HCT116 shows more overlap between samples than TOV21G.
   
   
2) Peak distribution:

   * ARID1A KO in WT HCT116 cells dramatically altered overall chromatin accessibility resulting in thousands of increased and decreased sites.
   * ARID1B KD in WT HCT116 cells had little effect on chromatin accessibility.
   * In contrast, ARID1B KD in ARID1A-/- HCT116 cells resulted in hundreds of changed sites, primarily at regions where accessibility was lost.
   * ARID1A-mutant TOV21G cell line infected with shRNAs to ARID1B (ARID1B KD) showed little effect of chromatin accessibility.

     
3) Genome location: Majority of loss of accessibility occur at distal intergenic regions. Decreased sites are enriched at intronic regions regions while increased sites are enriched at promoters.

  
     
4) Motif analysis:
   *
   *
   * Accessibility at AP-1 motif sites was also decreased in TOV21G cells following ARID1B knockdown
     
5) Footprinting





