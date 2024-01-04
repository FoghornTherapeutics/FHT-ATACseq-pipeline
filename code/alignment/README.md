<img src="../../diagrams/ATAC_alignment_diagram.drawio.png" alt="image" style="width:750px;height:auto;">


**trim_fastqs**: remove adapter from sequence coming from ATAC-seq process <br/>

**compute sh512sum**: QC measure <br/>

**run fastqc**: quality control checks on raw sequence data <br/> 

**bowtie_align**: index the genome to the hg reference into reads <br/>
**filter_sort**: sort by read names <br/>
**select_top_alignment**: select reads with the least amount of mismatch of bp between our sequence and the hg reference <br/>
**filter_sort_top_alignment**: sort reads by chromosomal coordinates <br/>

**mark_duplicates**: mark duplicates <br/>
**remove_duplicates**: remove duplicates <br/>

**QC_fragment_size**: remove reads that are longer that a nucleotide  <br/>
**flasgstat**: report number of reads that are alignment with the genome  <br/>
**mutliQC**: report flagstat into html report  <br/>

**create_nucleosome_free**: additional filter on reads depending on the size of the insert (~120 to obtain more meaningful reads). Removes reads that are too long and that have insert right before and after nucleosomes.  <br/>

**small_final_sort**: sort by read names  <br/>
**index_small_bam**: sort reads by chromosomal coordinates and index BAM files  <br/>
**Convert_to_BEDPE_TA**: get PE (paired end) only focusing on the start and end, add a shift (2 bed columns)  <br/>
**call_peaks**: filter denser peaks  <br/>
**filter_peaks**: remove chromosome names that are chr1:22 or chrX and chrY  <br/>

**big_wig_norm_using_None**: compressed version of BAM file that is representing how dense the signal is in IGV Browser Portal  <br/>
