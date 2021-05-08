# bf528_project5

This repository contains the scripts used to process and analyze the data to reproduce a small portion of the analysis in Baron et al.

Baron, M., Veres, A., Wolock, S. L., Faust, A. L., Gaujoux, R., Vetere, A., Ryu, J. H., Wagner, B. K., Shen-Orr, S. S., Klein, A. M., Melton, D. A., & Yanai, I. (2016). A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. Cell systems, 3(4), 346–360.e4. https://doi.org/10.1016/j.cels.2016.08.011

Written for Boston University class BF528
Allison Nau
Scripts modified from scripts previously written by Professor Adam Labadorf, Mae Rose Gott, and Sheila Yee.

# Repository Contents

#### run_tophat_anau.qsub ####
Aligns reads using samtools bowtie2 boost and tophat.

#### rseqc_metrics_anau.qsub ####
Script to samtools index accepted_hits.bam, then run geneBody_coverage.py, inner_distance.py, and bam_stat on the indexed bam file of the accepted hits. 
Outputs two pdfs and a text file on the read stats.

#### run_cufflinks_anau.qsub ####
Finds FPKM values on a bam file of accepted hits. Uses cufflinks.

#### cufflink_hist_anau.R ####
R script to create density histograms of FPKM values.

#### run_cuffdiff_anau.qsub ####
Identifies differentially expressed genes taking in several bam files. Uses cuffdiff.

#### analyst_anau.R ####
R script to identify which differentially expressed genes generated by cuffdiff were significantly up- or down-regulated.
