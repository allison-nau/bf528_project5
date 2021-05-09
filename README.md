# bf528_project5

This repository contains the scripts used to process and analyze the data to reproduce a small portion of the analysis in Baron et al.

Baron, M., Veres, A., Wolock, S. L., Faust, A. L., Gaujoux, R., Vetere, A., Ryu, J. H., Wagner, B. K., Shen-Orr, S. S., Klein, A. M., Melton, D. A., & Yanai, I. (2016). A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. Cell systems, 3(4), 346–360.e4. https://doi.org/10.1016/j.cels.2016.08.011

Written for Boston University class BF528  
Allison Nau  
Scripts modified from scripts previously written by Professor Adam Labadorf, Mae Rose Gott, and Sheila Yee.  

# Repository Contents

### run_tophat_anau.qsub ###
Aligns reads using samtools bowtie2 boost and tophat.

Submit as qsub job to cluster:
```
dos2unix run_tophat_anau.qsub
qsub -P <project> run_tophat_anau.qsub
```


### rseqc_metrics_anau.qsub ###
Script to samtools index accepted_hits.bam, then run geneBody_coverage.py, inner_distance.py, and bam_stat on the indexed bam file of the accepted hits. 
Outputs two pdfs and a text file on the read stats.

Submit as qsub job to cluster:
```
dos2unix rseqc_metrics_anau.qsub
qsub -P <project> rseqc_metrics_anau.qsub
```


### run_cufflinks_anau.qsub ###
Finds FPKM values on a bam file of accepted hits. Uses cufflinks.

Submit as qsub job to cluster:
```
dos2unix run_cufflinks_anau.qsub
qsub -P <project> run_cufflinks_anau.qsub
```


### cufflink_hist_anau.R ###
R script to create density histograms of FPKM values.

Recommend running through Rstudio. To run through command line instead:
```
module load R/4.0.2
Rscript cufflink_hist_anau.R
```


### run_cuffdiff_anau.qsub ###
Identifies differentially expressed genes taking in several bam files. Uses cuffdiff.

Submit as qsub job to cluster:
```
dos2unix run_cuffdiff_anau.qsub
qsub -P <project> run_cuffdiff_anau.qsub
```


### analyst_anau.R ###
R script to identify which differentially expressed genes generated by cuffdiff were significantly up- or down-regulated. 
Outputs csv of up regulated genes, csv of down regulated genes, csv of 10 of the most significant differentially expressed genes, 
csvs with significant gene count for q-value<0.01 and q-value<0.05, 
histograms of log2 fold change.

Recommend running through Rstudio. To run through command line instead:
```
module load R/4.0.2
Rscript analyst_anau.R
```


Note: files were written on Windows. Recommend running dos2unix before submitting as a qsub job on a computing cluster.
