# title: Project 2 analyst
# Modified from Sheila Yee's script
# Allison Nau

# Scientific notations options:
options(scipen=1)

# load differential expression file analysis
diff_exp <- read.table("/projectnb2/bf528/users/dachshund/project5_anau/cuffdiff_out/gene_exp.diff",
                             header = TRUE)
# Old version:
#diff_exp <- read.table("/projectnb/bf528/users/dachshund/project_2/project-2-project-2-dachsund/programmer/cuffdiff_out/gene_exp.diff",
#                             header = TRUE)

# sort the above data table so that the smallest q_values are at the top 
# Could also consider Fold change and q-value
diff_exp <- diff_exp[order(diff_exp$q_value, diff_exp$p_value, -abs(diff_exp$log2.fold_change.)),]

# create a table with top ten differentially expressed genes, with their names, 
# FPKM values, log fold change, p-value, and q-value
# Disregard inf fold change
diff_exp_filter <- diff_exp[abs(diff_exp$log2.fold_change.)!=Inf,]
diff_exp_top10 <- diff_exp_filter[1:10, c("gene", "value_1", "value_2", "log2.fold_change.", "p_value", "q_value")]
write.csv(diff_exp_top10, file="diff_genes_top10.csv")

# create a new data frame that contains only the genes where the last column named significant is equal to yes
diff_exp_sig <- subset(diff_exp, significant=="yes")

# Create and Save Histograms
png(filename="Fold Change histograms.png", width=800, height=400)
par(mfrow=c(1,2))
# produce a histogram of the log2.foldchange column for all genes
hist(diff_exp$log2.fold_change, breaks=70, main="A) Fold Change", 
     ylim=c(0, 20000), xlim=c(-10, 10),
     xlab="Fold Change (log2)", col=c("lightblue"))
# create a second histogram of the log2 fold change values only for significant genes
hist(diff_exp_sig$log2.fold_change., breaks=70, 
     main = "B) Fold Change (Significant)", 
     ylim = c(0, 500), xlim = c(-10, 10),
     xlab = "Fold Change (log2)", col = c("lightblue"))
dev.off()

# create separate dataframes for only significant genes that had positive and negative log fold change
diff_exp_pos <-subset(diff_exp_sig, log2.fold_change. > 0 )
diff_exp_neg <-subset(diff_exp_sig, log2.fold_change. < 0 )

# Get and save the significant count:
sign_count_pos <- length(diff_exp_pos$gene_id)
sign_count_pos
sign_count_neg <- length(diff_exp_neg$gene_id)
sign_count_neg
sign_count <- data.frame(sign_count_neg, sign_count_pos)
write.csv(sign_count, file="sign_gene_count.csv")

# write the up and down regulated gene names to separate files
write.csv(diff_exp_pos$gene, "up_regulated_genes.csv")
write.csv(diff_exp_neg$gene, "down_regulated_genes.csv")

# Determine how many equally significant genes there were:
most_sign <- diff_exp[diff_exp$q_value==min(diff_exp$q_value),]
dim(most_sign)