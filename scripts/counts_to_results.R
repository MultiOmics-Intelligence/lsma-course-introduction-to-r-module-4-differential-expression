# Module 4 project · Counts to results
# One script that reads the files, runs every step and saves the outputs, so you can rerun it
# from the top. Run it in the project folder (open the folder as an RStudio Project),
# so that "data/..." and "results/..." point to the right places.
# Predict before you run: write down which genes you expect to change, and in which direction.

library(DESeq2)
library(ggplot2)
library(ggrepel)
dir.create("results", showWarnings = FALSE)


# ---- Task 1 · Counts (lesson 1) ----
# TODO: read data/practice_counts.txt with read.delim(..., comment.char = "#"),
#       using the gene ID column as row names.
# TODO: keep the four count columns (the last four, after the six gene-information columns)
#       as a matrix, and rename them to match the sample names in data/practice_samples.csv,
#       in the same order.
# TODO: check the size with dim(), and the library size of each sample with colSums().
# TODO: save results/library_sizes.csv with the columns sample, library_size


# ---- Task 2 · Differential expression (lessons 1 and 2) ----
# TODO: read data/practice_samples.csv and build the sample table (colData) with the sample
#       names as row names. Check its row order matches the count columns.
# TODO: make condition a factor with control as the reference level (relevel(), as for dex in lesson 1).
# TODO: build the dataset with DESeqDataSetFromMatrix(..., design = ~ condition), run DESeq(),
#       and extract the treated-versus-control results with results().
# TODO: save the full results table to results/deseq2_results.csv, one row per gene,
#       with a gene column first:
#       res_df <- as.data.frame(res)
#       write.csv(data.frame(gene = rownames(res_df), res_df), "results/deseq2_results.csv", row.names = FALSE)
# TODO: list the genes with padj < 0.05 (exclude NA padj values explicitly) and save them to
#       results/significant_genes.csv with the columns gene, log2FoldChange, padj


# ---- Task 3 · Compare with the truth (lesson 2) ----
# The true changes: PRG003 x4, PRG004 x3, PRG010 x2 and PRG007 x1/4.
# TODO: turn DESeq2's log2 fold changes back into plain ratios with 2^log2FoldChange.
# TODO: compute CPM (each count / the sample's library size x 1,000,000), average the two control
#       and the two treated samples for each gene, and divide treated by control.
# TODO: save results/truth_comparison.csv with one row per built-in gene and the columns
#       gene, truth, deseq2, cpm      (all three as plain ratios, for example 4 and 0.25)
# Then ask why the estimates differ (look at sizeFactors(dds)).


# ---- Task 4 · A volcano plot and a methods paragraph (lesson 3, and Module 3) ----
# TODO: draw a volcano plot of your results, labelling only the genes with padj < 0.05
#       (geom_text_repel(data = ...)). Check each label sits on the side its fold change says.
# TODO: save it as results/volcano.png with ggsave().
# TODO: write results/methods.md by hand: see the README for what it should contain.
#       packageVersion("DESeq2") and R.version.string give the versions.
