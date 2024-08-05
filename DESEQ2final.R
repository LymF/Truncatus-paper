if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

install.packages("ggplot2")


library("DESeq2")

# import counts table per gene - TRANSCRIPTOME
counts <- read.table("C:/Users/Administrator/Desktop/Ortogolos/countsdrosophilagenefinal.tsv", 
                     header = TRUE, row.names = 1, sep = "\t")

sample_info <- read.csv("C:/Users/Administrator/Desktop/ArquivosR/sample_infowolbachia.csv",
                        header = TRUE, row.names = 1)


sample_info$Condition <- factor(sample_info$Condition)
summary(sample_info)

#checking if data is annotated correctly 
all(rownames(sample_info) == colnames(counts)) # must be TRUE
all(rownames(sample_info) %in% colnames(counts)) # must be TRUE



### CREATING DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = sample_info,
                              design = ~ Condition)

dds

# using relevel, specify the REFERENCE condition:
dds$Condition <- relevel(dds$Condition, ref = "un")


########################################################
#### create DESeq2 objects andestimate size factors ####
########################################################

dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

resultsall <- results(dds)

res <- results(dds)
res

write.table(res, "C:/Users/Administrator/Desktop/ArquivosR/resultsIs.tsv", sep = "\t", 
            quote = FALSE, row.names = TRUE)

res <- results(dds, contrast=c("Condition","Iws","un"))

res
write.table(res, "C:/Users/Administrator/Desktop/ArquivosR/results_Is_un.tsv", sep = "\t", 
            quote = FALSE, row.names = TRUE)


# p-values and adjusted p-values
resOrdered <- res[order(res$padj),]
resultado_filtrado <- as.data.frame(resOrdered)

resultado_filtrado[resultado_filtrado$padj < 0.01,]

na.omit(resultado_filtrado[resultado_filtrado$padj < 0.01
                           & resultado_filtrado$log2FoldChange < -2 , ])


# We can summarize some basic tallies using the summary function.
summary(res)

# How many adjusted p-values were less than 0.01?
sum(res$padj < 0.01, na.rm=TRUE)

# Filtering only genes with p-value less than 0.01?
res001 <- results(dds, alpha=0.01)
summary(res001)

sum(res001$padj < 0.05, na.rm=TRUE)

res_diff_expr <- subset(res, padj < 0.05 & !is.na(padj))
write.table(res_diff_expr, "C:/Users/Administrator/Desktop/ArquivosR/differential_genescoinfection.tsv", sep = "\t", 
            quote = FALSE, row.names = TRUE)

library("EnhancedVolcano")

#https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

EnhancedVolcano(res_diff_expr,
                lab = rownames(res_diff_expr),
                title = "ABMPOS vs ABMNEG",
                x = 'log2FoldChange',
                y = 'pvalue')





if (!requireNamespace("plyr", quietly = TRUE))
  install.packages("plyr")
library(plyr)

# Get results for all contrasts
resultsall <- results(dds)


# Extract the unique levels of the conditions in the results
unique_conditions <- unique(resultsall$Condition)

# Create a vector of colors corresponding to the unique levels
colors <- c("Is" = "blue", "Iw" = "red", "Iws" = "green", "un" = "black")

# Map the conditions in the results to the corresponding colors
resultsall$color <- colors[as.character(resultsall$Condition)]

# Plot Volcano plot with differentiated colors
EnhancedVolcano(resultsall,
                lab = rownames(resultsall),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Volcano Plot - All Contrasts",
                labSize = 3,
                xlim = c(-10, 10),
                ylim = c(0, 50),
                pCutoff = 0.05,
                FCcutoff = 2,
                legendLabels = names(colors),
                legendLabSize = 10,
                legendPos = "top")

