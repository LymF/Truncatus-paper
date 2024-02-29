if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

if (!requireNamespace("patchwork", quietly = TRUE))+
  install.packages("patchwork")

library(patchwork)

library(DESeq2)
library(ggplot2)


count_data <- read.csv("C:\\Users\\Lucas Melo\\Desktop\\count_data.csv", header = TRUE, row.names = 1, sep = ",")


rounded_count_data <- round(count_data)
transposed_count_data <-t(rounded_count_data)

sample_info <- read.csv("C:\\Users\\Lucas Melo\\Desktop\\sample_info.csv", header = TRUE, row.names = 1)


dds <- DESeqDataSetFromMatrix(countData = transposed_count_data, colData = sample_info, design = ~Condition)

dds <- DESeq(dds)



results <- results(dds)

#Change the condition levels
results_contrast <- results(dds, contrast = c("Condition", "Iws", "un"))
results_contrast
 #Save results_contrast to a TSV file
write.table(results_contrast, "results_ABMNEG_ABMPOS.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

#summary(results_contrast)

p_value_cutoff <- 0.05


viruses <- c("TtDV-1", "AVT", "BVT-1", "BVT-2", "TtOV", "TtTV", "CVT", "TtNoV", "GVT-1", "ANT", "PVY", "CVA", "TtDV-2", "GVT 2", "GVT 3", "TtUV", "GVT 4", "eve1", "eve2", "eve3")


volcano_data <- data.frame(
  Log2FC = results_contrast$log2FoldChange,
  P_value = results_contrast$pvalue,
  Highlight = ifelse(rownames(results_contrast) %in% viruses, "Virus", "Normal")
)


virus_labels <- data.frame(
  Log2FC = volcano_data$Log2FC[volcano_data$Highlight == "Virus"],
  P_value = volcano_data$P_value[volcano_data$Highlight == "Virus"],
  VirusName = as.character(viruses[volcano_data$Highlight == "Virus"])
)


volcano_data <- data.frame(
  Log2FC = results_contrast$log2FoldChange,
  P_value = results_contrast$pvalue,
  Highlight = ifelse(rownames(results_contrast) %in% viruses, "Virus", "Normal")
)

p_value_cutoff <- 0.05

virus_labels_filtered <- subset(virus_labels, Log2FC != 0)


volcano_data <- data.frame(
  Log2FC = results_contrast$log2FoldChange,
  P_value = results_contrast$pvalue,
  Highlight = ifelse(rownames(results_contrast) %in% viruses, "Virus", "Normal")
)

p_value_cutoff <- 0.05

virus_labels_filtered <- subset(virus_labels, Log2FC != 0)


volcano_data$DEStatus <- ifelse(
  volcano_data$Highlight == "Normal" &
    ((volcano_data$Log2FC >= 2 * log(2) & volcano_data$P_value <= p_value_cutoff)),
  "Upregulated",
  ifelse(
    volcano_data$Highlight == "Normal" &
      ((volcano_data$Log2FC <= -2 * log(2) & volcano_data$P_value <= p_value_cutoff)),
    "Downregulated",
    "Not DE"
  )
)

volcano_data$PointSize <- ifelse(volcano_data$Highlight == "Virus", "Virus", volcano_data$DEStatus)




# Defina as coordenadas para o zoom
x_zoom <- c(-5, 5)
y_zoom <- c(0, 10)

# Crie o gráfico principal (volcano plot)
volcano_plot <- ggplot(data = volcano_data, aes(x = Log2FC, y = -log10(P_value))) +
  geom_point(aes(color = PointSize, size = PointSize), alpha = 0.7) +
  geom_text(data = virus_labels_filtered, aes(x = Log2FC, y = -log10(P_value), label = VirusName), hjust = 0, vjust = 0) +
  scale_color_manual(values = c("Virus" = "red", "Not DE" = "grey", "Upregulated" = "blue", "Downregulated" = "green")) +
  scale_size_manual(values = c("Virus" = 4, "Not DE" = 2, "Upregulated" = 2, "Downregulated" = 2)) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-log10(p-value)") +
  geom_hline(yintercept = -log10(p_value_cutoff), linetype = "dashed", color = "red") +
  ggtitle("Volcano Plot of DESeq2 Results") +
  labs(
    subtitle = paste("Contrast: ", "HighTemp vs OrdTemp")
  ) +
  theme(legend.position="none")  # Remover a legenda para simplificar

# Crie o gráfico de zoom
zoom_plot <- ggplot(data = volcano_data, aes(x = Log2FC, y = -log10(P_value))) +
  geom_point(aes(color = PointSize, size = PointSize), alpha = 0.7) +
  geom_text(data = virus_labels_filtered, aes(x = Log2FC, y = -log10(P_value), label = VirusName), hjust = 0, vjust = 0) +
  scale_color_manual(values = c("Virus" = "red", "Not DE" = "grey", "Upregulated" = "blue", "Downregulated" = "green")) +
  scale_size_manual(values = c("Virus" = 4, "Not DE" = 2, "Upregulated" = 2, "Downregulated" = 2)) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-log10(p-value)") +
  ggtitle("Zoom Inset") +
  xlim(x_zoom) + ylim(y_zoom)  # Defina os limites para o zoom

# Combine os dois gráficos usando patchwork
final_plot <- volcano_plot / zoom_plot + plot_layout(guides = "collect")

# Exiba o gráfico final
print(final_plot)
print(zoom_plot)
