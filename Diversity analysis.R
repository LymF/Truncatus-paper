library(vegan)

library(ggplot2)

# Read the data
dados <- read.table("C:\\Users\\Lucas Melo\\Desktop\\quanttablewolbachia.tsv", header = TRUE, sep = "\t")

# Round the data
dadoround <- round(dados)
transposeddadoround <- t(dados)
sampleinfo <- read.table("C:\\Users\\Lucas Melo\\Desktop\\sampleinfowolbachia.tsv", header = TRUE, row.names = 1, sep = "\t")

resultados_diversidade <- data.frame(biblioteca = colnames(dadoround), diversity = NA)

# Loop to calculate alpha diversity for each virus under each condition
for (i in 1:ncol(dadoround)) {
  diversity_alpha <- diversity(dadoround[, i, drop = FALSE])
  ##Print the output, add to excel and make another table to plot the boxplot!!
  print(diversity_alpha)
}

# Remove rows with NA values
resultados_diversidade <- na.omit(resultados_diversidade)

# Create a boxplot with points colored by virus
ggplot(sampleinfo, aes(x = sampleinfo$condicaosuperior, y = Diversity)) +
  geom_boxplot(fill = "lightblue", color = "blue") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 2) +  # Add jittered points
  labs(title = "Alpha Diversity by Condition", x = "Condition", y = "Alpha Diversity")

# Identify rows containing only zeros
non_zero_rows <- apply(transposeddadoround != 0, 1, any)

# Remove rows containing only zeros
data_without_zeros <- transposeddadoround[non_zero_rows, ]

# Calculate the distance matrix between samples
distances <- vegdist(data_without_zeros, method = "bray")
# Remove NA values
distances_no_na <- na.omit(distances)

# Perform an NMDS analysis
nmds <- metaMDS(distances)

# Add the NMDS results to your information table
sampleinfo$NMDS1 <- nmds$points[, 1]
sampleinfo$NMDS2 <- nmds$points[, 2]

# Create an NMDS scatter plot
ggplot(sampleinfo, aes(x = NMDS1, y = NMDS2, color = condicaosuperior)) +
  geom_point(size=4) +
  labs(title = "Beta Diversity Analysis (NMDS)", x = "NMDS1", y = "NMDS2", color = "Condition")




