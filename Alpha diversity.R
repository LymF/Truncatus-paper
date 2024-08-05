# Load libraries
library(ggplot2)
library(dplyr)
library(rstatix)

# Read data
dados <- read.csv("C:/Users/Administrator/Desktop/truncatusquant.tsv", sep = ",")

# Define colors
cores_condicao <- c("A" = "lightblue", "B" = "salmon", "C" = "lightyellow")

# Calculate p-values for comparisons between conditions

wilcox_results <- dados %>%
  group_by(Superior) %>%
  summarize(
    p_value_AB = ifelse(sum(Condition == "A") & sum(Condition == "B"), 
                        wilcox.test(TPM[Condition == "A"], TPM[Condition == "B"], alternative = "two.sided")$p.value, NA),
    p_value_AC = ifelse(sum(Condition == "A") & sum(Condition == "C"), 
                        wilcox.test(TPM[Condition == "A"], TPM[Condition == "C"], alternative = "two.sided")$p.value, NA),
    p_value_BC = ifelse(sum(Condition == "B") & sum(Condition == "C"), 
                        wilcox.test(TPM[Condition == "B"], TPM[Condition == "C"], alternative = "two.sided")$p.value, NA)
  )
write.table(wilcox_results, "C:/Users/Administrator/Desktop/lucas/wilcoxtwosidedalphadiversity.tsv", sep = "\t")

# Join p-values with original data
dados <- merge(dados, wilcox_results, by = "Superior")

# Plot with significance labels
ggplot(dados, aes(x = Condition, y = TPM, fill = Condition)) +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  scale_fill_manual(values = cores_condicao) +
  geom_boxplot(position = position_dodge(width = 0.1), color = "black") +
  geom_jitter(color = 1) +
  geom_line(aes(group = Virus), position = position_dodge(width = 0.1), color = "black") +
  geom_text(aes(label = ifelse(!is.na(p_value_AB), paste0("AB p-value: ", round(p_value_AB, 3)), "")),
            x = 1.5, y = max(dados$TPM), hjust = 0, vjust = -1, size = 4) +
  geom_text(aes(label = ifelse(!is.na(p_value_AC), paste0("AC p-value: ", round(p_value_AC, 3)), "")),
            x = 2.5, y = max(dados$TPM), hjust = 0, vjust = -1, size = 4) +
  geom_text(aes(label = ifelse(!is.na(p_value_BC), paste0("BC p-value: ", round(p_value_BC, 3)), "")),
            x = 3.5, y = max(dados$TPM), hjust = 0, vjust = -1, size = 4) +
  labs(x = "Condition", y = "TPM") +
  facet_wrap(~Superior, ncol = 4) +
  theme(legend.position = "bottom", panel.spacing = unit(0.1, "lines"))
