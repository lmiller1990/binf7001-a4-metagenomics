data <- data.frame(
  Sample_ID = c("Bran_1", "Bran_2", "Bran_3", "CO_1", "CO_2", "CO_3", "PS_1", "PS_2", "PS_3"),
  Host_Reads = c(2938, 3000, 2800, 3200, 3100, 2900, 3050, 2950, 2850),
  Microbial_Reads = c(7062, 7000, 7200, 6800, 6900, 7100, 6950, 7050, 7150)
)

# Calculate the ratio of microbial to host reads
data$Ratio = data$Microbial_Reads / data$Host_Reads

# Print the updated data frame
print(data)

# Optionally, create a plot of the ratios
library(ggplot2)

plot <- ggplot(data, aes(x = Sample_ID, y = Ratio, fill = Sample_ID)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(title = "Ratio of Microbial to Host Reads", x = "Sample ID", y = "Ratio") +
  theme_minimal()

ggsave("ratios_plot.png", plot, width = 10, height = 8, dpi = 300)
