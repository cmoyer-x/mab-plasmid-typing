library(ggplot2)
library(dplyr)
library(readr)
library(forcats)

BG <- "#fdf4ff"
DARK <- "#4a044e"
MID <- "#9d4edd"
ACCENT <- "#d4537e"
LIGHT <- "#c9a0dc"

matrix_data <- read_tsv("plasmid_dcc_matrix.tsv", show_col_types = FALSE)
assoc_data <- read_tsv("plasmid_dcc_association.tsv", show_col_types = FALSE)

sig_lookup <- assoc_data %>%
  select(plasmid_type, most_associated_dcc, direction, fisher_p_fdr) %>%
  rename(dcc = most_associated_dcc)

plot_data <- matrix_data %>%
  left_join(sig_lookup, by = c("plasmid_type", "dcc")) %>%
  mutate(
    sig_label = case_when(
      !is.na(fisher_p_fdr) & fisher_p_fdr < 0.001 ~ "***",
      !is.na(fisher_p_fdr) & fisher_p_fdr < 0.01 ~ "**",
      !is.na(fisher_p_fdr) & fisher_p_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

type_order <- assoc_data %>%
  arrange(fisher_p_fdr) %>%
  pull(plasmid_type)

plot_data$plasmid_type <- factor(plot_data$plasmid_type, levels = rev(type_order))

p <- ggplot(plot_data, aes(x = dcc, y = plasmid_type, fill = frac_carriers_of_type_in_dcc)) +
  geom_tile(color = BG, linewidth = 0.6) +
  geom_text(aes(label = sig_label), color = "white", size = 4, fontface = "bold", vjust = 0.7) +
  scale_fill_gradient(
    low = LIGHT, high = DARK,
    name = "% of type's\ncarriers in DCC",
    labels = scales::percent
  ) +
  labs(
    x = "Clonal Complex (DCC)",
    y = "Plasmid Type",
    title = "Plasmid type distribution across clonal complexes",
    subtitle = "* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001 (Fisher's exact, vs most-associated DCC)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = BG, color = NA),
    panel.background = element_rect(fill = BG, color = NA),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = DARK),
    axis.text.y = element_text(color = DARK, face = "bold"),
    axis.title = element_text(color = DARK),
    plot.title = element_text(color = DARK, face = "bold", size = 14),
    plot.subtitle = element_text(color = MID, size = 9),
    legend.title = element_text(color = DARK, size = 9),
    legend.text = element_text(color = DARK)
  )

ggsave("plasmid_dcc_heatmap.png", p, width = 9, height = 6, dpi = 300, bg = BG)
ggsave("plasmid_dcc_heatmap.pdf", p, width = 9, height = 6, bg = BG)

print(p)
