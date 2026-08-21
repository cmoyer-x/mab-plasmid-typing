library(ggplot2)
library(dplyr)
library(readr)
library(forcats)

BG <- "#fdf4ff"
DARK <- "#4a044e"
MID <- "#9d4edd"
ACCENT <- "#d4537e"
LIGHT <- "#c9a0dc"

types_data <- read_tsv("plasmid_final_types.tsv", show_col_types = FALSE)

plot_data <- types_data %>%
  filter(source != "excluded") %>%
  mutate(
    is_singleton = startsWith(final_type_label, "singleton_p"),
    category = case_when(
      source %in% c("paper_cluster_anchor", "paper_cluster_nearest", "paper_cluster_ambiguous") ~ "Paper cluster (2021)",
      TRUE ~ "Novel cluster"
    ),
    display_label = if_else(is_singleton, "Singletons (combined)", final_type_label)
  ) %>%
  count(display_label, name = "n_members") %>%
  arrange(desc(n_members))

category_lookup <- types_data %>%
  filter(source != "excluded", !startsWith(final_type_label, "singleton_p")) %>%
  mutate(
    category = if_else(
      source %in% c("paper_cluster_anchor", "paper_cluster_nearest", "paper_cluster_ambiguous"),
      "Paper cluster (2021)", "Novel cluster"
    )
  ) %>%
  distinct(final_type_label, category) %>%
  rename(display_label = final_type_label)

singleton_row <- tibble::tibble(display_label = "Singletons (combined)", category = "Mixed singletons")
category_lookup <- bind_rows(category_lookup, singleton_row)

plot_data <- plot_data %>% left_join(category_lookup, by = "display_label")

plot_data$display_label <- factor(plot_data$display_label, levels = rev(plot_data$display_label))

p <- ggplot(plot_data, aes(x = display_label, y = n_members, fill = category)) +
  geom_col(color = BG, linewidth = 0.4) +
  geom_text(aes(label = n_members), hjust = -0.3, color = DARK, size = 3.2) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Paper cluster (2021)" = DARK,
      "Novel cluster" = ACCENT,
      "Mixed singletons" = LIGHT
    ),
    name = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL,
    y = "Number of genomes",
    title = "Plasmid cluster sizes across the cohort",
    subtitle = "Paper clusters (Dedrick et al. 2021) vs novel cohort-specific clusters"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = BG, color = NA),
    panel.background = element_rect(fill = BG, color = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = LIGHT, linewidth = 0.2),
    axis.text = element_text(color = DARK),
    axis.title = element_text(color = DARK),
    plot.title = element_text(color = DARK, face = "bold", size = 14),
    plot.subtitle = element_text(color = MID, size = 9),
    legend.position = "top",
    legend.text = element_text(color = DARK)
  )

ggsave("plasmid_cluster_sizes.png", p, width = 8, height = 9, dpi = 300, bg = BG)
ggsave("plasmid_cluster_sizes.pdf", p, width = 8, height = 9, bg = BG)

print(p)
