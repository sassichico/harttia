library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(scales)
library(seqinr)

csv_files <- list.files(pattern = "^Chr\\d+\\.csv$")
slm_data <- map_dfr(csv_files, ~ {
  data <- read_csv(.x, show_col_types = FALSE)
  data$chromosome <- tools::file_path_sans_ext(.x)
  data <- data %>% filter(grepl("^SLMs_", Name))
  return(data)
})

slm_data <- slm_data %>%
  mutate(Position = as.numeric(gsub(",", "", Minimum)))

fasta_seqs <- read.fasta("Pterygoplichthys.fasta", seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
chr_lengths_df <- data.frame(
  full_name = names(fasta_seqs),
  length = getLength(fasta_seqs)
)

chr_mapping <- setNames(paste0("Chr", 1:26), chr_lengths_df$full_name[1:26])
chr_lengths_df$chromosome <- chr_mapping[chr_lengths_df$full_name]
chr_lengths_df <- chr_lengths_df %>% filter(!is.na(chromosome))

chrom_order <- paste0("Chr", 1:26)
chr_lengths_df$chromosome <- factor(chr_lengths_df$chromosome, levels = chrom_order)
slm_data$chromosome <- factor(slm_data$chromosome, levels = chrom_order)

ideogram_data <- chr_lengths_df %>% rename(end = length) %>% mutate(start = 0)

final_plot <- ggplot() +
  geom_segment(data = ideogram_data,
               aes(x = start, xend = end, y = chromosome, yend = chromosome),
               color = "grey80", size = 5, lineend = "round") +
  geom_point(data = slm_data,
             aes(x = Position, y = chromosome, color = Type),
             alpha = 0.7, size = 1.5) +
  scale_x_continuous(labels = comma) +
  labs(title = "Genomic Distribution of Sex-Linked Markers",
       x = "Genomic Position (bp)", y = "Chromosome", color = "Feature Type") +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

print(final_plot)
ggsave("SLM_Distribution.png", width = 10, height = 8, dpi = 300)




