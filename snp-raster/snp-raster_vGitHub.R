# Load libraries
library(readxl)
library(dplyr)
library(ggplot2)

# Read Excel file
df <- read_excel("snp-sites-deidentified-final.xlsx")

# Clean numeric columns
df$segment  <- as.numeric(df$segment)
df$length   <- as.numeric(df$length)

df$position <- gsub(",", "", df$position)
df$position <- trimws(df$position)
df$position[df$position == "n/a"] <- NA
df$position <- as.numeric(df$position)

df <- df %>% filter(!is.na(position))

# Build segment reference
segment_lengths <- df %>%
  group_by(segment) %>%
  summarize(length = max(length), .groups = "drop") %>%
  arrange(segment) %>%
  mutate(
    seg_start = lag(cumsum(length), default = 0),
    seg_end   = seg_start + length
  )

# Concatenated genome coordinates
df2 <- df %>%
  left_join(segment_lengths, by = "segment") %>%
  mutate(genome_pos = seg_start + position)

# Update SNP type
df2 <- df2 %>%
  # Convert SNP to proper case to match colors
  mutate(snp_type = case_when(
    SNP == "synonymous"    ~ "Synonymous",
    SNP == "nonsynonymous" ~ "Nonsynonymous",
    TRUE                    ~ "Other"
  ))

# Order samples
df2 <- df2 %>% arrange(serotype, sample_ID)
df2$sample_ID <- factor(df2$sample_ID, levels = unique(df2$sample_ID))

# Filter out "Other" (i.e., "n/a")
df2 <- df2 %>% filter(snp_type != "Other")

# Row boundaries
row_bounds <- df2 %>%
  distinct(sample_ID) %>%
  mutate(y = as.numeric(sample_ID),
         ymin = y - 0.5,
         ymax = y + 0.5)

# Make sure sample_ID is a factor in the correct order
df2$sample_ID <- factor(df2$sample_ID, levels = unique(df2$sample_ID))

# Base plot with y-axis labels
p <- ggplot() +
  
  # SNP vertical bars
  geom_linerange(
    data = df2,
    aes(
      x = genome_pos,
      ymin = as.numeric(sample_ID) - 0.45,
      ymax = as.numeric(sample_ID) + 0.45,
      color = snp_type
    ),
    linewidth = 0.6
  ) +
  
  # Sample row borders
  geom_hline(data = row_bounds, aes(yintercept = ymin),
             color = "grey80", linewidth = 0.3) +
  geom_hline(data = row_bounds, aes(yintercept = ymax),
             color = "grey80", linewidth = 0.3) +
  
  # Segment boundaries
  geom_vline(data = segment_lengths, aes(xintercept = seg_start),
             linetype = "dashed", color = "grey60") +
  geom_vline(data = segment_lengths, aes(xintercept = seg_end),
             linetype = "dashed", color = "grey60") +
  
  # SNP colors
  scale_color_manual(values = c(
    "Synonymous" = "#1a80bb",
    "Nonsynonymous" = "#ea801c"
  )) +
  
  # Remove axis labels
  labs(x = NULL, y = NULL, color = NULL) +
  
  theme_classic() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 15),
    axis.ticks.y = element_line(),
    legend.position = "none"
  ) +
  
  # Explicit x-axis limits
  scale_x_continuous(
    limits = c(0, max(segment_lengths$seg_end)),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = 1:length(unique(df2$sample_ID)),
    labels = levels(df2$sample_ID),
    expand = c(0, 0)
  )

p
