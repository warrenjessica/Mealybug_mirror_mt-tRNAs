#R script to filter the mitochondrial YAMAT reads and then plot their abundance in CPM across the length of the mitogenome
#Script below is for the first replicate, and was repeated for the other two replicates with idential code

library(stringr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(scales) 

rep1<-read.delim("Mito_1_CKDL240032380-1A_22GKFMLT4_L2_collapsed.ID_above1.CCA_blast_above1_CCA_fmt6.txt", header = F, sep = "\t")
colnames(rep1) <- c("query_ID", "subject_strand", "query_length", "percent_identity", "alignment_length", 
                    "query_start", "query_end", "subject_start", "subject_end", 
                    "e_value", "bit_score", "query_coverage")
# Extract ID (number after '_collapsed_')
rep1$ID <- as.integer(sub(".*_collapsed_(\\d+)-.*", "\\1", rep1$query_ID))
# Extract read count (number after last hyphen)
rep1$read_count <- as.integer(sub(".*-(\\d+)$", "\\1", rep1$query_ID))
rep1$replicate <- "Mito_1"

hsp_filtered_rep1 <- subset(rep1, query_coverage > 95 & percent_identity > 95)
unique(hsp_filtered_rep1$ID[duplicated(hsp_filtered_rep1$ID)])
mito1<-hsp_filtered_rep1[c(15,13,14,2,8,9)]

aggregated_rep1 <- aggregate(read_count ~ subject_start + subject_strand, data = mito1, FUN = sum)
total_reads <- sum(aggregated_rep1$read_count)
aggregated_rep1$cpm <- (aggregated_rep1$read_count / total_reads) * 1e6

aggregated_rep1$cpm_signed <- ifelse(aggregated_rep1$subject_strand == "minus",
                                     -aggregated_rep1$cpm,
                                     aggregated_rep1$cpm)

# Mitogenome coding annotation bounds
annots <- tribble(
  ~start, ~end,  ~gene,          ~color,
  0, 1511,  "cox1",       "blue",
  1510, 1574,  "trnL-TTA", "green",
  1571, 2224,  "cox2",     "blue",
  2237, 2290,  "trnK-TTT", "green",
  2289, 2351,  "trnD-GTC", "green",
  2466, 3078,  "atp6",     "orange",
  3078, 3831,  "cox3",     "blue",
  3839, 3896,  "trnG-TCC", "green",
  3885, 4218,  "nad3",     "purple",
  4236, 4293,  "trnC-GCA", "green",
  4244, 4302,  "trnA-TGC", "green",
  4302, 4361,  "trnR-TCG", "green",
  4370, 4425,  "trnN-GTT", "green",
  4422, 4478,  "trnS2-TCT","green",
  4479, 4537,  "trnE-TTC", "green",
  4525, 4589,  "trnF-GAA", "green",
  4588, 6106,  "nad5",     "purple",
  6103, 6169,  "trnH-GTG", "green",
  6195, 7331,  "nad4",     "purple",
  7478, 7685,  "nad4L",    "purple",
  7723, 7777,  "trnT-TGT", "green",
  7776, 7832,  "trnP-TGG", "green",
  8014, 8299,  "nad6",     "purple",
  8300, 9392,  "cytb",     "pink",
  9414, 9467,  "trnS1-TGA","green",
  9696, 10556, "nad1",     "purple",
  10550, 10611,"trnL1-TAG","green",
  12499,11696, "12S",      "green",
  11690,10614, "12S",      "green",
  2351, 2461,  "atp8",     "orange",
  12588,13380, "Control_region","green",
  13644,13706, "trnM-CAT","green",
  13707,13765, "trnV-TAC","green",
  13707,13760, "trnY-GTA","green",
  13757,13813, "trnQ-TTG","green",
  13807,13872, "trnI-GAT","green",
  14073,14820, "nad2",     "purple",
  14821,14881, "trnW-TCA","green"
)

# Normalize intervals 
annots2 <- annots %>%
  mutate(
    xmin = pmin(start, end),
    xmax = pmax(start, end),
    xmid = (xmin + xmax) / 2
  )

fill_values <- c(
  plus  = "#762a83",
  minus = "#b58dc4",
  gray  = "gray",
  brown = "brown",
  green = "green",
  blue  = "blue",
  orange= "orange",
  purple= "purple",
  pink  = "pink"
)

y_limits <- c(-70000, 170000) 
label_y  <- y_limits[2] * 0.97

p <- ggplot(aggregated_rep1, aes(x = subject_start, y = cpm_signed, fill = subject_strand)) +
  geom_rect(
    data = annots2,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = color),
    inherit.aes = FALSE,
    alpha = 0.12,
    color = NA,
    show.legend = FALSE
  ) +
  geom_col(width = 50, position = "identity") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = fill_values, breaks = c("plus","minus"), name = "Strand") +
  scale_y_continuous(
    limits = y_limits,
    breaks = seq(-170000, 150000, by = 20000),
    labels = scales::label_number()
  ) +
  labs(
    title = "CPM by Genomic Position _ Mito1",
    x = "Subject Start Position",
    y = "Counts Per Million (Signed)"
  ) +
  coord_cartesian(xlim = c(0, 15000)) +
  theme_minimal()

#add gene labels across the top
p <- p +
  geom_text(
    data = annots2,
    aes(x = xmid, y = label_y, label = gene),
    inherit.aes = FALSE,
    angle = 90,
    vjust = 0,
    size = 2.8
  )

p
