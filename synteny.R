library(dplyr)
library(stringr)
library(ggtree)
library(ggplot2)
library(cowplot)
library(purrr)

setwd("/storage/zhenyingLab/kourongao/Gmel_synteny")
# blocks <- read.delim("LOC113519197_local_block.tsv")
blocks <- read.delim("all_species_final_blocks.tsv")
target_gene <- "LOC113519197"
target_row <- grep(target_gene, blocks$Galleria_mellonella)
context <- 5
ref_sp <- "Galleria_mellonella"
sub_blocks <- blocks[seq(target_row - context, target_row + context),]
gene_location <- read.delim("all_species.bed", header = F)
colnames(gene_location) <- c("Genome", "Start", "End", "Gene", "Score", "Strand")
gene_location$Species <- paste0(str_split_fixed(gene_location$Genome, '_', 3)[, 1], '_', str_split_fixed(gene_location$Genome, '_', 3)[, 2])
gene_location <- subset(gene_location, Gene %in% unlist(as.vector(sub_blocks)))
gene_location[duplicated(gene_location$Gene),] # 708 items, 684 genes, all duplications are tRNA genes

species_list <- c(ref_sp, "Achroia_grisella", "Plodia_interpunctella", "Bombyx_mori", "Danaus_plexippus", "Plutella_xylostella")
species_num <- length(species_list)

sub_blocks <- sub_blocks[species_list]

reverse <- function(x, start, end, strand){
  val <- c(x[start], x[end])
  if (x[strand] == "-") {
    val <- rev(val)
  }
  x[start] <- as.numeric(val[1])
  x[end] <- as.numeric(val[2])
  return(x)
}

pair_orthologs <- data.frame()
for (i in seq_along(species_list[-1])) {
  sp_a <- species_list[i]
  sp_b <- species_list[i+1]
  sp_pair <- sub_blocks[c(sp_a, sp_b)]
  colnames(sp_pair) <- c("sp_A", "sp_B")
  sp_pair <- sp_pair %>% 
    subset(sp_A != '.' & sp_B != '.') %>%
    left_join(gene_location %>% select(-Genome, -Score) %>% rename(sp_A = Gene), by = 'sp_A') %>%
    rename(Start_A = Start, End_A = End, Strand_A = Strand, Species_A = Species) %>%
    left_join(gene_location %>% select(-Genome, -Score) %>% rename(sp_B = Gene), by = 'sp_B') %>% 
    rename(Start_B = Start, End_B = End, Strand_B = Strand, Species_B = Species) 
  pair_orthologs <- rbind(pair_orthologs, sp_pair)
}


par(oma=c(0,0,0,0),mar=c(0,0,0,0))

width <- 0.9
sv_width <- 1.6

left_margin <- min(gene_location$Start)
right_margin <- max(gene_location$End)

group_stats <- gene_location %>%
  group_by(Species) %>%
  summarise(
    initial_offset = 0,
    aligned_min = min(Start) - initial_offset,
    aligned_max = max(End) - initial_offset,
  ) %>%
  mutate(
    centre_offset = (right_margin - left_margin)/2 - aligned_min - (aligned_max - aligned_min)/2
  ) %>%
  mutate(
    total_shift = centre_offset - initial_offset,
    aligned_min = aligned_min + total_shift,
    aligned_max = aligned_max + total_shift
  )

left_min_list <- setNames(as.list(group_stats$aligned_min), group_stats$Species)
right_max_list <- setNames(as.list(group_stats$aligned_max), group_stats$Species)
mid_offset_list <- setNames(as.list(group_stats$centre_offset), group_stats$Species)

gene_location <- gene_location %>%
  left_join(group_stats %>% select(Species, total_shift), by = "Species") %>%
  mutate(
    Start = Start + total_shift,
    End = End + total_shift
  ) %>%
  group_by(Species) %>%
  mutate(overlap = map2_lgl(Start, End, function(s, e) {
    any(Start < s & End > e)
  })
  ) %>%
  ungroup()

## important to update margin value here!
left_margin <- min(gene_location$Start) # this code must before reverse()
right_margin <- max(gene_location$End)
inter <- (right_margin - left_margin)/30

gene_location <- as.data.frame(t(apply(gene_location, 1, reverse, start = "Start", end = "End", strand = "Strand")))
target_gene_ortho <- blocks[target_row,][grep('gene-', blocks[target_row,])]
gene_location$color <- ifelse(gene_location$Gene %in% target_gene_ortho, yes = "red", no = '#84C9F7')

shift_map <- group_stats %>% select(Species, total_shift) %>%
  rename(shift_A = total_shift)

pair_orthologs <- pair_orthologs %>%
  left_join(shift_map %>% rename(Species_A = Species, shift_A = shift_A), by = 'Species_A') %>%
  left_join(shift_map %>% rename(Species_B = Species, shift_B = shift_A), by = 'Species_B') %>%
  mutate(
    Start_A = Start_A + shift_A,
    End_A = End_A + shift_A,
    Start_B = Start_B + shift_B,
    End_B = End_B + shift_B
  ) %>%
  select(-shift_A, -shift_B)

pair_orthologs <- t(apply(pair_orthologs, 1, reverse, start = "Start_A", end = "End_A", strand = "Strand_A"))
pair_orthologs <- as.data.frame(t(apply(pair_orthologs, 1, reverse, start = "Start_B", end = "End_B", strand = "Strand_B")))

draw_gene <- function(x, start, end, track, width, color='#84C9F7', border="NA"){
  if (color %in% names(x)) {
    color = x[color]
  }
  start <- as.numeric(x[start])
  end <- as.numeric(x[end])
  # start and end are absolute coordinates in genome, don't care about strand
  # if start > end, it automatically indicates the gene locates on neg strand
  polygon(c(start, start, start + 0.85*(end - start), end, start + 0.85*(end - start), start), 
          c(track*10+width, track*10-width, track*10-width, track*10, track*10+width, track*10+width),
          col=color, border=border)
}

draw_bezier <- function(x, start_a_col, end_a_col, track, start_b_col, end_b_col, width, t=seq(0, 1, length=100), synteny_col=rgb(220,220,220,max=255), highlight = F) {
  a_y <- track * 10
  b_y <- (track - 1) * 10
  start_a <- as.numeric(x[start_a_col])
  start_b <- as.numeric(x[start_b_col])
  end_a <- as.numeric(x[end_a_col])
  end_b <- as.numeric(x[end_b_col])
  wid_a <- a_y-width-0.2
  wid_b <- b_y+width+0.2
  p1 <- matrix(c(start_a, wid_a, 
                 start_a + (start_b - start_a) / 4, (wid_a + wid_b) / 2,
                 start_a + 3 * (start_b - start_a)/4, (wid_a + wid_b) / 2,
                 start_b, wid_b),
               nrow=4, ncol=2, byrow=TRUE)
  p2 <- matrix(c(start_a + 0.85 * (end_a - start_a), 
                 wid_a,
                 start_a + 0.85 * (end_a - start_a) + (start_b + 0.85 * (end_b - start_b) - (start_a + 0.85 * (end_a - start_a))) / 4,
                 (wid_a + wid_b) / 2,  
                 start_a + 0.85 * (end_a - start_a) + 3 * (start_b + 0.85 * (end_b - start_b) - (start_a + 0.85 * (end_a - start_a))) / 4,
                 (wid_a + wid_b) / 2, 
                 start_b + 0.85 * (end_b - start_b),
                 wid_b),
               nrow=4, ncol=2, byrow=TRUE)
  x1 <- bezier(t=t, p=p1)[,1]
  y1 <- bezier(t=t, p=p1)[,2]
  x2 <- rev(bezier(t=t, p=p2)[,1])
  y2 <- rev(bezier(t=t, p=p2)[,2])
  cx <- c(x1, x2, start_a)
  cy <- c(y1, y2, wid_a)
  
  if (is.character(highlight)) {
    synteny_col = x[highlight]
    }
  polygon(cx, cy, col=synteny_col, border=synteny_col)
}


#### plot axis and genes
# pdf("temp.pdf", width = 12, height = 10)
blank <- 150000
plot(0, 0, type='n', axes=F, main='', xlab='', ylab='', 
     xlim=c(left_margin - inter - blank, right_margin + inter + blank), 
     ylim=c(1*10-5,species_num*10+5)) # ylim: min-5=1*10-5, max+5=(species_num*10)+5

for (i in 1:species_num){
  sp <- species_list[i]
  left <- left_min_list[[sp]]
  right <- right_max_list[[sp]]
  y_track <- species_num - i + 1
  segments(x0 = left - inter, 
           y0 = y_track * 10, 
           x1 = right + inter, 
           y1 = y_track * 10, lwd=4, lend=1, col=rgb(128,128,128,max=255))
  
  text(x = left - inter * 2, 
       y = y_track * 10,
       labels = sp,
       adj = c(1, 0.5), 
       cex = 0.8)
  sub <- subset(gene_location, Species == sp)
  apply(sub, 1, draw_gene, start = 'Start', end = 'End', track = y_track, width = width, color = 'color')
  dup <- subset(sub, overlap == T) 
  if (nrow(dup) > 0) {
    dup$overlap <- ifelse(dup$Gene %in% target_gene_ortho, yes = "red", no = '#84A9A9')
    apply(dup, 1, draw_gene, start = 'Start', end = 'End', track = y_track, width = width, color = 'overlap')
  }
}

#### plot synteny relationship curved polygon
grey_synteny_col = rgb(220,220,220,max=255)
t <- seq(0, 1, length=100)
for (i in 2:species_num){
  sp <- species_list[i]
  y_track <- species_num - i + 2
  sub <- subset(pair_orthologs, Species_B == sp)
  sub$highlight <- ifelse(sub$sp_B %in% target_gene_ortho, "darkblue", 
                          ifelse(sub$sp_B %in% subset(gene_location, overlap == TRUE)$Gene, "#9D9D9D", "#DCDCDC"))
  apply(sub, 1, draw_bezier, start_a = 'Start_A', end_a = 'End_A', start_b = 'Start_B', end_b = 'End_B', track = y_track, width = width, t=t, synteny_col=grey_synteny_col, highlight = 'highlight')
}
for (i in seq_along(target_gene_ortho)) {
  corr <- as.list(subset(gene_location, Gene %in% target_gene_ortho)[i, ])
  midcorr <- (as.numeric(corr[2]) + as.numeric(corr[3])) / 2
  gene <- gsub('gene-', '', corr[4])
  print(gene)
  text(midcorr + 100000, (species_num - grep(corr[7], species_list) + 1) * 10 + 2, labels = gene, adj = c(1, 0.5), cex = 0.8)
}

# legend rule
segments(left, 4, left + 50000, 4)
text(x = 35000 + left, y = 6, labels = "50kb", adj = c(1, 0.5), cex = 0.8)

# dev.off()

tree_text <- paste0(
  "((((((((Achroia_grisella:0.1,Galleria_mellonella:0.1):0.1,Plodia_interpunctella:0.2):0.1),Bombyx_mori:0.3):0.1),Danaus_plexippus:0.4):0.1),Plutella_xylostella:0.5);"
)

tree <- read.tree(text = tree_text)

tree_plot <- ggtree(tree, ladderize = TRUE) + 
  geom_tiplab(size = 4, offset = 0.02, fontface = "italic") +
  geom_nodelab(size = 3, color = "darkred", hjust = -0.1) + 
  xlim(0, 10) +
  theme_tree()
ggsave("species_tree.pdf", tree_plot, height = 10, width = 8)
synteny_img <- ggdraw() + 
  draw_image("temp.pdf")

combined_plot <- plot_grid(
  tree_plot, 
  synteny_img,
  ncol = 2,
  rel_widths = c(0.3, 0.7),
  align = "h"
)

