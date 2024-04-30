library(tidyverse)
library(gtsummary)
theme_set(theme_bw())
metadata <- read.csv("sytycs_Enterococcus-AYTGGGYDTAAAGNG-CCGTCAATTYHTTTRAGT/tmp.queryhits.tsv", sep="\t", 
                     col.names = c("Name",
                                   "TaxID",
                                   "BioProject Accession",
                                   "BioProject ID",   "Group",   "SubGroup",        "Size (Mb)",       "GC%",     "Replicons",       "WGS",     "Scaffolds",       "Genes",   "Proteins",        "Release Date",    "Modify Date", "Status",   
                                   "Center",  "BioSample Accession",     "Assembly Accession",      "Reference",       "FTP Path",        "Pubmed ID",       
                                   "Strain"))
dat_pre <- read.csv("sytycs_Enterococcus-AYTGGGYDTAAAGNG-CCGTCAATTYHTTTRAGT/clusters.uc", sep="\t",
         col.names=c("type", "cluster_id", "size", "percid", "strand", "Y", "Z", "is_identical", "label", "centroid_label"))   %>% 
  filter(type != "C") %>% 
  mutate(
    label=gsub("-RC", "", label),
    centroid_label=gsub("-RC", "", centroid_label)
    )

(key<- data.frame(label=dat_pre$label) %>% 
  mutate(
    cleanlabel = gsub(":", "-", gsub("(.*?)--.*", "\\1", label)),
    tmpA= gsub("(.*?)--.*", "\\1", label),
    tmpB= gsub("(.*?)--(.*)", "\\2", label),
    acc=gsub("(.*?):([0-9]+?)-([0-9]+)", "\\1", tmpA),
    start=gsub("(.*?):([0-9]+?)-([0-9]+)", "\\2", tmpA),
    stop=gsub("(.*?):([0-9]+?)-([0-9]+)", "\\3", tmpA),
    species=gsub("(.*?)--(.*?)--.*", "\\1_\\2", tmpB),
    strain=gsub("(.*?),--complete--genome", "\\1", tmpB)
  ) %>% 
    select(-tmpA, -tmpB) %>% 
    group_by(acc, species, strain) %>% 
    mutate(id=paste0(cur_group_id(), "--", row_number())) %>% 
    ungroup() %>% 
    select(cleanlabel, everything())
)
dat <- dat_pre %>% left_join(key) %>% 
  left_join(key %>% rename(centroid_label=label, centroid_acc=acc, centroid_species=species) %>% select(centroid_label, centroid_acc, centroid_species)) %>% 
  group_by(centroid_label) %>% 
  mutate(
    same_accession = acc==centroid_acc,
    centroid_present_in_multiple_accessions=!all(same_accession),
    centroid_present_in_multiple_species=!all(species==centroid_species),
  )
dat %>%
  group_by(acc, species) %>% 
  summarize(length_diff = max(size) - min(size)) %>%
  ggplot(aes(y=species, x=length_diff)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_point() + labs("Within-genome ASV length differences")

tbl_summary(dat %>% group_by(species, acc) %>% 
              summarize(`N 16S`=n()) %>% select(-acc) %>% 
              select(species, `N 16S`), by=`N 16S`) 
# the +1 is because grouping by centoid and getting row number doesnt incorportate the center itself
(smry <- c(
  "N Genomes" = n_distinct(dat$acc),
  "N Species" = n_distinct(dat$species),
  "N 16S amplicons" = nrow(dat),
  "N unique amplicons" = dat %>% filter(type=="S") %>% nrow(),
  "Median amplicons per cluster" = dat %>% filter(type=="H") %>%  group_by(centroid_label) %>% summarize(n=n() + 1) %>% pull(n) %>% median(),
  "Mean amplicons per cluster" = dat %>% filter(type=="H") %>%  group_by(centroid_label) %>% summarize(n=n() + 1) %>% pull(n) %>% mean(),
  "Max amplicons per cluster" = dat %>% filter(type=="H") %>%  group_by(centroid_label) %>% summarize(n=n() + 1) %>% pull(n) %>% max(),
  "N cluster centroids matching multiple genomes" = dat %>% filter(type=="H") %>% filter(centroid_present_in_multiple_accessions) %>%  pull(centroid_label) %>% n_distinct(),
  "N cluster centroids spanning multiple species" = dat %>% filter(type=="H") %>% filter(centroid_present_in_multiple_species) %>%  pull(centroid_label) %>% n_distinct(),
  "N amplicons ambiguous at the species level" = dat %>% filter(type=="H") %>% filter(centroid_present_in_multiple_species) %>% pull(label) %>% n_distinct()
))
nrow(dat)
nrow(dat %>% filter(centroid_label != "*"))

tre <- treeio::read.newick("sytycs_Enterococcus-AYTGGGYDTAAAGNG-CCGTCAATTYHTTTRAGT/alignment/amplicons.mafft.aln.nwk")
tre$tip.label <- gsub("(.*?)--.*", "\\1", tre$tip.label)

library(ggtree)
ggtree::ggtree(tre) %<+% (key %>% select(-label)) +
  ggtree::geom_tiplab(size=2, aes(color=species), hjust = -.1) + 
  ggtree::geom_tippoint(size=4, aes(color=species)) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 12, name = "Set3") , RColorBrewer::brewer.pal(n = 12, name = "Paired")))

ggplot2::ggsave(width = 30, height = 60, filename = "tmp.pdf", limitsize = FALSE)
