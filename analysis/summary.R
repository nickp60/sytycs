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
  filter(type != "C") 
(key<- data.frame(label=dat_pre$label) %>% 
  mutate(
    tmpA= gsub("(.*?)--.*", "\\1", label),
    tmpB= gsub("(.*?)--(.*)", "\\2", label),
    acc=gsub("(.*?)_([0-9]+?)_([0-9]+)", "\\1", tmpA),
    start=gsub("(.*?)_([0-9]+?)_([0-9]+)", "\\2", tmpA),
    stop=gsub("(.*?)_([0-9]+?)_([0-9]+)", "\\3", tmpA),
    species=gsub("(.*?)--(.*?)--.*", "\\1_\\2", tmpB),
    strain=gsub("(.*?),--complete--genome", "\\1", tmpB)
  ) %>% 
    select(-tmpA, -tmpB) %>% 
    group_by(acc, species, strain) %>% 
    mutate(id=paste0(cur_group_id(), "--", row_number())) %>% 
    ungroup()
)
dat <- dat_pre %>% left_join(key) %>% 
  left_join(key %>% rename(centroid_label=label, centroid_acc=acc, centroid_species=species) %>% select(centroid_label, centroid_acc, centroid_species)) %>% 
  group_by(centroid_label) %>% 
  mutate(
    cross_acc_centroid=!all(acc==centroid_acc),
    cross_species_centroid=!all(species==centroid_species),
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
smry <- c(
  "N Genomes" = n_distinct(dat$acc),
  "N Species" = n_distinct(dat$species),
  "N 16S amplicons" = nrow(dat),
  "N unique amplicons" = dat %>% filter(type=="S") %>% nrow(),
  "Median amplicons per cluster" = dat %>% filter(type=="H") %>%  group_by(centroid_label) %>% summarize(n=n() + 1) %>% pull(n) %>% max(),
  "Max amplicons per cluster" = dat %>% filter(type=="H") %>%  group_by(centroid_label) %>% summarize(n=n() + 1) %>% pull(n) %>% median(),
  "N amplicons matching species genomes" = dat %>% filter(type=="H") %>% filter(cross_acc_centroid) %>%  pull(centroid_label) %>% n_distinct(),
  "N clusters spanning multiple species" = dat %>% filter(type=="H") %>% filter(cross_species_centroid) %>%  pull(centroid_label) %>% n_distinct()
)
nrow(dat)
nrow(dat %>% filter(centroid_label != "*"))

tre <- treeio::read.newick("sytycs_Enterococcus-AYTGGGYDTAAAGNG-CCGTCAATTYHTTTRAGT/alignment/amplicons.mafft.aln.nwk")
ggtree::ggtree(tre) + 
  ggtree::geom_tiplab(size=2)
ggplot2::ggsave(width = 30, height = 60, filename = "tmp.pdf", limitsize = FALSE)
