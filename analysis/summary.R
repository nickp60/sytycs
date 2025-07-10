library(tidyverse)
library(ggtree)


sytycs_summary <- function(sytycspath){
# sytycspath = "sytycs_Enterococcus-AYTGGGYDTAAAGNG-CCGTCAATTYHTTTRAGT/"
metadata <- read.csv(file.path(sytycspath, "tmp.queryhits.tsv"), sep="\t", 
                     col.names = c("Name",
                                   "TaxID",
                                   "BioProject Accession",
                                   "BioProject ID",   "Group",   "SubGroup",        "Size (Mb)",       "GC%",     "Replicons",       "WGS",     "Scaffolds",       "Genes",   "Proteins",        "Release Date",    "Modify Date", "Status",   
                                   "Center",  "BioSample Accession",     "Assembly Accession",      "Reference",       "FTP Path",        "Pubmed ID",       
                                   "Strain"))
dat_pre <- read.csv(file.path(sytycspath, "clusters.uc"), sep="\t",
                    col.names=c("type", "cluster_id", "size", "percid", "strand", "Y", "Z", "is_identical", "label", "centroid_label"), header = FALSE)   %>% 
  filter(type != "C") %>% 
  mutate(
    label=gsub("-RC", "", label),
    centroid_label=gsub("-RC", "", centroid_label)
  ) 

(key <- dat_pre %>% select(cluster_id, label) %>% 
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
  left_join(key %>% select(centroid_label=label, centroid_acc=acc, centroid_species=species, centroid_strain=strain)) %>% 
  group_by(centroid_label) %>% 
  mutate(
    same_accession = acc==centroid_acc,
    centroid_present_in_multiple_accessions=!all(same_accession),
    centroid_present_in_multiple_species=!all(species==centroid_species),
  ) %>% 
  ungroup()


write.table(dat, file = file.path(sytycspath, "clean_data.tsv"), sep="\t")
write.table(key, file = file.path(sytycspath, "key.tsv"), sep="\t")
(plength <- dat %>%
  group_by(acc, species) %>% 
  summarize(length_diff = max(size) - min(size)) %>%
  ggplot(aes(y=species, x=length_diff)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_point() + labs("Within-genome ASV length differences"))
ggsave(plength, filename = file.path(sytycspath, "length.pdf"))
saveRDS(object = plength, file = file.path(sytycspath, "lengthplot.rds"))

# make square matrix the lazy way.  I'm sure theres an elegant way out there...
distmat <- matrix(nrow = nrow(dat), ncol = nrow(dat), dimnames = list(dat$label,dat$label), data = 1)
for (i in 1:nrow(dat)){
  a = dat$label[i]
  b = dat$centroid_label[i]
  if(b != "*"){
    # populate above and below diagonals
    distmat[a, b] = 0
    distmat[b, a] = 0
  }
}

distmat = as.dist(distmat)

saveRDS(object = distmat, file = file.path(sytycspath, "distance_matrix.rds"))

# tbl_summary(dat %>% filter(!grepl("asv", label)) %>% 
#               group_by(species, acc) %>% 
#               summarize(`N 16S`=n()) %>% select(-acc) %>% 
#               select(species, `N 16S`), by=`N 16S`) 
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
saveRDS(object = smry, file = file.path(sytycspath, "summary.rds"))


species_ambiguity <- dat %>%
  mutate(species= ifelse(endsWith(suffix = "_sp.", species), strain, species)) %>% 
  mutate(centroid_species= ifelse(endsWith(suffix = "_sp.",  centroid_species), centroid_strain, centroid_species)) %>% 
  filter(centroid_label != "*") %>% select(species, centroid_species ) %>% distinct() %>% filter(species != centroid_species) 
saveRDS(object = species_ambiguity, file = file.path(sytycspath, "species_ambiguity.rds"))

pdf(file = file.path(sytycspath, "sequence_clusters.pdf"), width = 14, height = 14)
plot(igraph::graph_from_data_frame(species_ambiguity, directed = FALSE))
dev.off()


tre <- treeio::read.newick(file.path(sytycspath, "alignment/amplicons.mafft.aln.nwk"))
tre$tip.label <- gsub("-RC", "", gsub("(.*?)--.*", "\\1", tre$tip.label))


ggtree::ggtree(tre) %<+% (key %>% select(-label)) +
  ggtree::geom_tiplab(size=2, aes(color=species), hjust = -.1) + 
#  ggtree::geom_tippoint(size=2, aes(shape=species, color=species))
  ggtree::geom_tippoint(size=2, aes(color=species))

ggplot2::ggsave(width = 15, height = 200, filename = file.path(sytycspath, "tree_full.pdf"), limitsize = FALSE)





representatives = dat %>% group_by(cluster_id, species)  %>% slice(1) %>%
  group_by(cluster_id) %>% filter(n() > 1) %>% pull(cleanlabel) 

table(representatives %in% tre$tip.label)
representatives[!representatives  %in%tre$tip.label]



ggtree::ggtree(treeio::drop.tip(tre, tre$tip.label[!tre$tip.label %in% representatives] )) %<+% (key %>% select(-label)) +
  ggtree::geom_tiplab(size=3, aes(color=as.character(cluster_id), label=species), hjust = -.1) + 
#  ggtree::geom_tippoint(size=3, aes(shape=species, color=as.character(cluster_id))) + 
  ggtree::geom_tippoint(size=3, aes(color=as.character(cluster_id))) + 
  guides(color="none") + coord_cartesian(clip="off")
ggplot2::ggsave(width = 15, height = 20, filename = file.path(sytycspath, "tree_clusters.pdf"), limitsize = FALSE)

}

if (!interactive()){
  theme_set(theme_bw())
  args=commandArgs(trailingOnly = TRUE)
  sytycspath = args[1]
  sytycs_summary(sytycspath = sytycspath)
}



