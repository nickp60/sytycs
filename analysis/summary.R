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
  left_join(key %>% rename(centroid_label=label, centroid_acc=acc, centroid_species=species, centroid_strain=strain) %>% select(centroid_label, centroid_acc, centroid_species, centroid_strain)) %>% 
  group_by(centroid_label) %>% 
  mutate(
    same_accession = acc==centroid_acc,
    centroid_present_in_multiple_accessions=!all(same_accession),
    centroid_present_in_multiple_species=!all(species==centroid_species),
  ) %>% 
  ungroup()
dat %>%
  group_by(acc, species) %>% 
  summarize(length_diff = max(size) - min(size)) %>%
  ggplot(aes(y=species, x=length_diff)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_point() + labs("Within-genome ASV length differences")

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
species_ambiguity <- dat %>%
  mutate(species= ifelse(endsWith(suffix = "_sp.", species), strain, species)) %>% 
  mutate(centroid_species= ifelse(endsWith(suffix = "_sp.",  centroid_species), centroid_strain, centroid_species)) %>% 
  filter(centroid_label != "*") %>% select(species, centroid_species ) %>% distinct() %>% filter(species != centroid_species) 
plot(igraph::graph_from_data_frame(species_ambiguity, directed = FALSE))

species_ambiguity %>% mutate(tmp=label, label=centroid_label, centroid_label=tmp) %>% select(-tmp)

tre <- treeio::read.newick("sytycs_Enterococcus-AYTGGGYDTAAAGNG-CCGTCAATTYHTTTRAGT/alignment/amplicons.mafft.aln.nwk")
tre$tip.label <- gsub("-RC", "", gsub("(.*?)--.*", "\\1", tre$tip.label))

library(ggtree)
enterococcus_cols <-  c(
  "avium"="#8DD3C7", casseliflavus="yellow", cecorum="#BEBADA",dispar= "#FB8072",durans="#80B1D3",faecalis= "#FF7F00", 
  faecium="forestgreen", gallinarum="firebrick",gilvus= "blue",hirae= "#BC80BD",innesii= "#CCEBC5",lactis= "#FFED6F", 
  montenedrensis="#A6CEE3", mundtii="#1F78B4",raffinosus= "#B2DF8A",saigonensis= "purple","sp."= "#FB9A99",thailandicus= "#6A3D9A", 
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
)
names(enterococcus_cols) <- paste0("Enterococcus_", names(enterococcus_cols))
enterococcus_shapes=1:19
names(enterococcus_shapes) <- names(enterococcus_cols)[1:19]

ggtree::ggtree(tre) %<+% (key %>% select(-label)) +
  ggtree::geom_tiplab(size=2, aes(color=species), hjust = -.1) + 
  ggtree::geom_tippoint(size=2, aes(shape=species, color=species)) + 
  scale_color_manual(values =enterococcus_cols) + 
  scale_shape_manual(values=enterococcus_shapes)

ggplot2::ggsave(width = 15, height = 200, filename = "tmp.pdf", limitsize = FALSE)
