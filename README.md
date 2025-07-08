# So you think you can (trust your amplicon) species (annotations)?

# Prereqs

- clone this repo
- set environment variable SYTYCS_CACHEDIR to some location so that genomes can be cached between runs: `export SYTYCS_CACHEDIR=$HOME/sytycs_cache`
- a conda env with emboss so we can use primersearch, eg `mamba create -n sytycs emboss mafft fasttree vsearch biopython bioconductor-biostrings bioconductor-ggtree -c conda-forge libiconv -c  r-tidyverse`

## Usage:

Download the desired genus's genomes and run primer search specifying the maximum expected amplicon size :

```
bash run.sh Enterococcus AYTGGGYDTAAAGNG CCGTCAATTYHTTTRAGT 400 
```

Then run the clustering and tree building.  You can add any additional sequences (eg from your own ASVs).

> the header for your additional sequences must be <UNIQUE_ID>:1-<sequence length>--<Genus>--<species>--<strain> . Anything after a comma will be ignored.  this format allows for nice consistent plots


```
bash cluster.sh sytycs_Enterococcus-AYTGGGYDTAAAGNG-CCGTCAATTYHTTTRAGT
```

Lastly, summarize the results:
```
Rscript analysis/summary.R  sytycs_Enterococcus-AYTGGGYDTAAAGNG-CCGTCAATTYHTTTRAGT
```

## V1 Goals:

- add some ANI-based approach to pre-screen for poor homology that would ruin alignments
- Swap primersearch for bbmap  msa + bbmap extract
- implement ability to run on higher taxonomic ranks
- implement subsets
- switch to infernal for alignment
- arg to provide outgroup
