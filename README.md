# So you think you can (trust your) species (amplicons)?

# Prereqs

- clone this repo
- set environment variable SYTYCS_CACHEDIR to some location so that genomes can be cached between runs: `export SYTYCS_CACHEDIR=$HOME/sytycs_cache`
- a conda env with emboss so we can use primersearch, eg `mamba create -n sytycs emboss mafft fasttree vsearch biopython -c conda-forge libiconv`

## Usage:

Download the desired genus's genomes and run primer search:

```
bash run.sh Enterococcus AYTGGGYDTAAAGNG CCGTCAATTYHTTTRAGT
```

Then run the clustering and tree building

```
bash cluster.sh sytycs_Enterococcus-AYTGGGYDTAAAGNG-CCGTCAATTYHTTTRAGT
```

## V1 Goals:

- add some ANI-based approach to pre-screen for poor homology that would ruin alignments
- Swap primersearch for bbmap  msa + bbmap extract
- implement ability to run on higher taxonomic ranks
- implement subsets
- switch to infernal for alignment
- arg to provide outgroup
