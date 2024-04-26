set -xeu
set -o pipefail
indir=$1
cat $indir/pcr/*.fasta | sed "s| |--|g" > $indir/cluster_input.fasta
vsearch --cluster_fast $indir/cluster_input.fasta --threads 16  --id 1 --strand both  -uc $indir/clusters.uc
